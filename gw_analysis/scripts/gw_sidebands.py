#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gw_sidebands.py
Анализ данных GWOSC (O4a 16kHz): поиск побочных полос (sidebands) у прецессирующих бинаров.

Функции:
- загрузка strain из HDF5 (GWOSC)
- обрезка окна по UTC или GPS
- базовая предобработка (bandpass, notch, whitening)
- поиск временного пика (excess power)
- спектр/спектрограмма вокруг пика
- оценка сайдбэндов через огибающую (Hilbert) узкополосного сигнала
- сравнение H1/L1 по частотам модуляции

Зависимости: gwpy, numpy, scipy, matplotlib, astropy, h5py
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py

from gwpy.timeseries import TimeSeries
from gwpy.signal.filter_design import bandpass, notch
from astropy.time import Time
from scipy.signal import hilbert, butter, sosfiltfilt
from scipy.signal.windows import tukey
from scipy.fft import rfft, rfftfreq

def parse_args():
    p = argparse.ArgumentParser(description="GWOSC sideband scan (H1/L1)")
    p.add_argument("--h1", type=str, default=None, help="H1 HDF5 file (GWOSC)")
    p.add_argument("--l1", type=str, default=None, help="L1 HDF5 file (GWOSC)")
    timegrp = p.add_mutually_exclusive_group(required=True)
    timegrp.add_argument("--utc", type=str, help="UTC time like 2023-05-29T18:15:00.7")
    timegrp.add_argument("--gps", type=float, help="GPS time (seconds)")
    p.add_argument("--pre", type=float, default=200.0, help="seconds before time")
    p.add_argument("--post", type=float, default=200.0, help="seconds after time")
    p.add_argument("--fmin", type=float, default=20.0, help="bandpass low Hz")
    p.add_argument("--fmax", type=float, default=1024.0, help="bandpass high Hz")
    p.add_argument("--line-50hz", action="store_true", help="notch 50/100/150 Hz lines")
    p.add_argument("--line-60hz", action="store_true", help="notch 60/120/180 Hz lines")
    p.add_argument("--carrier-margin", type=float, default=20.0,
                   help="±Hz вокруг несущей для узкополосного фильтра")
    p.add_argument("--outdir", type=str, default="gw_out", help="output folder")
    return p.parse_args()

def utc_to_gps(utc_str: str) -> float:
    # Astropy: by default .gps gives float seconds
    return Time(utc_str, format="isot", scale="utc").gps

def load_segment_from_file(path: str, gps_center: float, pre: float, post: float):
    # GWOSC файлы длинные (4096с); читаем участок, попадающий в файл
    # TimeSeries.read читает весь файл; фильтруем .crop
    ts = TimeSeries.read(path, "strain")
    # если окно выходит за пределы файла — всё равно crop обрежет по границам
    return ts.crop(gps_center - pre, gps_center + post)

def apply_notches(ts: TimeSeries, mains_hz: int):
    # вырез сетевых линий и гармоник до 5-й
    freqs = [mains_hz * k for k in range(1, 6)]
    out = ts.copy()
    for f0 in freqs:
        # узкая вырезка ±0.5 Гц
        out = out.filter(notch(f0, out.sample_rate.value, quality=30))
    return out

def preprocess(ts: TimeSeries, fmin: float, fmax: float,
               notch50: bool=False, notch60: bool=False):
    # Use scipy butter directly for compatibility
    fs = ts.sample_rate.value
    sos = butter(8, [fmin, fmax], btype='band', fs=fs, output='sos')
    filtered_data = sosfiltfilt(sos, ts.value)
    out = TimeSeries(filtered_data, t0=ts.t0, dt=ts.dt)

    if notch50:
        out = apply_notches(out, 50)
    if notch60:
        out = apply_notches(out, 60)
    # whitening: оценим PSD на всём окне и whiten()
    psd = out.psd(fftlength=8, overlap=4)
    w = out.whiten(asd=np.sqrt(psd))
    # лёгкое косинусное окно для краёв
    win = tukey(len(w.value), alpha=0.05)
    # Create new TimeSeries with windowed data
    w2 = TimeSeries(w.value * win, t0=w.t0, dt=w.dt)
    return w2

def find_peak_window(tsw: TimeSeries, span=8.0):
    # найдём максимум энергии через полосатый Q-transform proxy: используем скользящее RMS
    sr = tsw.sample_rate.value
    n = int(span*sr)
    if n < 1: n = 1
    x = tsw.value
    # скользящее среднее мощности
    powv = np.convolve(x*x, np.ones(n)/n, mode='same')
    idx = int(np.argmax(powv))
    t0 = tsw.t0.value + idx/sr
    return (t0 - span/2, t0 + span/2), powv[idx]

def plot_spectrogram(tsw: TimeSeries, tpair, outpng: str):
    sub = tsw.crop(*tpair)
    fig = sub.q_transform(outseg=(tpair[0], tpair[1]), qrange=(8, 64),
                          frange=(20, 1024), whiten=False).plot(figsize=(9,4))
    fig.axes[0].set_title("Q-transform (whitened)")
    fig.savefig(outpng, dpi=140, bbox_inches="tight"); plt.close(fig)

def one_sided_fft(x, fs):
    n = len(x)
    X = rfft(x)
    f = rfftfreq(n, 1.0/fs)
    return f, np.abs(X)

def estimate_carrier_and_sidebands(tsw: TimeSeries, tpair, carrier_margin=20.0, outprefix=""):
    # 1) оценим несущую частоту как пик амплитудного спектра в 30–300 Hz
    sub = tsw.crop(*tpair)
    fs = sub.sample_rate.value
    f, A = one_sided_fft(sub.value, fs)
    mask = (f >= 30) & (f <= 300)
    fpk = f[mask][np.argmax(A[mask])]
    # 2) узкополосный фильтр вокруг несущей
    # делаем sos фильтр Баттерворта 4-го порядка
    flo = max(20.0, fpk - carrier_margin)
    fhi = fpk + carrier_margin
    sos = butter(4, [flo, fhi], btype='band', fs=fs, output='sos')
    x_nb = sosfiltfilt(sos, sub.value)
    # 3) огибающая (Hilbert) и её спектр => частоты модуляции
    env = np.abs(hilbert(x_nb))
    fe, Ae = one_sided_fft(env - np.mean(env), fs)
    # фильтруем до 50 Гц — модуляции обычно низкие
    m = (fe > 0.1) & (fe < 50.0)
    fe2, Ae2 = fe[m], Ae[m]
    # найдём до 3 пиков модуляции
    topk = np.argsort(Ae2)[-3:][::-1]
    mod_peaks = list(zip(fe2[topk], Ae2[topk]))
    # 4) сохраняем диагностические графики
    fig, ax = plt.subplots(1,2, figsize=(11,4))
    ax[0].plot(f, A); ax[0].set_xlim(10, 512); ax[0].set_xlabel("Hz"); ax[0].set_ylabel("|FFT|")
    ax[0].set_title(f"Carrier scan (peak ~ {fpk:.1f} Hz)")
    ax[1].plot(fe2, Ae2); ax[1].set_xlabel("Hz"); ax[1].set_ylabel("|FFT(env)|")
    if len(mod_peaks):
        txt = ", ".join([f"{p[0]:.2f}" for p in mod_peaks])
        ax[1].set_title(f"Envelope spectrum (mods ~ {txt} Hz)")
    else:
        ax[1].set_title("Envelope spectrum")
    fig.suptitle(outprefix)
    fig.tight_layout()
    fig.savefig(f"{outprefix}_carrier_envelope.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    # отдельный график временных форм
    t = np.linspace(tpair[0], tpair[1], len(x_nb))
    fig2, ax2 = plt.subplots(figsize=(10,3))
    ax2.plot(t, x_nb, lw=0.6, alpha=0.8, label="narrowband")
    ax2.set_xlabel("GPS time (s)"); ax2.set_ylabel("strain (whitened)")
    ax2.legend()
    fig2.savefig(f"{outprefix}_narrowband_wave.png", dpi=140, bbox_inches="tight")
    plt.close(fig2)
    return fpk, mod_peaks

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.utc:
        gps_center = utc_to_gps(args.utc)
    else:
        gps_center = float(args.gps)

    print(f"[i] Center time (GPS) = {gps_center:.3f}")

    results = {}
    for det, path in (("H1", args.h1), ("L1", args.l1)):
        if path is None: 
            continue
        print(f"[i] Loading {det} from {os.path.basename(path)} ...")
        # Try to read HDF5 - handle both GWOSC and gwpy formats
        with h5py.File(path, 'r') as f:
            # Check which format this is
            if 'strain/Strain' in f:
                # GWOSC format
                strain_data = f['strain/Strain'][()]
                attrs = dict(f['strain/Strain'].attrs)
                t0_gps = attrs['Xstart']
                sample_rate = attrs['Xspacing']
            else:
                # gwpy format - find the TimeSeries dataset
                # Usually stored as detector name
                keys = list(f.keys())
                if len(keys) > 0:
                    ds = f[keys[0]]
                    strain_data = ds[()]
                    t0_gps = ds.attrs.get('x0', ds.attrs.get('t0', 0))
                    sample_rate = ds.attrs.get('dx', ds.attrs.get('dt', 1/4096.0))
                else:
                    raise ValueError("Cannot find strain data in HDF5 file")

        # Create TimeSeries manually
        ts = TimeSeries(strain_data, t0=t0_gps, dt=sample_rate)
        print(f"[i] {det}: file span = {ts.t0.value:.1f} .. {(ts.t0.value + ts.duration.value):.1f} GPS")

        # окно вокруг события
        seg = ts.crop(gps_center - args.pre, gps_center + args.post)
        if len(seg) == 0:
            print(f"[!] {det}: window outside file range — check file/interval")
            continue
        print(f"[i] {det}: samples = {len(seg)}, fs = {seg.sample_rate.value:.1f} Hz")
        # предобработка
        seg_bp = preprocess(seg, args.fmin, args.fmax,
                            notch50=args.line_50hz, notch60=args.line_60hz)
        # найдём локальный пик энергии и узкое окно вокруг него
        (t0, t1), peakpow = find_peak_window(seg_bp, span=8.0)
        print(f"[i] {det}: peak window {t0:.3f} .. {t1:.3f} (power={peakpow:.3e})")

        # спектрограмма (Q-transform)
        spng = os.path.join(args.outdir, f"{det}_qtransform.png")
        plot_spectrogram(seg_bp, (t0, t1), spng)
        print(f"[i] {det}: saved {spng}")

        # оценка несущей и модуляций (сайдбэнды)
        outprefix = os.path.join(args.outdir, f"{det}")
        fcarrier, mods = estimate_carrier_and_sidebands(
            seg_bp, (t0, t1), carrier_margin=args.carrier_margin, outprefix=outprefix
        )
        print(f"[i] {det}: carrier ~ {fcarrier:.2f} Hz; modulation peaks (Hz): {mods}")
        results[det] = {"carrier": fcarrier, "mods": mods, "peakwin": (t0, t1)}

        # базовый спектр вокруг пика (для отчёта)
        sub = seg_bp.crop(t0, t1)
        f, A = one_sided_fft(sub.value, sub.sample_rate.value)
        fig, ax = plt.subplots(figsize=(8,3))
        ax.plot(f, A); ax.set_xlim(10, 1024)
        ax.set_xlabel("Hz"); ax.set_ylabel("|FFT|")
        ax.set_title(f"{det} spectrum in peak window")
        spfile = os.path.join(args.outdir, f"{det}_spectrum.png")
        fig.savefig(spfile, dpi=140, bbox_inches="tight"); plt.close(fig)
        print(f"[i] {det}: saved {spfile}")

    # сводка по совпадению модчастот H1/L1
    if "H1" in results and "L1" in results:
        def mod_list(x): return np.array([m[0] for m in x]) if x else np.array([])
        hmods = mod_list(results["H1"]["mods"])
        lmods = mod_list(results["L1"]["mods"])
        if hmods.size and lmods.size:
            # наивное сопоставление по ближайшему соседу
            pairs = []
            for hm in hmods:
                idx = np.argmin(np.abs(lmods - hm))
                pairs.append((hm, lmods[idx], float(np.abs(lmods[idx]-hm))))
            print("[i] H1 vs L1 modulation matches (Hz, Hz, |Δ|):")
            for p in pairs:
                print("    ", p)

    print(f"[i] Done. Outputs in: {args.outdir}")

if __name__ == "__main__":
    main()
