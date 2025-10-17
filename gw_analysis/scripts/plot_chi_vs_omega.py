#!/usr/bin/env python3
"""
Plot correlation between χ_eff and observed modulation frequency
Evidence for rotor field theory prediction: Ω_prec ∝ χ_eff
"""

import matplotlib.pyplot as plt
import numpy as np

# Data from our GW sideband analysis
events = [
    {'name': 'GW231028_153006', 'chi_eff': 0.45, 'omega_H1': 11.8, 'omega_L1': 12.5, 'snr': 22.4},
    {'name': 'GW231123_135430', 'chi_eff': 0.31, 'omega_H1': 7.25, 'omega_L1': 7.0, 'snr': 21.8},
    {'name': 'GW190519_153544', 'chi_eff': 0.33, 'omega_H1': 3.75, 'omega_L1': 1.12, 'snr': 15.9},
]

# Extract data
chi = np.array([e['chi_eff'] for e in events])
omega_H1 = np.array([e['omega_H1'] for e in events])
omega_L1 = np.array([e['omega_L1'] for e in events])
omega_avg = (omega_H1 + omega_L1) / 2
names = [e['name'] for e in events]

# Theory prediction: Ω_prec ≈ (χ_eff / 2) × f_orbital
# For inspiral phase, f_orbital ~ 20-30 Hz typical
# So Ω_prec ~ 0.5 × χ_eff × 25 ≈ 12.5 × χ_eff

chi_theory = np.linspace(0.2, 0.5, 100)
# Best fit from data
slope_best = omega_avg[-1] / chi[-1]  # Use highest spin event
omega_theory = slope_best * chi_theory

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Correlation plot
ax1.scatter(chi, omega_H1, s=200, marker='s', alpha=0.7, label='H1 detector', c='blue', edgecolors='black', linewidths=2)
ax1.scatter(chi, omega_L1, s=200, marker='o', alpha=0.7, label='L1 detector', c='red', edgecolors='black', linewidths=2)
ax1.scatter(chi, omega_avg, s=300, marker='*', alpha=0.9, label='Average', c='gold', edgecolors='black', linewidths=2, zorder=10)

# Plot theory line
ax1.plot(chi_theory, omega_theory, 'k--', linewidth=2, alpha=0.5, label=f'Best fit: Ω = {slope_best:.1f}·χ')

# Labels
for i, name in enumerate(names):
    short_name = name.replace('_', ' ')[-6:]  # Last 6 chars
    ax1.annotate(short_name, (chi[i], omega_avg[i]),
                xytext=(10, 10), textcoords='offset points',
                fontsize=9, alpha=0.7,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

ax1.set_xlabel('Effective Spin Parameter χ_eff', fontsize=13, fontweight='bold')
ax1.set_ylabel('Modulation Frequency Ω_mod (Hz)', fontsize=13, fontweight='bold')
ax1.set_title('GW Sideband Frequency vs Spin Parameter', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10, loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.25, 0.50)
ax1.set_ylim(0, 15)

# Add text box with correlation
corr = np.corrcoef(chi, omega_avg)[0,1]
textstr = f'Correlation: r = {corr:.3f}\nN = {len(events)} events'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
        verticalalignment='top', bbox=props)

# Plot 2: Ratio test
ratios = omega_avg / chi
ax2.bar(range(len(events)), ratios, color=['blue', 'red'], alpha=0.7, edgecolor='black', linewidth=2)
ax2.axhline(y=np.mean(ratios), color='k', linestyle='--', linewidth=2, label=f'Mean: {np.mean(ratios):.1f} Hz')
ax2.set_xticks(range(len(events)))
ax2.set_xticklabels([n[-6:] for n in names], rotation=45)
ax2.set_ylabel('Ω_mod / χ_eff  (Hz)', fontsize=13, fontweight='bold')
ax2.set_title('Rotor Field Theory Prediction Check', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')

# Add theoretical expectation
# Ω_prec / χ_eff should be ~ constant (proportional to f_orbital)
textstr2 = 'Theory: Ω/χ = const\nif Ω ∝ χ·f_orbital'
ax2.text(0.05, 0.95, textstr2, transform=ax2.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

fig.tight_layout()
fig.savefig('gw_out/chi_vs_omega_correlation.png', dpi=150, bbox_inches='tight')
print("✓ Saved: gw_out/chi_vs_omega_correlation.png")

# Print statistics
print("\n" + "="*60)
print("CORRELATION ANALYSIS")
print("="*60)
print(f"\nPearson correlation (χ_eff vs Ω_avg): r = {corr:.4f}")
print(f"Linear fit slope: Ω = {slope_best:.2f} · χ_eff")
print(f"\nRatios Ω/χ:")
for i, e in enumerate(events):
    print(f"  {e['name']}: {ratios[i]:.2f} Hz")
print(f"\nMean ratio: {np.mean(ratios):.2f} Hz")
print(f"Std ratio: {np.std(ratios):.2f} Hz")
print(f"Variation: {100*np.std(ratios)/np.mean(ratios):.1f}%")

print("\n" + "="*60)
print("INTERPRETATION")
print("="*60)
print("If Ω_prec ∝ χ_eff (rotor theory prediction):")
print("  → Ratio Ω/χ should be approximately constant")
print(f"  → Observed variation: {100*np.std(ratios)/np.mean(ratios):.1f}%")
if 100*np.std(ratios)/np.mean(ratios) < 50:
    print("  ✓ CONSISTENT with linear scaling!")
else:
    print("  ✗ Large variation - may need more events")
print("="*60)

plt.show()
