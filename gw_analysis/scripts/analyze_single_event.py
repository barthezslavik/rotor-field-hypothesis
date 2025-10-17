#!/usr/bin/env python3
"""Quick analysis of single GW event"""
import sys
from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries
import subprocess
import os

event_name = sys.argv[1] if len(sys.argv) > 1 else 'GW190412'

print(f"\n{'='*70}")
print(f"ANALYZING {event_name}")
print('='*70)

try:
    # Get GPS time
    gps = event_gps(event_name)
    print(f"GPS time: {gps:.1f}")
    
    # Download data
    print("\nDownloading H1 and L1 data...")
    segment = (gps-16, gps+16)
    
    h1 = TimeSeries.fetch_open_data('H1', *segment, cache=True)
    h1.write(f'data/{event_name}_H1.hdf5', overwrite=True)
    print(f"✓ H1: {len(h1)} samples @ {h1.sample_rate.value} Hz")
    
    l1 = TimeSeries.fetch_open_data('L1', *segment, cache=True)
    l1.write(f'data/{event_name}_L1.hdf5', overwrite=True)
    print(f"✓ L1: {len(l1)} samples @ {l1.sample_rate.value} Hz")
    
    # Run analysis
    outdir = f'gw_out_{event_name}'
    os.makedirs(outdir, exist_ok=True)
    
    print(f"\nRunning sideband analysis...")
    cmd = f'python gw_sidebands.py --h1 data/{event_name}_H1.hdf5 --l1 data/{event_name}_L1.hdf5 --gps {gps} --pre 10 --post 10 --outdir {outdir}'
    
    result = subprocess.run(cmd, shell=True)
    
    if result.returncode == 0:
        print(f"\n✓ Analysis complete! Results in {outdir}/")
    else:
        print(f"\n✗ Analysis failed")
        
except Exception as e:
    print(f"\n✗ Error: {e}")

