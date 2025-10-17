#!/usr/bin/env python3
"""
Batch analysis of gravitational wave events for Rotor Field Theory sideband testing
"""
import os
import sys
from gwpy.timeseries import TimeSeries
from gwosc.datasets import event_gps
import subprocess

# High-priority events (χ_eff > 0.25) from GWTC-2.1 and GWTC-3
HIGH_PRIORITY = [
    ('GW200208_222617', 0.45),  # Highest χ_eff in GWTC-3
    ('GW190805_211137', 0.37),  # Second highest
    ('GW200306_093714', 0.32),  # Third
    ('GW200322_091133', 0.27),  # Fourth
]

# Medium-priority events (0.15 < χ_eff < 0.25)
MEDIUM_PRIORITY = [
    ('GW190916_200658', 0.20),
    ('GW190930_133541', 0.19),
    ('GW200316_215756', 0.13),
    ('GW200308_173609', 0.16),
    ('GW190728_064510', 0.13),
]

def analyze_event(event_name, chi_eff):
    """Download and analyze a single GW event"""
    print(f"\n{'='*70}")
    print(f"ANALYZING {event_name} (χ_eff = {chi_eff})")
    print('='*70)

    try:
        # Get GPS time
        gps = event_gps(event_name)
        print(f"GPS time: {gps:.1f}")

        # Check if already analyzed
        outdir = f'gw_analysis/outputs/gw_out_{event_name}'
        if os.path.exists(f'{outdir}/analysis_complete.txt'):
            print(f"✓ Already analyzed, skipping...")
            return True

        # Download data
        print("\nDownloading H1 and L1 data...")
        segment = (gps-16, gps+16)

        h1_path = f'data/{event_name}_H1.hdf5'
        l1_path = f'data/{event_name}_L1.hdf5'

        if not os.path.exists(h1_path):
            h1 = TimeSeries.fetch_open_data('H1', *segment, cache=True)
            h1.write(h1_path, overwrite=True)
            print(f"✓ H1: {len(h1)} samples @ {h1.sample_rate.value} Hz")
        else:
            print(f"✓ H1 data already exists")

        if not os.path.exists(l1_path):
            l1 = TimeSeries.fetch_open_data('L1', *segment, cache=True)
            l1.write(l1_path, overwrite=True)
            print(f"✓ L1: {len(l1)} samples @ {l1.sample_rate.value} Hz")
        else:
            print(f"✓ L1 data already exists")

        # Run analysis
        os.makedirs(outdir, exist_ok=True)

        print(f"\nRunning sideband analysis...")
        cmd = [
            'python', 'gw_sidebands.py',
            '--h1', h1_path,
            '--l1', l1_path,
            '--gps', str(gps),
            '--pre', '10',
            '--post', '10',
            '--outdir', outdir
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            # Mark as complete
            with open(f'{outdir}/analysis_complete.txt', 'w') as f:
                f.write(f"Event: {event_name}\n")
                f.write(f"χ_eff: {chi_eff}\n")
                f.write(f"GPS: {gps}\n")
            print(f"\n✓ Analysis complete! Results in {outdir}/")
            return True
        else:
            print(f"\n✗ Analysis failed:")
            print(result.stderr)
            return False

    except Exception as e:
        print(f"\n✗ Error: {e}")
        return False

def main():
    """Run batch analysis"""

    # Determine which events to analyze
    if len(sys.argv) > 1:
        priority = sys.argv[1].lower()
    else:
        priority = 'high'

    if priority == 'all':
        events = HIGH_PRIORITY + MEDIUM_PRIORITY
    elif priority == 'medium':
        events = MEDIUM_PRIORITY
    else:  # high
        events = HIGH_PRIORITY

    print(f"\n{'#'*70}")
    print(f"# BATCH ANALYSIS: {priority.upper()} PRIORITY")
    print(f"# Total events: {len(events)}")
    print(f"{'#'*70}\n")

    successful = 0
    failed = 0

    for event_name, chi_eff in events:
        if analyze_event(event_name, chi_eff):
            successful += 1
        else:
            failed += 1

    print(f"\n{'#'*70}")
    print(f"# BATCH ANALYSIS COMPLETE")
    print(f"# Successful: {successful}/{len(events)}")
    print(f"# Failed: {failed}/{len(events)}")
    print(f"{'#'*70}\n")

if __name__ == '__main__':
    main()
