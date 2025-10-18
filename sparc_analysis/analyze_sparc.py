#!/usr/bin/env python3
"""
SPARC Rotation Curve Analysis - Rotor Field Theory Test
=========================================================

Tests the Rotor Field Theory prediction:
    ΔV² = V_obs² - V_Newton² should correlate with disk geometry

Dark Matter interpretation: ΔV² ~ M_DM (dark matter needed)
Rotor Field interpretation: ΔV² ~ f(h/R, j, tidal) (geometry-dependent)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from scipy import stats

def parse_sparc_file(filepath):
    """Parse a SPARC rotation curve file"""
    data = []
    galaxy_name = Path(filepath).stem.replace('_rotmod', '')

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    rad = float(parts[0])      # kpc
                    vobs = float(parts[1])     # km/s observed
                    verr = float(parts[2])     # km/s error
                    vgas = float(parts[3])     # km/s gas
                    vdisk = float(parts[4])    # km/s disk
                    vbul = float(parts[5])     # km/s bulge

                    data.append({
                        'rad': rad,
                        'vobs': vobs,
                        'verr': verr,
                        'vgas': vgas,
                        'vdisk': vdisk,
                        'vbul': vbul
                    })
                except ValueError:
                    continue

    return galaxy_name, data

def compute_dark_matter_component(galaxy_data):
    """
    Compute "dark matter" velocity component

    V_obs² = V_gas² + V_disk² + V_bul² + V_DM²

    Therefore:
    V_DM² = V_obs² - V_gas² - V_disk² - V_bul²
    """
    results = []

    for point in galaxy_data:
        vobs = point['vobs']
        vgas = point['vgas']
        vdisk = point['vdisk']
        vbul = point['vbul']

        # Newtonian velocity (baryons only)
        v_newton_sq = vgas**2 + vdisk**2 + vbul**2
        v_obs_sq = vobs**2

        # "Dark matter" component
        v_dm_sq = v_obs_sq - v_newton_sq

        # Skip negative (unphysical) values
        if v_dm_sq < 0:
            continue

        v_dm = np.sqrt(v_dm_sq)
        v_newton = np.sqrt(v_newton_sq)

        results.append({
            'rad': point['rad'],
            'vobs': vobs,
            'v_newton': v_newton,
            'v_dm': v_dm,
            'v_dm_sq': v_dm_sq,
            'fraction_dm': v_dm_sq / v_obs_sq if vobs > 0 else 0
        })

    return results

def estimate_disk_thickness(galaxy_data):
    """
    Estimate disk scale height from surface brightness profile

    For exponential disks: h/R ~ 0.1-0.2 typically
    We'll use the variance of the rotation curve as a proxy
    """
    if len(galaxy_data) < 3:
        return 0.15  # typical value

    radii = np.array([p['rad'] for p in galaxy_data])
    vobs = np.array([p['vobs'] for p in galaxy_data])

    # Disk scale length (where V starts to flatten)
    r_flat = radii[np.argmax(vobs)]

    # Estimate h/R from velocity gradient
    # Steeper gradient → thinner disk
    if r_flat > 0:
        gradient = np.gradient(vobs, radii)
        avg_gradient = np.mean(np.abs(gradient[:len(gradient)//2]))

        # Normalize: high gradient → low h/R
        h_over_R = 0.05 + 0.2 * np.exp(-avg_gradient / 50)
        return np.clip(h_over_R, 0.05, 0.3)

    return 0.15

def analyze_galaxy(filepath):
    """Analyze single galaxy"""
    name, data = parse_sparc_file(filepath)

    if len(data) < 5:  # need at least 5 points
        return None

    dm_data = compute_dark_matter_component(data)

    if len(dm_data) < 5:
        return None

    # Compute galaxy properties
    radii = np.array([p['rad'] for p in dm_data])
    v_dm = np.array([p['v_dm'] for p in dm_data])
    v_dm_sq = np.array([p['v_dm_sq'] for p in dm_data])
    v_obs = np.array([p['vobs'] for p in dm_data])

    # Total "dark matter" (area under curve)
    total_dm_sq = np.trapz(v_dm_sq, radii)

    # Maximum radius
    r_max = radii[-1]

    # Estimate disk thickness ratio
    h_over_R = estimate_disk_thickness(data)

    # Average DM fraction
    avg_dm_fraction = np.mean([p['fraction_dm'] for p in dm_data])

    return {
        'name': name,
        'n_points': len(dm_data),
        'r_max': r_max,
        'total_dm_sq': total_dm_sq,
        'avg_dm_fraction': avg_dm_fraction,
        'h_over_R': h_over_R,
        'v_max': np.max(v_obs),
        'data': dm_data
    }

def main():
    print("="*70)
    print("SPARC ROTATION CURVE ANALYSIS")
    print("Rotor Field Theory vs Dark Matter")
    print("="*70)

    # Find all SPARC files
    sparc_dir = Path(__file__).parent / 'Rotmod_LTG'
    files = sorted(glob.glob(str(sparc_dir / '*_rotmod.dat')))

    print(f"\nFound {len(files)} galaxies in SPARC database")

    # Analyze all galaxies
    galaxies = []
    for filepath in files:
        result = analyze_galaxy(filepath)
        if result:
            galaxies.append(result)

    print(f"Successfully analyzed {len(galaxies)} galaxies")

    # Extract data for correlation analysis
    h_over_R = np.array([g['h_over_R'] for g in galaxies])
    total_dm_sq = np.array([g['total_dm_sq'] for g in galaxies])
    avg_dm_fraction = np.array([g['avg_dm_fraction'] for g in galaxies])
    v_max = np.array([g['v_max'] for g in galaxies])

    # Normalize for comparison
    dm_normalized = total_dm_sq / v_max**2

    print("\n" + "="*70)
    print("ROTOR FIELD THEORY PREDICTION TEST")
    print("="*70)

    print("\nHypothesis:")
    print("  Dark Matter: ΔV² is independent of disk geometry")
    print("  Rotor Field: ΔV² ∝ f(h/R) - geometry matters!")

    # Statistical test
    slope, intercept, r_value, p_value, std_err = stats.linregress(h_over_R, dm_normalized)

    print(f"\nCorrelation Analysis:")
    print(f"  Pearson r = {r_value:.3f}")
    print(f"  p-value = {p_value:.4f}")
    print(f"  Slope = {slope:.2f} ± {std_err:.2f}")

    if p_value < 0.05:
        print(f"\n✅ STATISTICALLY SIGNIFICANT (p < 0.05)!")
        if r_value < 0:
            print(f"   → Thicker disks (higher h/R) have LESS 'dark matter'")
            print(f"   → Consistent with Rotor Field Theory!")
        else:
            print(f"   → Thicker disks have MORE 'dark matter'")
    else:
        print(f"\n⚠️  Not statistically significant (p > 0.05)")
        print(f"   → Need more analysis or different proxy for h/R")

    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: h/R vs "Dark Matter"
    ax = axes[0, 0]
    ax.scatter(h_over_R, dm_normalized, alpha=0.6, s=50)
    ax.plot(h_over_R, slope * h_over_R + intercept, 'r--',
            label=f'r={r_value:.3f}, p={p_value:.4f}')
    ax.set_xlabel('Disk Thickness Ratio (h/R)', fontsize=12)
    ax.set_ylabel('Normalized "Dark Matter" ΔV²', fontsize=12)
    ax.set_title('Rotor Field Test: Geometry vs "Dark Matter"', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: V_max vs DM fraction
    ax = axes[0, 1]
    ax.scatter(v_max, avg_dm_fraction, alpha=0.6, s=50, c=h_over_R, cmap='viridis')
    ax.set_xlabel('Maximum Velocity V_max (km/s)', fontsize=12)
    ax.set_ylabel('Average DM Fraction', fontsize=12)
    ax.set_title('DM Fraction vs Rotation Speed', fontsize=13)
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('h/R', fontsize=10)

    # Plot 3: Distribution of h/R
    ax = axes[1, 0]
    ax.hist(h_over_R, bins=20, alpha=0.7, edgecolor='black')
    ax.axvline(np.median(h_over_R), color='r', linestyle='--',
               label=f'Median = {np.median(h_over_R):.3f}')
    ax.set_xlabel('Disk Thickness Ratio (h/R)', fontsize=12)
    ax.set_ylabel('Number of Galaxies', fontsize=12)
    ax.set_title('Distribution of h/R (N={})'.format(len(galaxies)), fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Sample rotation curves
    ax = axes[1, 1]
    # Show 5 example galaxies
    for i, gal in enumerate(galaxies[:5]):
        radii = [p['rad'] for p in gal['data']]
        vobs = [p['vobs'] for p in gal['data']]
        ax.plot(radii, vobs, '-o', alpha=0.7, label=gal['name'], markersize=3)

    ax.set_xlabel('Radius (kpc)', fontsize=12)
    ax.set_ylabel('Rotation Velocity (km/s)', fontsize=12)
    ax.set_title('Sample Rotation Curves', fontsize=13)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    output_file = Path(__file__).parent / 'sparc_rotor_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved: {output_file}")

    # Save results
    results_file = Path(__file__).parent / 'sparc_results.txt'
    with open(results_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("SPARC ROTATION CURVE ANALYSIS RESULTS\n")
        f.write("="*70 + "\n\n")
        f.write(f"Sample Size: N = {len(galaxies)} galaxies\n\n")
        f.write("ROTOR FIELD THEORY TEST\n")
        f.write("-"*70 + "\n")
        f.write(f"Correlation (h/R vs ΔV²): r = {r_value:.3f}\n")
        f.write(f"Statistical significance: p = {p_value:.4f}\n")
        f.write(f"Linear fit: ΔV² = {slope:.2f}×(h/R) + {intercept:.2f}\n\n")

        if p_value < 0.05:
            f.write("✅ RESULT: Statistically significant correlation found!\n")
            if r_value < 0:
                f.write("   → Geometry affects 'dark matter' as predicted by Rotor Theory\n")
        else:
            f.write("⚠️  RESULT: No significant correlation (need better h/R proxy)\n")

        f.write("\n" + "="*70 + "\n")
        f.write("GALAXY SAMPLE STATISTICS\n")
        f.write("="*70 + "\n")
        f.write(f"h/R range: {h_over_R.min():.3f} - {h_over_R.max():.3f}\n")
        f.write(f"h/R median: {np.median(h_over_R):.3f}\n")
        f.write(f"V_max range: {v_max.min():.1f} - {v_max.max():.1f} km/s\n")
        f.write(f"Avg DM fraction: {avg_dm_fraction.mean():.2f} ± {avg_dm_fraction.std():.2f}\n")

    print(f"✓ Results saved: {results_file}")

    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    print("""
If Rotor Field Theory is correct:
  • Thicker disks (high h/R) → more 'volume' → stronger rotor field → LESS need for DM
  • We should see NEGATIVE correlation: r < 0

If Dark Matter is correct:
  • h/R should be IRRELEVANT to dark matter amount
  • We should see NO correlation: r ≈ 0
    """)

    print("Analysis complete! Check the plots and results files.")

if __name__ == '__main__':
    main()
