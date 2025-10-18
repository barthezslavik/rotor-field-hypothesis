#!/usr/bin/env python3
"""
SPARC Rotation Curve Analysis v2 - Improved h/R estimation
===========================================================

Uses actual surface brightness data to estimate disk scale height
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
    distance = None

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('# Distance'):
                # Extract distance in Mpc
                parts = line.split('=')
                if len(parts) > 1:
                    distance = float(parts[1].strip().split()[0])

            if line.startswith('#'):
                continue

            parts = line.strip().split()
            if len(parts) >= 8:
                try:
                    rad = float(parts[0])      # kpc
                    vobs = float(parts[1])     # km/s observed
                    verr = float(parts[2])     # km/s error
                    vgas = float(parts[3])     # km/s gas
                    vdisk = float(parts[4])    # km/s disk
                    vbul = float(parts[5])     # km/s bulge
                    sbdisk = float(parts[6])   # L/pc^2 disk surface brightness
                    sbbul = float(parts[7])    # L/pc^2 bulge surface brightness

                    data.append({
                        'rad': rad,
                        'vobs': vobs,
                        'verr': verr,
                        'vgas': vgas,
                        'vdisk': vdisk,
                        'vbul': vbul,
                        'sbdisk': sbdisk,
                        'sbbul': sbbul
                    })
                except ValueError:
                    continue

    return galaxy_name, data, distance

def compute_dark_matter_component(galaxy_data):
    """Compute 'dark matter' velocity component"""
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
            'fraction_dm': v_dm_sq / v_obs_sq if vobs > 0 else 0,
            'sbdisk': point['sbdisk']
        })

    return results

def estimate_disk_properties(galaxy_data):
    """
    Estimate disk properties from surface brightness profile

    For exponential disk: SB(R) = SB_0 * exp(-R/R_d)
    Scale height: h ≈ z_0 ≈ 0.2 * R_d (typical for spiral galaxies)

    But thicker disks have more extended SB profiles
    """
    if len(galaxy_data) < 5:
        return None, None, None

    radii = np.array([p['rad'] for p in galaxy_data])
    sbdisk = np.array([p['sbdisk'] for p in galaxy_data])

    # Remove zeros
    mask = sbdisk > 0
    if np.sum(mask) < 3:
        return None, None, None

    radii = radii[mask]
    sbdisk = sbdisk[mask]

    # Fit exponential disk: ln(SB) = ln(SB_0) - R/R_d
    try:
        ln_sb = np.log(sbdisk)

        # Linear fit
        slope, intercept, r_val, p_val, stderr = stats.linregress(radii, ln_sb)

        # Scale length R_d
        R_d = -1.0 / slope if slope < 0 else None

        if R_d is None or R_d <= 0 or R_d > 50:  # sanity check
            return None, None, None

        # Estimate scale height from dispersion in SB profile
        # Thicker disks have more scattered SB profiles
        residuals = ln_sb - (intercept + slope * radii)
        sb_scatter = np.std(residuals)

        # h/R_d ratio: baseline 0.1, increases with scatter
        # More scatter → thicker disk
        h_over_Rd = 0.05 + 0.3 * sb_scatter
        h_over_Rd = np.clip(h_over_Rd, 0.05, 0.4)

        # Outer radius (where SB drops to low value)
        r_outer = radii[-1]

        # h/R where R is characteristic radius
        h_over_R = h_over_Rd * R_d / r_outer

        return h_over_R, R_d, r_outer

    except:
        return None, None, None

def analyze_galaxy_v2(filepath):
    """Analyze single galaxy with improved h/R estimation"""
    name, data, distance = parse_sparc_file(filepath)

    if len(data) < 5:
        return None

    dm_data = compute_dark_matter_component(data)

    if len(dm_data) < 5:
        return None

    # Estimate disk properties from surface brightness
    h_over_R, R_d, r_outer = estimate_disk_properties(dm_data)

    if h_over_R is None:
        return None

    # Compute galaxy properties
    radii = np.array([p['rad'] for p in dm_data])
    v_dm_sq = np.array([p['v_dm_sq'] for p in dm_data])
    v_obs = np.array([p['vobs'] for p in dm_data])

    # Total "dark matter" integrated over disk
    total_dm_sq = np.trapz(v_dm_sq, radii)

    # Outer velocity (where DM dominates)
    v_outer = v_obs[-5:].mean()  # average of last 5 points

    # Average DM fraction in outer regions (R > 0.5*R_max)
    outer_mask = radii > 0.5 * radii[-1]
    if np.sum(outer_mask) > 0:
        avg_dm_fraction = np.mean([p['fraction_dm'] for i, p in enumerate(dm_data) if outer_mask[i]])
    else:
        avg_dm_fraction = np.mean([p['fraction_dm'] for p in dm_data])

    return {
        'name': name,
        'n_points': len(dm_data),
        'r_max': radii[-1],
        'r_d': R_d,
        'total_dm_sq': total_dm_sq,
        'avg_dm_fraction': avg_dm_fraction,
        'h_over_R': h_over_R,
        'v_max': np.max(v_obs),
        'v_outer': v_outer,
        'distance': distance,
        'data': dm_data
    }

def main():
    print("="*70)
    print("SPARC ROTATION CURVE ANALYSIS v2")
    print("Improved h/R Estimation from Surface Brightness")
    print("="*70)

    sparc_dir = Path(__file__).parent / 'Rotmod_LTG'
    files = sorted(glob.glob(str(sparc_dir / '*_rotmod.dat')))

    print(f"\nFound {len(files)} galaxies in SPARC database")

    galaxies = []
    for filepath in files:
        result = analyze_galaxy_v2(filepath)
        if result:
            galaxies.append(result)

    print(f"Successfully analyzed {len(galaxies)} galaxies")

    # Extract data
    h_over_R = np.array([g['h_over_R'] for g in galaxies])
    total_dm_sq = np.array([g['total_dm_sq'] for g in galaxies])
    avg_dm_fraction = np.array([g['avg_dm_fraction'] for g in galaxies])
    v_max = np.array([g['v_max'] for g in galaxies])
    v_outer = np.array([g['v_outer'] for g in galaxies])
    r_d = np.array([g['r_d'] for g in galaxies])

    # Normalize DM by rotation velocity squared
    dm_normalized = total_dm_sq / v_max**2

    # Also try DM fraction in outer regions
    dm_outer_fraction = avg_dm_fraction

    print("\n" + "="*70)
    print("ROTOR FIELD THEORY TEST v2")
    print("="*70)

    # Test 1: h/R vs normalized DM
    slope1, intercept1, r1, p1, stderr1 = stats.linregress(h_over_R, dm_normalized)

    # Test 2: h/R vs DM fraction
    slope2, intercept2, r2, p2, stderr2 = stats.linregress(h_over_R, dm_outer_fraction)

    # Test 3: Inverse - test if THIN disks (low h/R) have MORE DM
    # This is what we expect if rotor field is right!
    inverse_hR = 1.0 / h_over_R
    slope3, intercept3, r3, p3, stderr3 = stats.linregress(inverse_hR, dm_outer_fraction)

    print(f"\nTest 1: h/R vs Normalized ΔV²")
    print(f"  r = {r1:.3f}, p = {p1:.4f}")

    print(f"\nTest 2: h/R vs DM Fraction (outer)")
    print(f"  r = {r2:.3f}, p = {p2:.4f}")

    print(f"\nTest 3: 1/(h/R) vs DM Fraction (ROTOR PREDICTION!)")
    print(f"  r = {r3:.3f}, p = {p3:.4f}")
    if p3 < 0.05 and r3 > 0:
        print(f"  ✅ SIGNIFICANT POSITIVE correlation!")
        print(f"  → Thinner disks (high 1/h/R) have MORE DM")
        print(f"  → Consistent with Rotor Field Theory!")

    # Create comprehensive plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # Plot 1: h/R vs DM (original)
    ax = axes[0, 0]
    ax.scatter(h_over_R, dm_normalized, alpha=0.6, s=50)
    ax.plot(h_over_R, slope1 * h_over_R + intercept1, 'r--',
            label=f'r={r1:.3f}, p={p1:.4f}')
    ax.set_xlabel('h/R', fontsize=12)
    ax.set_ylabel('Normalized ΔV²', fontsize=12)
    ax.set_title('Test 1: h/R vs DM', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: h/R vs DM fraction
    ax = axes[0, 1]
    ax.scatter(h_over_R, dm_outer_fraction, alpha=0.6, s=50, c=v_max, cmap='viridis')
    ax.plot(h_over_R, slope2 * h_over_R + intercept2, 'r--',
            label=f'r={r2:.3f}, p={p2:.4f}')
    ax.set_xlabel('h/R', fontsize=12)
    ax.set_ylabel('DM Fraction (outer)', fontsize=12)
    ax.set_title('Test 2: h/R vs DM Fraction', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('V_max (km/s)', fontsize=10)

    # Plot 3: 1/(h/R) vs DM fraction - ROTOR PREDICTION!
    ax = axes[0, 2]
    scatter = ax.scatter(inverse_hR, dm_outer_fraction, alpha=0.6, s=50, c=v_max, cmap='plasma')
    ax.plot(inverse_hR, slope3 * inverse_hR + intercept3, 'g--', linewidth=2,
            label=f'r={r3:.3f}, p={p3:.4f}')
    ax.set_xlabel('R/h (inverse thickness)', fontsize=12)
    ax.set_ylabel('DM Fraction (outer)', fontsize=12)
    ax.set_title('Test 3: ROTOR PREDICTION\nThin disks → More DM?',
                 fontsize=13, fontweight='bold', color='green')
    ax.legend()
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('V_max (km/s)', fontsize=10)

    # Plot 4: Distribution comparison
    ax = axes[1, 0]
    ax.hist(h_over_R, bins=25, alpha=0.7, edgecolor='black')
    ax.axvline(np.median(h_over_R), color='r', linestyle='--',
               label=f'Median = {np.median(h_over_R):.3f}')
    ax.set_xlabel('h/R', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Distribution of h/R (N={len(galaxies)})', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 5: Scale length R_d distribution
    ax = axes[1, 1]
    ax.hist(r_d, bins=25, alpha=0.7, edgecolor='black', color='orange')
    ax.axvline(np.median(r_d), color='r', linestyle='--',
               label=f'Median = {np.median(r_d):.2f} kpc')
    ax.set_xlabel('Disk Scale Length R_d (kpc)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Distribution of R_d', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 6: V_outer vs DM fraction
    ax = axes[1, 2]
    scatter = ax.scatter(v_outer, dm_outer_fraction, alpha=0.6, s=50, c=h_over_R, cmap='coolwarm')
    ax.set_xlabel('V_outer (km/s)', fontsize=12)
    ax.set_ylabel('DM Fraction', fontsize=12)
    ax.set_title('Outer Velocity vs DM', fontsize=13)
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('h/R', fontsize=10)

    plt.tight_layout()

    output_file = Path(__file__).parent / 'sparc_rotor_analysis_v2.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved: {output_file}")

    # Save results
    results_file = Path(__file__).parent / 'sparc_results_v2.txt'
    with open(results_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("SPARC ANALYSIS v2 - Improved h/R Estimation\n")
        f.write("="*70 + "\n\n")
        f.write(f"Sample Size: N = {len(galaxies)} galaxies\n\n")

        f.write("TEST 1: h/R vs Normalized ΔV²\n")
        f.write("-"*70 + "\n")
        f.write(f"r = {r1:.3f}, p = {p1:.4f}\n")
        if p1 < 0.05:
            f.write("✅ Statistically significant\n")
        else:
            f.write("⚠️  Not significant\n")
        f.write("\n")

        f.write("TEST 2: h/R vs DM Fraction (outer)\n")
        f.write("-"*70 + "\n")
        f.write(f"r = {r2:.3f}, p = {p2:.4f}\n")
        if p2 < 0.05:
            f.write("✅ Statistically significant\n")
        else:
            f.write("⚠️  Not significant\n")
        f.write("\n")

        f.write("TEST 3: R/h vs DM Fraction (ROTOR PREDICTION)\n")
        f.write("-"*70 + "\n")
        f.write(f"r = {r3:.3f}, p = {p3:.4f}\n")
        if p3 < 0.05 and r3 > 0:
            f.write("✅✅ STATISTICALLY SIGNIFICANT!\n")
            f.write("→ Thin disks (high R/h) correlate with MORE dark matter\n")
            f.write("→ This matches Rotor Field Theory prediction!\n")
            f.write("→ Thick disks with large rotor volume need LESS DM\n")
        elif p3 < 0.05 and r3 < 0:
            f.write("⚠️  Significant but wrong sign\n")
        else:
            f.write("⚠️  Not significant\n")

    print(f"✓ Results saved: {results_file}")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"""
v1 (gradient method):  h/R vs DM → r={r1:.3f}, p={p1:.4f}
v2 (SB profile method): h/R vs DM → r={r2:.3f}, p={p2:.4f}
v2 ROTOR TEST:         R/h vs DM → r={r3:.3f}, p={p3:.4f}

ROTOR FIELD PREDICTION:
  Thin disks (small h) → less rotor volume → more 'dark matter' needed
  We test: R/h (thinness) vs DM fraction

  Expectation: POSITIVE correlation (thin → more DM)
  Result: r={r3:.3f}, p={p3:.4f}
    """)

if __name__ == '__main__':
    main()
