#!/usr/bin/env python3
"""
SPARC Extended Rotation Curves Analysis
Test Rotor Field Prediction #6: Extended correlations with angular momentum
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import glob

print("="*70)
print("SPARC EXTENDED ANALYSIS: Angular Momentum Correlation")
print("Testing Rotor Field Prediction #6")
print("="*70)

# Directory with SPARC rotation model data
data_dir = '/Users/slavik/Rotor/data/Rotmod_LTG'

# ========================================================================
# Load SPARC Data
# ========================================================================

def load_galaxy_data(filename):
    """Load SPARC rotation model data file"""
    try:
        with open(filename, 'r') as f:
            # Read distance from header
            first_line = f.readline()
            if 'Distance' in first_line:
                distance = float(first_line.split('=')[1].split('Mpc')[0].strip())
            else:
                distance = None

            # Skip second header line
            f.readline()

            # Read data
            data = []
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        rad = float(parts[0])  # kpc
                        vobs = float(parts[1])  # km/s
                        vgas = float(parts[3])  # km/s
                        vdisk = float(parts[4])  # km/s
                        vbul = float(parts[5])  # km/s
                        sbdisk = float(parts[6])  # L/pc^2

                        data.append({
                            'rad': rad,
                            'vobs': vobs,
                            'vgas': vgas,
                            'vdisk': vdisk,
                            'vbul': vbul,
                            'sbdisk': sbdisk
                        })
                    except ValueError:
                        continue

            return distance, data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None, []

def estimate_disk_properties(galaxy_data):
    """
    Estimate disk thickness h/R from surface brightness profile

    For exponential disk: SB(R) = SB_0 * exp(-R/R_d)
    """
    if len(galaxy_data) < 5:
        return None, None, None

    radii = np.array([p['rad'] for p in galaxy_data])
    sbdisk = np.array([p['sbdisk'] for p in galaxy_data if p['sbdisk'] > 0])
    radii_sb = np.array([p['rad'] for p in galaxy_data if p['sbdisk'] > 0])

    if len(sbdisk) < 5:
        return None, None, None

    # Fit exponential: ln(SB) = ln(SB_0) - R/R_d
    ln_sb = np.log(sbdisk)
    slope, intercept, r_val, p_val, stderr = stats.linregress(radii_sb, ln_sb)

    # R_d = -1/slope (scale length)
    R_d = -1.0 / slope if slope < 0 else None

    if R_d is None or R_d <= 0:
        return None, None, None

    # Estimate h/R from surface brightness scatter
    # Higher scatter → thicker disk
    residuals = ln_sb - (intercept + slope * radii_sb)
    sb_scatter = np.std(residuals)

    # Empirical relation: h/R_d ~ 0.05 + 0.3 * scatter
    h_over_Rd = 0.05 + 0.3 * sb_scatter

    # Outer radius (for angular momentum calculation)
    r_outer = np.max(radii)

    return h_over_Rd, R_d, r_outer

def calculate_angular_momentum(galaxy_data, distance_Mpc):
    """
    Calculate specific angular momentum L_spin = M × V × R

    Using baryonic mass from rotation curve components
    """
    if len(galaxy_data) < 3:
        return None, None

    radii = np.array([p['rad'] for p in galaxy_data])  # kpc
    vobs = np.array([p['vobs'] for p in galaxy_data])  # km/s
    vgas = np.array([p['vgas'] for p in galaxy_data])
    vdisk = np.array([p['vdisk'] for p in galaxy_data])
    vbul = np.array([p['vbul'] for p in galaxy_data])

    # Baryonic mass component: M_baryon ~ V_baryon^2 * R / G
    # V_baryon^2 = V_gas^2 + V_disk^2 + V_bul^2
    v_baryon_sq = vgas**2 + vdisk**2 + vbul**2

    # Total baryonic mass (integrated over disk)
    # Using simple trapezoid integration
    # M_total ~ ∫ (V^2 * R / G) dR
    # In units of (km/s)^2 * kpc ~ 10^10 M_sun
    mass_integrand = v_baryon_sq * radii
    M_baryon = np.trapz(mass_integrand, radii)  # Arbitrary units

    # Specific angular momentum j = L / M
    # j ~ V * R (characteristic)
    # Use value at half-light radius (median radius)
    r_median_idx = len(radii) // 2
    v_char = vobs[r_median_idx]
    r_char = radii[r_median_idx]

    j_spin = v_char * r_char  # km/s * kpc

    # Total angular momentum L = M * j
    L_spin = M_baryon * j_spin  # Arbitrary units: 10^10 M_sun * km/s * kpc

    return L_spin, j_spin

def calculate_dm_component(galaxy_data):
    """
    Calculate "dark matter" component ΔV^2 = V_obs^2 - V_baryon^2

    Averaged over the rotation curve
    """
    if len(galaxy_data) < 3:
        return None

    vobs = np.array([p['vobs'] for p in galaxy_data])
    vgas = np.array([p['vgas'] for p in galaxy_data])
    vdisk = np.array([p['vdisk'] for p in galaxy_data])
    vbul = np.array([p['vbul'] for p in galaxy_data])

    # Baryonic velocity: V_baryon^2 = V_gas^2 + V_disk^2 + V_bul^2
    v_baryon_sq = vgas**2 + vdisk**2 + vbul**2
    v_obs_sq = vobs**2

    # Dark matter component
    delta_v_sq = v_obs_sq - v_baryon_sq

    # Average over curve (excluding negative values which are unphysical)
    delta_v_sq_positive = delta_v_sq[delta_v_sq > 0]

    if len(delta_v_sq_positive) < 3:
        return None

    dm_avg = np.mean(delta_v_sq_positive)

    return dm_avg

# ========================================================================
# Analyze all galaxies
# ========================================================================

print("\nLoading SPARC galaxies...")

galaxy_files = sorted(glob.glob(os.path.join(data_dir, '*_rotmod.dat')))
print(f"Found {len(galaxy_files)} galaxy files")

results = []

for i, filename in enumerate(galaxy_files):
    galaxy_name = os.path.basename(filename).replace('_rotmod.dat', '')

    distance, data = load_galaxy_data(filename)

    if distance is None or len(data) < 5:
        continue

    # Calculate properties
    h_over_R, R_d, r_outer = estimate_disk_properties(data)
    L_spin, j_spin = calculate_angular_momentum(data, distance)
    dm_component = calculate_dm_component(data)

    if h_over_R is not None and L_spin is not None and dm_component is not None:
        results.append({
            'name': galaxy_name,
            'distance': distance,
            'h_over_R': h_over_R,
            'R_d': R_d,
            'r_outer': r_outer,
            'L_spin': L_spin,
            'j_spin': j_spin,
            'dm_component': dm_component
        })

    if (i+1) % 20 == 0:
        print(f"  Processed {i+1}/{len(galaxy_files)} galaxies...")

print(f"\n✓ Successfully analyzed {len(results)} galaxies")

# ========================================================================
# Statistical Analysis
# ========================================================================

print("\n" + "="*70)
print("STATISTICAL CORRELATIONS")
print("="*70)

# Extract arrays
h_over_R_arr = np.array([r['h_over_R'] for r in results])
L_spin_arr = np.array([r['L_spin'] for r in results])
j_spin_arr = np.array([r['j_spin'] for r in results])
dm_arr = np.array([r['dm_component'] for r in results])

# Normalize DM component (for better visualization)
dm_normalized = (dm_arr - np.mean(dm_arr)) / np.std(dm_arr)

# Log transform L_spin for better linear correlation
log_L_spin = np.log10(L_spin_arr)

print(f"\nSample size: N = {len(results)} galaxies")
print(f"\nParameter ranges:")
print(f"  h/R:     {np.min(h_over_R_arr):.3f} - {np.max(h_over_R_arr):.3f}")
print(f"  L_spin:  {np.min(L_spin_arr):.2e} - {np.max(L_spin_arr):.2e}")
print(f"  j_spin:  {np.min(j_spin_arr):.1f} - {np.max(j_spin_arr):.1f} km/s*kpc")
print(f"  ΔV²:     {np.min(dm_arr):.1f} - {np.max(dm_arr):.1f} (km/s)²")

# Test 1: ΔV² vs h/R (repeat from previous analysis)
corr_hR, p_hR = stats.pearsonr(h_over_R_arr, dm_normalized)
print(f"\n1. Correlation (h/R, ΔV²):")
print(f"   r = {corr_hR:.3f}, p = {p_hR:.4f}", "✅ SIGNIFICANT" if p_hR < 0.05 else "")

# Test 2: ΔV² vs log(L_spin) [NEW!]
corr_L, p_L = stats.pearsonr(log_L_spin, dm_normalized)
print(f"\n2. Correlation (log L_spin, ΔV²):")
print(f"   r = {corr_L:.3f}, p = {p_L:.4f}", "✅ SIGNIFICANT" if p_L < 0.05 else "")

# Test 3: ΔV² vs j_spin (specific angular momentum)
corr_j, p_j = stats.pearsonr(j_spin_arr, dm_normalized)
print(f"\n3. Correlation (j_spin, ΔV²):")
print(f"   r = {corr_j:.3f}, p = {p_j:.4f}", "✅ SIGNIFICANT" if p_j < 0.05 else "")

# Test 4: Multivariate regression
# ΔV² ~ β0 + β1*h/R + β2*log(L_spin)
from scipy.stats import linregress

# Prepare data for multivariate analysis
X = np.column_stack([h_over_R_arr, log_L_spin])
y = dm_normalized

# Simple bivariate correlation
from scipy.optimize import curve_fit

def multivariate_model(X, beta0, beta1, beta2):
    h_R, log_L = X
    return beta0 + beta1 * h_R + beta2 * log_L

try:
    popt, pcov = curve_fit(multivariate_model, (h_over_R_arr, log_L_spin), dm_normalized)
    beta0, beta1, beta2 = popt

    # Predicted values
    y_pred = multivariate_model((h_over_R_arr, log_L_spin), *popt)

    # R² statistic
    ss_res = np.sum((dm_normalized - y_pred)**2)
    ss_tot = np.sum((dm_normalized - np.mean(dm_normalized))**2)
    r_squared = 1 - (ss_res / ss_tot)

    print(f"\n4. Multivariate Regression: ΔV² ~ β₀ + β₁*(h/R) + β₂*log(L)")
    print(f"   β₀ = {beta0:.3f}")
    print(f"   β₁ = {beta1:.3f} (h/R coefficient)")
    print(f"   β₂ = {beta2:.3f} (log L coefficient)")
    print(f"   R² = {r_squared:.3f}")

except Exception as e:
    print(f"\n4. Multivariate regression failed: {e}")

# ========================================================================
# Visualization
# ========================================================================

print("\n" + "="*70)
print("CREATING VISUALIZATIONS")
print("="*70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: h/R vs ΔV²
ax = axes[0, 0]
ax.scatter(h_over_R_arr, dm_normalized, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
# Linear fit
slope_hR, intercept_hR = np.polyfit(h_over_R_arr, dm_normalized, 1)
x_fit = np.linspace(np.min(h_over_R_arr), np.max(h_over_R_arr), 100)
ax.plot(x_fit, slope_hR*x_fit + intercept_hR, 'r--', linewidth=2, label=f'Linear fit')
ax.set_xlabel('Disk Thickness h/R', fontsize=12)
ax.set_ylabel('Dark Matter Component ΔV² (normalized)', fontsize=12)
ax.set_title(f'h/R vs ΔV²: r={corr_hR:.3f}, p={p_hR:.4f}', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 2: log(L_spin) vs ΔV²
ax = axes[0, 1]
ax.scatter(log_L_spin, dm_normalized, alpha=0.6, s=50, color='green', edgecolors='black', linewidth=0.5)
slope_L, intercept_L = np.polyfit(log_L_spin, dm_normalized, 1)
x_fit_L = np.linspace(np.min(log_L_spin), np.max(log_L_spin), 100)
ax.plot(x_fit_L, slope_L*x_fit_L + intercept_L, 'r--', linewidth=2, label='Linear fit')
ax.set_xlabel('log₁₀(Angular Momentum L)', fontsize=12)
ax.set_ylabel('Dark Matter Component ΔV² (normalized)', fontsize=12)
ax.set_title(f'log(L_spin) vs ΔV²: r={corr_L:.3f}, p={p_L:.4f}', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 3: j_spin vs ΔV²
ax = axes[1, 0]
ax.scatter(j_spin_arr, dm_normalized, alpha=0.6, s=50, color='purple', edgecolors='black', linewidth=0.5)
slope_j, intercept_j = np.polyfit(j_spin_arr, dm_normalized, 1)
x_fit_j = np.linspace(np.min(j_spin_arr), np.max(j_spin_arr), 100)
ax.plot(x_fit_j, slope_j*x_fit_j + intercept_j, 'r--', linewidth=2, label='Linear fit')
ax.set_xlabel('Specific Angular Momentum j (km/s·kpc)', fontsize=12)
ax.set_ylabel('Dark Matter Component ΔV² (normalized)', fontsize=12)
ax.set_title(f'j_spin vs ΔV²: r={corr_j:.3f}, p={p_j:.4f}', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Plot 4: h/R vs log(L_spin) - check independence
ax = axes[1, 1]
ax.scatter(h_over_R_arr, log_L_spin, alpha=0.6, s=50, color='orange', edgecolors='black', linewidth=0.5)
corr_hL, p_hL = stats.pearsonr(h_over_R_arr, log_L_spin)
slope_hL, intercept_hL = np.polyfit(h_over_R_arr, log_L_spin, 1)
x_fit_hL = np.linspace(np.min(h_over_R_arr), np.max(h_over_R_arr), 100)
ax.plot(x_fit_hL, slope_hL*x_fit_hL + intercept_hL, 'r--', linewidth=2, label='Linear fit')
ax.set_xlabel('Disk Thickness h/R', fontsize=12)
ax.set_ylabel('log₁₀(Angular Momentum L)', fontsize=12)
ax.set_title(f'h/R vs log(L): r={corr_hL:.3f}, p={p_hL:.4f}', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig('sparc_analysis/sparc_extended_analysis.png', dpi=300, bbox_inches='tight')
print("\n✓ Plot saved: sparc_extended_analysis.png")

# ========================================================================
# Summary
# ========================================================================

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("\n🔍 ROTOR FIELD PREDICTION #6: Extended Rotation Curves")
print("\nPrediction: ΔV² should correlate with:")
print("  1. Disk geometry (h/R) - rotor volume effect")
print("  2. Angular momentum (L_spin) - rotor coupling to spin")
print("  3. Multivariate model incorporating both")

print(f"\n📊 RESULTS (N={len(results)} galaxies):")

print(f"\n1. h/R correlation:")
print(f"   r = {corr_hR:.3f}, p = {p_hR:.4f}")
if p_hR < 0.001:
    print(f"   ✅ HIGHLY SIGNIFICANT (p < 0.001)")
elif p_hR < 0.05:
    print(f"   ✅ SIGNIFICANT (p < 0.05)")
else:
    print(f"   ❌ Not significant")

print(f"\n2. log(L_spin) correlation:")
print(f"   r = {corr_L:.3f}, p = {p_L:.4f}")
if p_L < 0.001:
    print(f"   ✅ HIGHLY SIGNIFICANT (p < 0.001)")
elif p_L < 0.05:
    print(f"   ✅ SIGNIFICANT (p < 0.05)")
else:
    print(f"   ❌ Not significant")

print(f"\n3. j_spin correlation:")
print(f"   r = {corr_j:.3f}, p = {p_j:.4f}")
if p_j < 0.001:
    print(f"   ✅ HIGHLY SIGNIFICANT (p < 0.001)")
elif p_j < 0.05:
    print(f"   ✅ SIGNIFICANT (p < 0.05)")
else:
    print(f"   ❌ Not significant")

print("\n" + "="*70)

# Count significant correlations
n_significant = sum([p_hR < 0.05, p_L < 0.05, p_j < 0.05])

if n_significant >= 2:
    print("✅ ROTOR FIELD THEORY: MULTIPLE CORRELATIONS CONFIRMED!")
    print(f"   {n_significant}/3 predicted correlations are statistically significant")
elif n_significant == 1:
    print("⚠️  PARTIAL SUPPORT: One correlation significant")
else:
    print("❌ NO SIGNIFICANT CORRELATIONS FOUND")

print("="*70)

# Save results to file
output_file = 'sparc_analysis/extended_results.txt'
with open(output_file, 'w') as f:
    f.write("SPARC Extended Analysis Results\n")
    f.write("="*70 + "\n\n")
    f.write(f"Sample size: N = {len(results)}\n\n")
    f.write(f"1. h/R correlation: r = {corr_hR:.3f}, p = {p_hR:.4f}\n")
    f.write(f"2. log(L_spin) correlation: r = {corr_L:.3f}, p = {p_L:.4f}\n")
    f.write(f"3. j_spin correlation: r = {corr_j:.3f}, p = {p_j:.4f}\n")
    f.write(f"\nSignificant correlations: {n_significant}/3\n")

print(f"\n✓ Results saved: {output_file}")
