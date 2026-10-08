import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Enforce standard light publication style
plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.edgecolor'] = 'none'
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'

df = pd.read_csv("~/Projects/physics/DESI/clean_project/coma_galaxy_data.csv")

v_median = np.median(df['LOS_V'])
sigma_v = np.std(df['LOS_V'], ddof=1)

print(f"Sample size before filter: {len(df):,}")
print(f"Median Velocity: {v_median} km/s")
print(f"Velocity Dispersion (Sigma): {sigma_v} km/s")

mask_3sigma = np.abs(df['LOS_V'] - v_median) <= (3.0 * sigma_v)

filtered_df = df[mask_3sigma].copy().reset_index(drop=True)

filtered_df['REL_V'] = filtered_df['LOS_V'] - v_median

print(f"Confirmed Coma Members: {len(filtered_df)}")

print(f"Minimum LOS Velocity: {filtered_df['LOS_V'].min()} km/s")
print(f"Maximum LOS Velocity: {filtered_df['LOS_V'].max()} km/s")

bin_edges = np.linspace(0.0, 2.0, 9)
bin_labels = [f"{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}°" for i in range(len(bin_edges) - 1)]

filtered_df['RADIAL_BIN'] = pd.cut(
    filtered_df['ANGLR_SEP_DEG'],
    bins=bin_edges,
    labels=bin_labels,
    include_lowest=True
)

profile_records = []

for label in bin_labels:
    bin_data = filtered_df[filtered_df['RADIAL_BIN'] == label]
    n_gal = len(bin_data)
    
    if n_gal > 1:
        theta_mid = bin_data['ANGLR_SEP_DEG'].median()
        sigma_obs = bin_data['LOS_V'].std(ddof=1)
        mean_v_err_sq = np.mean(bin_data['V_ERR']**2)
        sigma_int = np.sqrt(max(0.0, sigma_obs**2 - mean_v_err_sq))
        sigma_err = sigma_int / np.sqrt(2.0 * (n_gal - 1))
    else:
        theta_mid, sigma_obs, sigma_int, sigma_err = np.nan, np.nan, np.nan, np.nan

    profile_records.append({
        'Bin_Range': label,
        'N_Galaxy': n_gal,
        'Theta_Deg': theta_mid,
        'Sigma_Obs_Kms': sigma_obs,
        'Sigma_Int_Kms': sigma_int,
        'Sigma_Err_Kms': sigma_err
    })

profile_df = pd.DataFrame(profile_records)
print(profile_df.to_string(index=False))

plt.figure(figsize=(8, 5.5))

plt.errorbar(
    profile_df['Theta_Deg'],
    profile_df['Sigma_Int_Kms'],
    yerr=profile_df['Sigma_Err_Kms'],
    fmt='o-',
    color='black',
    ecolor='crimson',
    elinewidth=1.5,
    capsize=4,
    label=r'DESI Coma Core Members ($\sigma_{\mathrm{int}}$)'
)

plt.xlabel(r'Projected Radius $\theta$ [degrees]', fontsize=12)
plt.ylabel(r'Velocity Dispersion $\sigma_v$ [$\mathrm{km\ s^{-1}}$]', fontsize=12)
plt.title('Coma Cluster Velocity Dispersion Profile', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
# plt.show()

theta = profile_df['Theta_Deg'].to_numpy()
sigma = profile_df['Sigma_Int_Kms'].to_numpy()
err = profile_df['Sigma_Err_Kms'].to_numpy()

# Isothermal Parameters
weights = 1.0 / err**2
sigma_0 = np.sum(weights * sigma) / np.sum(weights)
sigma_0_err = 1.0 / np.sqrt(np.sum(weights))

# Fit evaluation
dof = len(theta) - 1
chi2 = np.sum(((sigma - sigma_0) / err) ** 2)
chi2_red = chi2 / dof

print("---Singular Isothermal Sphere Fit---")
print(f"Best Fit Sigma 0: {sigma_0:.2f} ± {sigma_0_err:.2f} km/s")
print(f"Total Chi^2: {chi2:.2f}")
print(f"Degrees of Freedom: {dof}")
print(f"Reduced Chi^2: {chi2_red:.2f}")

plt.figure(figsize=(8, 5.5))

plt.errorbar(
    theta,
    sigma,
    yerr=err,
    fmt='o',
    color='black',
    ecolor='crimson',
    elinewidth=1.5,
    capsize=4,
    capthick=1.2,
    markersize=6,
    label=r'DESI Coma Profile ($\sigma_{\mathrm{int}}$)',
)

theta_grid = np.linspace(0.05, 2.0, 100)
plt.plot(
    theta_grid,
    np.full_like(theta_grid, sigma_0),
    'b--',
    linewidth=2,
    label=f'Isothermal Fit ($\\sigma_0 = {sigma_0:.1f}$ km/s, $\\chi_\\nu^2 = {chi2_red:.2f}$)',
)

plt.xlabel(r'Projected Radius $\theta$ [degrees]', fontsize=12)
plt.ylabel(
    r'Velocity Dispersion $\sigma_v$ [$\mathrm{km\ s^{-1}}$]', fontsize=12
)
plt.title('Coma Cluster: Observed Profile vs. Isothermal Sphere', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True, fontsize=10)
plt.tight_layout()
plt.show()

# NFW Profile
# 1. Define the empirical isotropic NFW velocity dispersion model
def sigma_nfw(theta, sigma_max, theta_s):
    x = theta / theta_s
    # Avoid zero-division at x -> 0 by enforcing numerical stability
    x = np.where(x < 1e-5, 1e-5, x)
    return sigma_max * np.sqrt(np.log(1.0 + x) / x)

# 2. Extract kinematic profile data from your profile_df
theta_data = profile_df['Theta_Deg'].to_numpy()
sigma_data = profile_df['Sigma_Int_Kms'].to_numpy()
err_data = profile_df['Sigma_Err_Kms'].to_numpy()

# 3. Fit the NFW model using non-linear least squares
# Initial guesses: sigma_max ~ 1100 km/s, theta_s ~ 0.4 degrees
p0 = [1100.0, 0.4]
bounds = ([500.0, 0.05], [2000.0, 2.0])

popt, pcov = curve_fit(
    sigma_nfw,
    theta_data,
    sigma_data,
    p0=p0,
    sigma=err_data,
    absolute_sigma=True,
    bounds=bounds
)

sigma_max_fit, theta_s_fit = popt
sigma_max_err, theta_s_err = np.sqrt(np.diag(pcov))

# 4. Goodness-of-fit (N = 8, p = 2 parameters -> DOF = 6)
dof_nfw = len(theta_data) - 2
residuals = (sigma_data - sigma_nfw(theta_data, *popt)) / err_data
chi2_nfw = np.sum(residuals**2)
chi2_red_nfw = chi2_nfw / dof_nfw

print("--- Navarro-Frenk-White (NFW) Fit ---")
print(f"Sigma_max        : {sigma_max_fit:.2f} ± {sigma_max_err:.2f} km/s")
print(f"Scale Radius (θs): {theta_s_fit:.3f} ± {theta_s_err:.3f} deg")
print(f"Total Chi^2      : {chi2_nfw:.2f}")
print(f"Degrees of Freedom: {dof_nfw}")
print(f"Reduced Chi^2    : {chi2_red_nfw:.2f}")
print("-" * 38)

# 5. Comparative Plot
plt.figure(figsize=(8.5, 5.5))

# Observed DESI data
plt.errorbar(
    theta_data,
    sigma_data,
    yerr=err_data,
    fmt='o',
    color='black',
    ecolor='crimson',
    elinewidth=1.5,
    capsize=4,
    capthick=1.2,
    markersize=6,
    label=r'DESI Coma Profile ($\sigma_{\mathrm{int}}$)'
)

theta_fine = np.linspace(0.05, 2.0, 200)

# SIS Baseline (Flat)
plt.plot(
    theta_fine,
    np.full_like(theta_fine, 946.99),
    'b--',
    linewidth=1.5,
    label=r'SIS Baseline ($\chi_\nu^2 = 3.23$)'
)

# NFW Fit (Declining)
plt.plot(
    theta_fine,
    sigma_nfw(theta_fine, *popt),
    'g-',
    linewidth=2.2,
    label=f'NFW Fit ($\\theta_s = {theta_s_fit:.2f}^\\circ$, $\\chi_\\nu^2 = {chi2_red_nfw:.2f}$)'
)

plt.xlabel(r'Projected Radius $\theta$ [degrees]', fontsize=12)
plt.ylabel(r'Velocity Dispersion $\sigma_v$ [$\mathrm{km\ s^{-1}}$]', fontsize=12)
plt.title('Coma Cluster: Isothermal Sphere vs. NFW Dark Matter Profile', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True, fontsize=10)
plt.tight_layout()
plt.show()

# Bayesian Information Criterion (BIC)
N = len(theta_data)

# SIS: k = 1 parameter (sigma_0)
chi2_sis = 22.61
bic_sis = chi2_sis + 1 * np.log(N)

# NFW: k = 2 parameters (sigma_max, theta_s)
chi2_nfw = 4.11
bic_nfw = chi2_nfw + 2 * np.log(N)

delta_bic = bic_sis - bic_nfw

print(f"BIC (Isothermal Sphere) : {bic_sis:.2f}")
print(f"BIC (NFW Profile)       : {bic_nfw:.2f}")
print(f"Delta BIC (SIS - NFW)   : {delta_bic:.2f}")



# Cored / Plummer Profile

# Define the soft-core / Plummer dispersion model
def sigma_core(theta, sigma_0, theta_c):
    return sigma_0 / (1.0 + (theta / theta_c)**2)**0.25

# Initial guess and bounds
p0_core = [1100.0, 1.0]
bounds_core = ([500.0, 0.1], [2000.0, 5.0])

popt_core, pcov_core = curve_fit(
    sigma_core,
    theta_data,
    sigma_data,
    p0=p0_core,
    sigma=err_data,
    absolute_sigma=True,
    bounds=bounds_core
)

sigma_0_c, theta_c_fit = popt_core
dof_core = len(theta_data) - 2
residuals_core = (sigma_data - sigma_core(theta_data, *popt_core)) / err_data
chi2_core = np.sum(residuals_core**2)
chi2_red_core = chi2_core / dof_core
bic_core = chi2_core + 2 * np.log(N)

print("--- Soft-Core Profile Fit ---")
print(f"Central Sigma_0  : {sigma_0_c:.2f} km/s")
print(f"Core Radius (θc) : {theta_c_fit:.3f} deg")
print(f"Total Chi^2      : {chi2_core:.2f}")
print(f"Reduced Chi^2    : {chi2_red_core:.2f}")
print(f"BIC (Core)       : {bic_core:.2f}")


# --- 3-way model comparison ---

theta_fine = np.linspace(0.05, 2.0, 200)

plt.figure(figsize=(9, 6))

# 1. DESI Data Points
plt.errorbar(
    theta_data,
    sigma_data,
    yerr=err_data,
    fmt='o',
    color='black',
    ecolor='crimson',
    elinewidth=1.5,
    capsize=4,
    capthick=1.2,
    markersize=6,
    label=r'DESI Coma Members ($\sigma_{\mathrm{int}}$)',
)

# 2. Singular Isothermal Sphere (Flat)
plt.plot(
    theta_fine,
    np.full_like(theta_fine, 946.99),
    'b--',
    linewidth=1.6,
    label=r'SIS Baseline (BIC = 24.69, $\chi_\nu^2 = 3.23$)',
)

# 3. Soft-Core Profile
sigma_core_line = 1027.35 / (1.0 + (theta_fine / 1.414) ** 2) ** 0.25
plt.plot(
    theta_fine,
    sigma_core_line,
    color='darkorange',
    linestyle='-.',
    linewidth=1.8,
    label=r'Soft-Core Fit (BIC = 11.63, $\chi_\nu^2 = 1.24$)',
)

# 4. Navarro-Frenk-White (Cuspy CDM)
x_fine = theta_fine / 1.047
sigma_nfw_line = 1093.30 * np.sqrt(np.log(1.0 + x_fine) / x_fine)
plt.plot(
    theta_fine,
    sigma_nfw_line,
    'g-',
    linewidth=2.3,
    label=r'NFW Cuspy Halo (BIC = 8.27, $\chi_\nu^2 = 0.69$)',
)

plt.xlabel(r'Projected Radius $\theta$ [degrees]', fontsize=12)
plt.ylabel(
    r'Velocity Dispersion $\sigma_v$ [$\mathrm{km\ s^{-1}}$]', fontsize=12
)
plt.title(
    'Coma Cluster: Model Selection (SIS vs. Cored vs. NFW)', fontsize=13
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True, fontsize=10)
plt.tight_layout()
plt.savefig('model_comparison.png', dpi=300)
plt.show()


# Total Virial Mass Estimation
# -------------------------------------------------------------
# 1. Physical Constants and Cluster Distance
# -------------------------------------------------------------
G = 4.30091e-9       # Gravitational constant in (km/s)^2 * Mpc / M_sun
c_light = 299792.458 # km/s
z_coma = 0.0231      # Mean redshift
H0 = 70.0            # km/s/Mpc

# Distance to Coma
D_Mpc = (c_light * z_coma) / H0  # ~98.9 Mpc

# -------------------------------------------------------------
# 2. Fitted NFW Parameters
# -------------------------------------------------------------
sigma_max = 1093.30  # km/s
theta_s_deg = 1.047  # degrees

# Convert scale radius from degrees to physical Mpc
r_s_Mpc = D_Mpc * np.deg2rad(theta_s_deg)  # ~1.81 Mpc

# NFW characteristic mass normalization M_0
# f_V ≈ 0.2162 is the maximum value of [ln(1+x) - x/(1+x)] / x
f_V = 0.2162
M_0 = (r_s_Mpc * (sigma_max**2)) / (G * f_V)

def mass_nfw(r_Mpc):
    """Enclosed mass M(<r) in Solar Masses (M_sun)"""
    x = r_Mpc / r_s_Mpc
    return M_0 * (np.log(1.0 + x) - (x / (1.0 + x)))

# -------------------------------------------------------------
# 3. Calculate Mass at Key Boundaries
# -------------------------------------------------------------
# Mass within our 2.0 degree survey aperture
r_max_deg = 2.0
r_max_Mpc = D_Mpc * np.deg2rad(r_max_deg)
M_2deg = mass_nfw(r_max_Mpc)

# Classical Virial Radius (R200 ~ 2.0 - 2.5 Mpc for Coma)
r_vir_Mpc = 2.0
M_vir = mass_nfw(r_vir_Mpc)

print(f"--- Coma Cluster Dynamical Mass Estimates ---")
print(f"Assumed Distance (D)       : {D_Mpc:.1f} Mpc")
print(f"Physical Scale Radius (rs) : {r_s_Mpc:.3f} Mpc ({r_s_Mpc*1000:.0f} kpc)")
print(f"Enclosed Mass within θ=2.0°: {M_2deg:.2e} M_sun")
print(f"Virial Mass (r ≈ 2.0 Mpc)  : {M_vir:.2e} M_sun")
print("-" * 45)

# -------------------------------------------------------------
# 4. Plot Enclosed Mass Profile M(<r)
# -------------------------------------------------------------
theta_eval = np.linspace(0.05, 2.0, 100)
r_eval = D_Mpc * np.deg2rad(theta_eval)
M_eval = mass_nfw(r_eval)

plt.figure(figsize=(8, 5.5))
plt.plot(theta_eval, M_eval / 1e14, 'purple', linewidth=2.2, label=r'$M_{\mathrm{NFW}}(<\theta)$')

plt.axvline(x=theta_s_deg, color='gray', linestyle=':', label=f'Scale Radius $r_s$ ({theta_s_deg:.2f}°)')
plt.xlabel(r'Projected Radius $\theta$ [degrees]', fontsize=12)
plt.ylabel(r'Enclosed Mass $M(<\theta)$ [$10^{14}\ M_\odot$]', fontsize=12)
plt.title('Coma Cluster Enclosed Mass Profile from DESI Kinematics', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(frameon=True, fontsize=11)
plt.tight_layout()
plt.savefig("coma_enclosed_mass_profile.png", dpi=300)
plt.show()

