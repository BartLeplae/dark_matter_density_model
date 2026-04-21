"""
=============================================================================
Linear Shear Vacuum Model - Virtual 3D Rotation Curve Generator
=============================================================================
Features included:
- Hardcoded virtual parameters (K, c, hv_fixed) driven by config.yml
- Bulletproof regex string matching to handle formatting errors
- Dynamic, all-inclusive parameter text box on generated plots
=============================================================================
"""

import os
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import re

warnings.filterwarnings('ignore')

# =============================================================================
# UNIVERSAL PHYSICS CONSTANTS (DEFAULTS)
# =============================================================================
GLOBAL_A = 1.74e-03  
GLOBAL_B = 1.0000    
G = 4.3009e-6        
UPSILON_DISK = 0.5   
UPSILON_BULGE = 0.7  

def load_config(config_path: str):
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: Could not find configuration file at '{config_path}'.")
        return None

def evaluate_galaxy(hv: float, R_grid: np.ndarray, z_grid: np.ndarray, 
                    R_2d: np.ndarray, z_2d: np.ndarray, V_midplane: np.ndarray, 
                    masks: dict, kernels: list, Rad_obs: np.ndarray, V_true: np.ndarray) -> tuple:
    
    V_2d = V_midplane * np.exp(-z_2d / hv)
    Omega_2d = V_2d / R_2d
    
    Shear_2d = np.sqrt(np.gradient(Omega_2d, R_grid, axis=0)**2 + np.gradient(Omega_2d, z_grid, axis=1)**2)
    Rho_2d = GLOBAL_A * (Shear_2d ** GLOBAL_B)
    
    Rho_2d[masks['outer']] *= masks['outer_taper']
    Rho_2d[np.sqrt(R_2d**2 + z_2d**2) > masks['edge']] = 0.0
    Rho_kpc = Rho_2d[:, :, None] * 1e9
    
    V_syn = np.array([np.sqrt(np.maximum(0, R * np.sum(Rho_kpc * kernels[i]))) for i, R in enumerate(Rad_obs)])
    mse = np.mean((V_syn - V_true)**2) + (1.0 * (hv - 3.0)**2)
    return mse, V_syn

def process_galaxy_3d(df_target: pd.DataFrame, target_config: dict, output_folder: str):
    galaxy_name = target_config['name']
    
    # Extract dynamic parameters (Assuming K, c, and hv_fixed are ALWAYS present)
    u_disk = target_config.get('upsilon_disk', UPSILON_DISK)
    u_bulge = target_config.get('upsilon_bulge', UPSILON_BULGE)
    description = target_config.get('description', 'Standard Virtual Galaxy')
    
    K_val = target_config['K']
    c_val = target_config['c']
    best_hv = float(target_config['hv_fixed'])
    
    df_target = df_target[df_target['Rad'] > 0.01].copy()
    df_target = df_target.sort_values(by='Rad').reset_index(drop=True)
    
    if len(df_target) < 3:
        print(f"  -> Skipping {galaxy_name}: Not enough data points.")
        return

    R_obs = df_target['Rad'].values
    V_obs = df_target['Vobs'].values
    
    V_disk_sq = u_disk * (pd.to_numeric(df_target['Vdisk']).values ** 2)
    V_bulge_sq = u_bulge * (pd.to_numeric(df_target['Vbulge']).values ** 2)
    V_gas_sq = pd.to_numeric(df_target['Vgas']).values * np.abs(pd.to_numeric(df_target['Vgas']).values)
    
    V_bar_sq = V_disk_sq + V_bulge_sq + V_gas_sq
    V_bar_sq[V_bar_sq < 0] = 0
    V_bar = np.sqrt(V_bar_sq)
    
    V_true_sq = (V_obs ** 2) - V_bar_sq
    V_true_sq[V_true_sq < 0] = 0
    V_true = np.sqrt(V_true_sq)
    
    # Pre-computation
    R_max = R_obs.max()
    R_grid = np.linspace(max(0.1, R_obs.min()), R_max, 50)
    z_grid = np.linspace(0.0, R_max * 1.5, 40)
    phi_grid = np.linspace(0, 2 * np.pi, 50)
    
    R_3d, z_3d, phi_3d = np.meshgrid(R_grid, z_grid, phi_grid, indexing='ij')
    R_2d, z_2d = np.meshgrid(R_grid, z_grid, indexing='ij')
    
    V_midplane = np.interp(R_grid, R_obs, V_obs)[:, None]
    
    r_sph2d = np.sqrt(R_2d**2 + z_2d**2)
    masks = {
        'outer': r_sph2d > R_max * 0.95,
        'outer_taper': np.exp(-(r_sph2d[r_sph2d > R_max * 0.95] - R_max * 0.95) / 10.0),
        'edge': R_max * 1.3
    }
    
    dV = R_3d * (R_grid[1]-R_grid[0]) * (z_grid[1]-z_grid[0]) * (phi_grid[1]-phi_grid[0])
    dr_max = max(R_grid[1]-R_grid[0], z_grid[1]-z_grid[0]) * 2.0
    
    kernels = []
    for R in R_obs:
        if R > 0:
            denom = (R_3d**2 + R**2 - 2*R_3d*R*np.cos(phi_3d) + z_3d**2 + dr_max**2)**1.5
            k = 2 * G * dV * (R - R_3d * np.cos(phi_3d)) / denom
            kernels.append(k)
        else:
            kernels.append(np.zeros_like(R_3d))
            
    # =========================================================
    # PHYSICS INTEGRATION (Using Forced hv)
    # =========================================================
    print(f"  -> Using fixed hv = {best_hv} for {galaxy_name}")
    _, V_syn = evaluate_galaxy(best_hv, R_grid, z_grid, R_2d, z_2d, V_midplane, masks, kernels, R_obs, V_true)
    V_tot_predicted = np.sqrt(V_bar_sq + V_syn**2)

    # =========================================================
    # PLOTTING
    # =========================================================
    plt.figure(figsize=(12, 7))
    plt.plot(R_obs, V_obs, 'ko-', linewidth=2, label='Observed Total Velocity ($V_{obs}$)', markersize=6)
    plt.plot(R_obs, V_bar, 'b--', linewidth=2, label='Visible Matter Contribution ($V_{bar}$)')
    plt.fill_between(R_obs, V_bar, V_obs, color='crimson', alpha=0.15, label='The "Missing Mass" Anomaly')
    
    plt.plot(R_obs, V_syn, 'm-.', linewidth=2.5, label=f'Model Predicted Vacuum Mass ($V_{{syn}}$)')
    plt.plot(R_obs, V_tot_predicted, 'g-', linewidth=3, alpha=0.7, label='Total Model Prediction ($V_{syn} + V_{bar}$)')

    # Add extra space at the top so the curves don't hit the larger text box
    max_y = max(V_obs) * 1.40 
    
    plt.title(f'{galaxy_name} 3D Rotation Curve Breakdown\nTesting the Linear Shear Vacuum Model', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Radial Distance (kpc)', fontsize=14)
    plt.ylabel('Orbital Velocity (km/s)', fontsize=14)
    
    # =========================================================
    # CALCULATE COSMOLOGICAL MASS & ENERGY
    # =========================================================
    # 1. Total Vacuum Energy (E = K^2 * c / G)
    E_total = ((K_val**2) * c_val) / G
    
    # 2. Apparent Missing Mass at R_max (M = V^2 * R / G)
    # We use V_true to isolate only the dark matter / vacuum velocity requirement
    M_missing_max_billions = ((V_true[-1]**2) * R_max / G) / 1e9

# --- DYNAMIC PARAMETER BOX BUILDER ---
    param_lines = [f"{description}"]
    
    # Line 2: Universal Physics (A stands alone)
    physics_str = f"Vacuum-Shear Coupling Constant ($A$): {GLOBAL_A:.2e} (Kinematic-to-Mass Conversion)"
    param_lines.append(physics_str)
    
    # Line 3: The Generation Variables & Geometry (Disk Thickness moved here)
    gen_str = (
        f"Terminal Velocity ($K$) = {K_val} km/s   |   "
        f"Energy Scale ($c$) = {c_val} kpc$\\cdot$(km/s)$^2$   |   "
        f"Disk Thickness ($hv$): {best_hv:.2f} kpc"
    )
    param_lines.append(gen_str)
    
    # Line 4: The Physical Consequences
    cosmo_str = (
        f"Total Vacuum Energy ($E = \\frac{{K^2 \\cdot c}}{{G}}$) = {E_total:.2e} $M_\\odot\\cdot$(km/s)$^2$   |   "
        f"Apparent Missing Mass (at {R_max:.1f} kpc) = {M_missing_max_billions:.1f} Billion $M_\\odot$"
    )
    param_lines.append(cosmo_str)
    
    context_text = "\n".join(param_lines)
    
    plt.text(0.5, 0.98, context_text, transform=plt.gca().transAxes, fontsize=11,
             ha='center', va='top', bbox=dict(facecolor='white', alpha=0.95, edgecolor='lightgrey', boxstyle='round,pad=0.6'))    
    plt.legend(loc='lower right', fontsize=11)
    plt.grid(True, alpha=0.3, linestyle='--')
    
    plt.xlim(0, max(R_obs) * 1.05)
    plt.ylim(0, max_y)
    plt.tight_layout()
    
    safe_name = re.sub(r'\W+', '_', galaxy_name)
    output_path = os.path.join(output_folder, f'{safe_name}_rotation_curves.png')
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    print(f"  -> Saved 3D plot for {galaxy_name}: {output_path}")

if __name__ == "__main__":
    import sys
    # Read the default config file, or allow passing one via command line
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'virtual_galaxies_to_generate.yml'
    
    config = load_config(config_file)
    if not config: exit()
        
    data_file = config.get('data_file', 'cdsarc_152_157_table2.xlsx')
    output_folder = config.get('output_folder', 'Shear_Analysis_Plots')
    target_galaxies = config.get('target_galaxies', [])
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    print(f"Loading Dataset from {data_file}...")
    try:
        df_sparc = pd.read_excel(data_file) if data_file.endswith('.xlsx') else pd.read_csv(data_file)
    except FileNotFoundError:
        print(f"Error: Data file '{data_file}' not found.")
        exit()
        
    col_name = 'Name' if 'Name' in df_sparc.columns else 'Galaxy' if 'Galaxy' in df_sparc.columns else 'Galaxy identifier' if 'Galaxy identifier' in df_sparc.columns else None
    
    if col_name:
        df_sparc['Clean_Galaxy'] = df_sparc[col_name].astype(str).str.replace(r'\W+', '', regex=True).str.upper()
    else:
        print("Error: Could not identify Galaxy Name column.")
        exit()
    
    print(f"Processing {len(target_galaxies)} virtual galaxies via 3D Integration...")
    for target in target_galaxies:
        target_config = {'name': target} if isinstance(target, str) else target
        target_name = target_config['name']
        
        clean_target = re.sub(r'\W+', '', str(target_name)).upper()
        df_target = df_sparc[df_sparc['Clean_Galaxy'] == clean_target].copy()
        
        if df_target.empty:
            print(f"  -> Warning: Galaxy '{target_name}' not found. Skipping.")
        else:
            process_galaxy_3d(df_target, target_config, output_folder)
            
    print("\nBatch processing complete!")