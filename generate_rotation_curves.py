"""
=============================================================================
Linear Shear Vacuum Model - 3D Batch Rotation Curve Generator
=============================================================================
Features included:
- Dynamic Mass-to-Light ratio overrides per galaxy via config.yml
- Manual Scale-Height (hv_fixed) overrides for extreme LSBs
- Bulletproof regex string matching to handle formatting errors in SPARC
=============================================================================
"""

import os
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
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
    
    # Extract dynamic parameters
    u_disk = target_config.get('upsilon_disk', UPSILON_DISK)
    u_bulge = target_config.get('upsilon_bulge', UPSILON_BULGE)
    description = target_config.get('description', 'Standard SPARC Galaxy')
    
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
    # NATIVE OVERRIDE & OPTIMIZATION
    # =========================================================
    if 'hv_fixed' in target_config:
        best_hv = float(target_config['hv_fixed'])
        print(f"  -> Manual Override: using hv = {best_hv} for {galaxy_name}")
    else:
        # Bounds expanded to 35 to allow organic spherical discovery
        res = minimize_scalar(
            lambda hv: evaluate_galaxy(hv, R_grid, z_grid, R_2d, z_2d, V_midplane, masks, kernels, R_obs, V_true)[0],
            bounds=(1.0, 35.0), method='bounded'
        )
        best_hv = res.x
    
    _, V_syn = evaluate_galaxy(best_hv, R_grid, z_grid, R_2d, z_2d, V_midplane, masks, kernels, R_obs, V_true)
    V_tot_predicted = np.sqrt(V_bar_sq + V_syn**2)

    # PLOTTING
    plt.figure(figsize=(12, 7))
    plt.plot(R_obs, V_obs, 'ko-', linewidth=2, label='Observed Total Velocity ($V_{obs}$)', markersize=6)
    plt.plot(R_obs, V_bar, 'b--', linewidth=2, label='Visible Matter Contribution ($V_{bar}$)')
    plt.fill_between(R_obs, V_bar, V_obs, color='crimson', alpha=0.15, label='The "Missing Mass" Anomaly')
    
    plt.plot(R_obs, V_syn, 'm-.', linewidth=2.5, label=f'Model Predicted Vacuum Mass ($V_{{syn}}$)')
    plt.plot(R_obs, V_tot_predicted, 'g-', linewidth=3, alpha=0.7, label='Total Model Prediction ($V_{syn} + V_{bar}$)')

    max_y = max(V_obs) * 1.3 
    
    plt.title(f'{galaxy_name} 3D Rotation Curve Breakdown\nTesting the Linear Shear Vacuum Model', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Radial Distance (kpc)', fontsize=14)
    plt.ylabel('Orbital Velocity (km/s)', fontsize=14)
    
    context_text = (
        f"{description}\n"
        f"Optimized Disk Thickness ($hv$): {best_hv:.2f} kpc   |   "
        f"Global Constant ($A$): {GLOBAL_A:.2e}   |   "
        f"Mass-to-Light: $\\Upsilon_{{disk}}$={u_disk}, $\\Upsilon_{{bulge}}$={u_bulge}"
    )
    
    plt.text(0.5, 0.96, context_text, transform=plt.gca().transAxes, fontsize=11,
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
    
    print(f"  -> Saved 3D plot for {galaxy_name} (hv = {best_hv:.2f} kpc): {output_path}")

if __name__ == "__main__":
    config = load_config('galaxies_to_generate.yml')
    if not config: exit()
        
    data_file = config.get('data_file', 'cdsarc_152_157_table2.xlsx')
    output_folder = config.get('output_folder', 'Shear_Analysis_Plots')
    target_galaxies = config.get('target_galaxies', [])
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        
    print(f"Loading SPARC Database from {data_file}...")
    try:
        df_sparc = pd.read_excel(data_file) if data_file.endswith('.xlsx') else pd.read_csv(data_file)
    except FileNotFoundError:
        print(f"Error: Data file '{data_file}' not found.")
        exit()
        
    col_name = 'Name' if 'Name' in df_sparc.columns else 'Galaxy' if 'Galaxy' in df_sparc.columns else 'Galaxy identifier' if 'Galaxy identifier' in df_sparc.columns else None
    
    df_sparc['Clean_Galaxy'] = df_sparc[col_name].astype(str).str.replace(r'\W+', '', regex=True).str.upper()

    df_sparc = df_sparc.rename(columns={
        'rad': 'Rad', 'vobs': 'Vobs', 'errv': 'errV', 
        'vgas': 'Vgas', 'vdisk': 'Vdisk', 'vbulge': 'Vbulge'
    })

    print(f"Processing {len(target_galaxies)} galaxies via 3D Integration...")
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