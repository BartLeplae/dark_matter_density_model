"""
Generates a multi-target pure-vacuum dataset to test 
kinematic and geometric parameters in the 3D Shear model.
Parameters are dynamically loaded from virtual_galaxy_to_generate.yml
"""
import pandas as pd
import numpy as np
import yaml
import sys

def load_config(config_path: str):
    """Loads settings from the YAML configuration file."""
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: Could not find configuration file at '{config_path}'.")
        return None

def generate_virtual_galaxy(name, K, c, max_radius=40.0):
    """Generates the perfect pure-vacuum equilibrium kinematic dataset."""
    V_ideal = np.linspace(1, K * 0.99, 500)
    R_ideal = c * V_ideal / (K - V_ideal)**3
    
    valid_indices = (R_ideal >= 0.1) & (R_ideal <= max_radius)
    R_sim = R_ideal[valid_indices]
    V_sim = V_ideal[valid_indices]
    
    return pd.DataFrame({
        'Galaxy': [name] * len(R_sim),
        'Rad': R_sim,
        'Vobs': V_sim,
        'errV': [1.0] * len(R_sim),
        'Vgas': [0.0] * len(R_sim),
        'Vdisk': [0.0] * len(R_sim),
        'Vbulge': [0.0] * len(R_sim)
    })

# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    config_file = 'virtual_galaxies_to_generate.yml'
    config = load_config(config_file)
    if not config:
        sys.exit()

    target_galaxies = config.get('target_galaxies', [])
    data_file = config.get('data_file', 'virtual_vacuum_galaxies.csv')

    print(f"Synthesizing Virtual Galaxies from {config_file}...")

    dfs = []
    for target in target_galaxies:
        name = target.get('name')
        K = target.get('K')
        c = target.get('c')
        
        # Only generate galaxies that have the required kinematic parameters
        if K is not None and c is not None:
            max_r = target.get('max_radius', 40.0) # Defaults to 40 kpc if not specified
            print(f"  -> Generating {name} (K={K}, c={c})")
            df_g = generate_virtual_galaxy(name, K, c, max_r)
            dfs.append(df_g)
        else:
            print(f"  -> Skipping '{name}' (Missing K or c generation parameters in YAML)")

    if dfs:
        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_csv(data_file, index=False)
        print(f"\nSuccess! Saved {len(df_all)} data points to '{data_file}'.")
    else:
        print("\nNo valid virtual galaxies generated. Check YAML parameters.")