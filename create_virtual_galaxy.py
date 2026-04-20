"""
Generates a pure-vacuum "Virtual Galaxy" dataset formatted for the 
Linear Shear 3D Batch Processor.
"""
import pandas as pd
import numpy as np

# --- VIRTUAL GALAXY PARAMETERS ---
K = 150.0      # The terminal flat-line velocity (km/s)
c = 150000     # The scale factor (how slowly it rises)
max_radius = 40.0 # Calculate out to 40 kpc

print("Synthesizing Virtual Pure-Vacuum Galaxy...")

# 1. Generate the perfect theoretical equilibrium curve
# We sweep V from 1 up to 99% of K to avoid dividing by zero
V_ideal = np.linspace(1, K * 0.99, 500)

# Calculate Radius using the inverse kinematic equilibrium equation
R_ideal = c * V_ideal / (K - V_ideal)**3

# Filter out extreme mathematically abstract bounds to fit a standard plot
valid_indices = (R_ideal >= 0.1) & (R_ideal <= max_radius)
R_sim = R_ideal[valid_indices]
V_sim = V_ideal[valid_indices]

# 2. Package it into a SPARC-compatible DataFrame
# Notice that ALL baryonic mass (gas, disk, bulge) is strictly zero.
df_virtual = pd.DataFrame({
    'Galaxy': ['VIRTUAL'] * len(R_sim),
    'Rad': R_sim,
    'Vobs': V_sim,
    'errV': [1.0] * len(R_sim),  # Dummy standard error
    'Vgas': [0.0] * len(R_sim),  # ZERO Baryons
    'Vdisk': [0.0] * len(R_sim), # ZERO Baryons
    'Vbulge': [0.0] * len(R_sim) # ZERO Baryons
})

# 3. Save as CSV
filename = 'virtual_vacuum_galaxy.csv'
df_virtual.to_csv(filename, index=False)
print(f"Success! Saved pure vacuum kinematics to '{filename}'.")