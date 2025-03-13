#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d

# Function to calculate free energy
def free_energy(a, b, T, y0, ymax, x0, xmax):
    hist, xedges, yedges = np.histogram2d(a, b, bins=30, range=[[y0, ymax], [x0, xmax]], density=True)
    free_energy = -0.008314 * T * np.log(hist + 1e-6)  # Avoid log(0)
    free_energy = free_energy - np.min(free_energy)  # Normalize to 0
    return free_energy, xedges, yedges

# Load data
def load_data(csv_file):
    df = pd.read_csv(csv_file)
    return df['Rg'], df['Ree']

# File paths
csv_files = {
    'C36m': 'rg_ree_C36m.csv',
    'a99SBdisp': 'rg_ree_a99SBdisp.csv',
    'a03ws': 'rg_ree_a03ws.csv'
}
from matplotlib import gridspec

fig = plt.figure(figsize=(15, 5))
gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 0.05])  # Last column for colorbar
axs = [plt.subplot(gs[i]) for i in range(3)]  # Create subplots for the data

T = 300  # Temperature in Kelvin

for i, (label, file) in enumerate(csv_files.items()):
    rg, ree = load_data(file)
    rg = rg.values[5000:] # skipping 50ns for equilibration
    ree = ree.values[5000:] # skipping 50ns for equilibration
    
    # Compute free energy
    dG, xedges, yedges = free_energy(ree, rg, T, 0.2, 5, 0.5, 1.8)
    
    # Create heatmap
    im = axs[i].imshow(dG.T, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                       cmap='jet', aspect='auto', vmax=14, interpolation='gaussian')
    
    axs[i].set_title(label, fontsize=20)
    axs[i].set_xlabel(r'$r_{ee}$ (nm)', fontsize=16)
    axs[i].set_ylabel(r'$r_g$ (nm)', fontsize=16)
    axs[i].tick_params(axis='both', which='major', labelsize=16)

    # plot a star at three points
    color1 = 'white'
    color2 = 'white'
    if i == 0:
        axs[i].plot(1.57, 0.84, '*', color=color1, markerfacecolor='none', markersize=10) # cluster 1
        axs[i].plot(2.91, 1.40, '*', color=color1, markerfacecolor='none', markersize=10) # cluster 2
        axs[i].plot(1.41, 0.98, '*', color=color1, markerfacecolor='none', markersize=10) # cluster 3
    else: 
        axs[i].plot(1.57, 0.85, '*', color=color2, markerfacecolor='none', markersize=10)
        axs[i].plot(2.91, 1.40, '*', color=color2, markerfacecolor='none', markersize=10)
        axs[i].plot(1.41, 0.98, '*', color=color2, markerfacecolor='none', markersize=10)

# Add colorbar in the 4th column
cbar_ax = plt.subplot(gs[3])  # Colorbar axis
cbar = fig.colorbar(im, cax=cbar_ax, pad=0.0)
cbar.set_label('free energy (kJ/mol)', fontsize=16)
cbar.ax.tick_params(labelsize=16)

plt.tight_layout()
plt.savefig('rg_ree.png', dpi=300, bbox_inches='tight')
plt.show()
