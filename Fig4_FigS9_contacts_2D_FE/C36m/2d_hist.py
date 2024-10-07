#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

df = pd.read_csv('contacts_dimer_C36m.csv')
intra_contacts = df['Intrachain']/2
inter_SPHF6_aligned = df['Interchain SPHF6 Aligned']/2

# make 2d histogram of intra and inter
def free_energy(a, b, T, y0, ymax, x0, xmax):
    free_energy, xedges, yedges = np.histogram2d(
        a, b, 30, [[y0, ymax], [x0, xmax]], density=True, weights=None)
    free_energy = np.log(np.flipud(free_energy)+.000001)
    free_energy = -(0.008314*T)*free_energy # for kJ/mol
    # use 0.008314 for kJ/mol
    return free_energy, xedges, yedges

fig,ax = plt.subplots()
fontsize=16

dG, xedges, yedges = free_energy(inter_SPHF6_aligned, intra_contacts, 300, 0, 25, 0, 65)
dG = dG - np.min(dG)
im = plt.imshow(dG, interpolation='gaussian', extent=[
            yedges[0], yedges[-1], xedges[0], xedges[-1]], cmap='jet',
            aspect='auto',vmax=20)

# FORMAT THE COLORBAR
cb = fig.colorbar(im)
cb.ax.get_yaxis().labelpad = 12
cb.ax.set_ylabel('kJ/mol',fontsize=fontsize)
cbar_ticks = [0, 5, 10, 15, 20]
cb.set_ticks(cbar_ticks)
cb.ax.tick_params(labelsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

ax.set_title('jR2R3 P301L dimer: CHARMM36m',fontsize=fontsize)
ax.set_ylabel('interchain aligned PHF6 contacts',fontsize=fontsize)
ax.set_xlabel('intrachain contacts',fontsize=fontsize)
plt.subplots_adjust(bottom=0.15)
fig.savefig('2d_hist_C36m.png', dpi=300)

# Add rectangles
regions = {
    'i': ((0, 0), 60, 2),
    'ii': ((21, 5), 11, 6),
    'iii': ((2.5, 7.5), 12.5, 3.5),
    'iv': ((2.5, 16), 27.5, 5),
    'v': ((40, 12.5), 15, 4.5)
}

for label, (bottom_left, width, height) in regions.items():
    rect = patches.Rectangle(bottom_left, width, height, linewidth=1,
                             edgecolor='gray', facecolor='none', alpha=0.5)
    ax.add_patch(rect)
    center_x = bottom_left[0] + width / 2
    center_y = bottom_left[1] + height / 2
    ax.annotate(label, xy=(center_x, center_y), fontsize=12, color='black',
                ha='center', va='center')

plt.show()
# save fig
fig.savefig('2d_hist_C36m_annotated.png', dpi=300)

#%%
# output masks for each region
# Initialize a dictionary to hold the masks
masks = {}

# Create masks based on the region specifications
for region, ((x0, y0), width, height) in regions.items():
    x1, y1 = x0 + width, y0 + height
    mask = (df['Interchain SPHF6 Aligned'] >= y0) & (df['Interchain SPHF6 Aligned'] <= y1) & \
           (df['Intrachain'] >= x0) & (df['Intrachain'] <= x1)
    masks[region] = mask

# Convert the masks dictionary to DataFrame and save to CSV
masks_df = pd.DataFrame(masks)
masks_df.sum(axis=0)
masks_df.to_csv('region_masks.csv', index=False)

#%%
