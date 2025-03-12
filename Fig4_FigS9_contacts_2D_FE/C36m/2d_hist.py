#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

df_all = pd.read_csv('contacts_dimer_C36m.csv')/2 # divide by 2 to not double count
df = df_all.iloc[6000:] # 60ns of equilibration for C36m
intra_contacts = df['Intrachain']
inter_SPHF6_aligned = df['Interchain SPHF6 Aligned']

# make 2d histogram of intra and inter
def free_energy(a, b, T, y0, ymax, x0, xmax):
    free_energy, xedges, yedges = np.histogram2d(
        a, b, 30, [[y0, ymax], [x0, xmax]], density=True, weights=None)
    free_energy = np.log(np.flipud(free_energy)+.000001)
    free_energy = -(0.008314*T)*free_energy # for kJ/mol
    # use 0.008314 for kJ/mol
    return free_energy, xedges, yedges

fig,ax = plt.subplots() # figsize=(5.5,4.35) for annotated final version
fontsize=14

dG, xedges, yedges = free_energy(inter_SPHF6_aligned, intra_contacts, 300, 0, 30, 0, 65)
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
plt.ylim(0,25)

ax.set_title('jR2R3 P301L dimer: CHARMM36m',fontsize=fontsize)
ax.set_ylabel('interchain aligned PHF6 contacts',fontsize=fontsize)
ax.set_xlabel('intrachain contacts',fontsize=fontsize)
plt.subplots_adjust(bottom=0.15)
fig.savefig('2d_hist_C36m.png', dpi=300)

# Add rectangles
regions = {
    'i': ((0, 0), 60, 2),
    'ii': ((2.5, 7.5), 9, 3.5),
    'iii': ((21, 7.5), 11, 3.5),
    'iv': ((2.5, 16.5), 18, 6),
    'v': ((37, 15), 10, 3.5)
}

for label, (bottom_left, width, height) in regions.items():
    rect = patches.Rectangle(bottom_left, width, height, linewidth=1,
                             edgecolor='gray', facecolor='none', alpha=0.5)
    ax.add_patch(rect)
    center_x = bottom_left[0] + width / 2
    center_y = bottom_left[1] + height / 2
    color = 'white'
    if label=='iv': color='black'
    if label=='v': offset = 2.25
    else: offset = 0
    ax.annotate(label, xy=(center_x, center_y+offset), fontsize=12, color=color,
                ha='center', va='center',alpha=0.8)

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
    mask = (df_all['Interchain SPHF6 Aligned'] >= y0) & (df_all['Interchain SPHF6 Aligned'] <= y1) & \
           (df_all['Intrachain'] >= x0) & (df_all['Intrachain'] <= x1)
    masks[region] = mask

# Convert the masks dictionary to DataFrame and save to CSV
masks_df = pd.DataFrame(masks)
masks_df.sum(axis=0)
masks_df.to_csv('region_masks.csv', index=False)

#%%
# sum each column of masks df
print(masks_df.sum(axis=0))
#%%
