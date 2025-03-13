#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('contacts_dimer_a99SBdisp.csv')/2 # avoid double counting
df = df.iloc[4000:] # 40ns of equilibration for a99SBdisp
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

fig,ax = plt.subplots()

dG, xedges, yedges = free_energy(inter_SPHF6_aligned, intra_contacts, 300, 0, 30, 0, 65)
dG = dG - np.min(dG)
im = plt.imshow(dG, interpolation='gaussian', extent=[
            yedges[0], yedges[-1], xedges[0], xedges[-1]], cmap='jet',
            aspect='auto',vmax=20)

# FORMAT THE COLORBAR
cb = fig.colorbar(im)
cb.ax.get_yaxis().labelpad = 12
cb.ax.set_ylabel('kJ/mol',fontsize=16)
cbar_ticks = [0, 5, 10, 15, 20]
cb.set_ticks(cbar_ticks)

ax.set_title('jR2R3 P301L dimer: a99SBdisp',fontsize=16)
ax.set_ylabel('interchain aligned PHF6 contacts',fontsize=16)
ax.set_xlabel('intrachain contacts',fontsize=16)
cb.ax.tick_params(labelsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylim([0,25])

plt.subplots_adjust(bottom=0.15)
plt.savefig('2dhist_a99SBdisp.png',dpi=300)
plt.show()
