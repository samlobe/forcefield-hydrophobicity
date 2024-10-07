#%%
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

forcefields = ['a03ws', 'a99SBdisp', 'C36m']
types = ['monomer', 'dimer']
colors = {'a03ws': 'red', 'a99SBdisp': 'blue', 'C36m': 'green'}
hatches = {'monomer': '', 'dimer': '//'}

# Create an empty DataFrame to store all data
all_data = pd.DataFrame()

# Read and concatenate the data
for forcefield in forcefields:
    for type in types:
        data = pd.read_csv(f'{forcefield}_{type}.csv', index_col=0)
        if type == 'dimer':
            data = data / 2
        # Add a column for forcefield and type
        data['Forcefield'] = forcefield
        data['Type'] = type
        all_data = pd.concat([all_data, data])

# Set up the plot
fig, ax = plt.subplots(figsize=(7, 5))

# Create a list to hold the positions of each boxplot
positions = []

# Create a mapping of forcefields to their base positions
base_positions = {'a03ws': 1, 'a99SBdisp': 3, 'C36m': 5}

# Create the boxplots
for i, forcefield in enumerate(forcefields):
    for j, type in enumerate(types):
        type_data = all_data[(all_data['Forcefield'] == forcefield) & (all_data['Type'] == type)]
        # Plot each type next to each other
        pos = base_positions[forcefield] + j * 0.5
        positions.append(pos)
        box = ax.boxplot(type_data.iloc[:, :-2].values, positions=[pos], widths=0.3, patch_artist=True,
                         showfliers=False)
        # ax.boxplot()
        # Set colors and hatches
        for patch in box['boxes']:
            patch.set_facecolor(colors[forcefield])
            patch.set_hatch(hatches[type])
        
        # Set the mean line color to black
        # for mean in box['means']:
        #     mean.set_color('black')

# Customize the plot
ax.set_xticks([1.25, 3.25, 5.25])
ax.set_xticklabels(['a03ws', 'a99SBdisp', 'C36m'],fontsize=14)
ax.set_xlim(0, 7)
ax.set_xlabel('Forcefield',fontsize=14)
ax.set_ylabel('# Hydration Waters',fontsize=14)

# Add a legend
legend_handles = [Patch(facecolor='white', edgecolor='black', hatch=hatches['monomer'], label='Monomer'),
                  Patch(facecolor='white', edgecolor='black', hatch=hatches['dimer'], label='Dimer')]
ax.legend(handles=legend_handles)

# get mean of cryo_dimer.csv
cryo_dimer = pd.read_csv('cryo_dimer.csv',index_col=0)
cryo_dimer_mean = np.mean(cryo_dimer.values) / 2 # per monomer

# plot dashed horizontal line
plt.axhline(y=cryo_dimer_mean,color='black')

# get mean of C36m_dimer_cluster.csv
C36m_dimer_cluster = pd.read_csv('C36m_dimer_cluster.csv',index_col=0)
C36m_dimer_cluster_mean = np.mean(C36m_dimer_cluster.values) / 2 # per monomer

# plot dasehd horizontal line 
plt.axhline(y=C36m_dimer_cluster_mean,color='green')
# actually I should make an arrow pointing to C36m dimer at this height
# and label it cluster iv

# Maybe instead of a box-and-whisker plot I should show the monomer and dimer distros side-by-side
# in the style of graphs that Evan had done.

# Indicate a per monomer # of waters for infinitely long filament (strand A)


# Show the plot
# plt.show()
