#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

forcefields = ['a03ws','a99SBdisp','C36m']
colors = {'a03ws':'tab:red','a99SBdisp':'#039dfc','C36m':'black'}

all_data = pd.DataFrame()

for forcefield in forcefields:
    for type in ['monomer','dimer']:
        data = pd.read_csv(f'{forcefield}_{type}.csv',index_col=0)
        if type == 'dimer':
            data = data / 2 # per monomer
        data['Forcefield'] = forcefield
        data['Type'] = type
        all_data = pd.concat([all_data, data])

melted_data = all_data.melt(id_vars=['Forcefield','Type'],var_name='Variable',value_name='Value')

fig,ax = plt.subplots(figsize=(6,4))

sns.violinplot(x='Forcefield', y='Value', hue='Type', data=melted_data,
               split=True, palette=[colors['a03ws'], 'red', colors['a99SBdisp'], 'cyan', colors['C36m'], 'cyan'],
               inner=None, ax=ax)

# Customize each violin color
for i, violin in enumerate(ax.collections):
    if i % 2 == 0:
        forcefield = forcefields[i // 2]
        violin.set_facecolor(colors[forcefield])
        violin.set_edgecolor('gray')
    else:
        forcefield = forcefields[i // 2]
        violin.set_facecolor(colors[forcefield])
        violin.set_edgecolor('gray')
        violin.set_hatch('//')


# get mean of cryo_dimer.csv
cryo_dimer = pd.read_csv('cryo_dimer.csv',index_col=0)
cryo_dimer_mean = np.mean(cryo_dimer.values) / 2 # per monomer

# plot dashed horizontal line
plt.axhline(y=cryo_dimer_mean,color='gray',linestyle='--')
# add annotation
ax.text(0.5,cryo_dimer_mean-20,'dimer from\ncryo-EM structure',color='gray',fontsize=14,ha='center',va='center')

# get mean of C36m_dimer_cluster.csv
C36m_dimer_cluster = pd.read_csv('C36m_dimer_cluster.csv',index_col=0)
C36m_dimer_cluster_mean = np.mean(C36m_dimer_cluster.values) / 2 # per monomer

plt.xlim(right=3.3); plt.ylim(bottom=90)
# Draw a star for the comopact cluster's avg # waters
plt.scatter(2,C36m_dimer_cluster_mean,marker='*',color='purple',s=100)

# add arrow pointing left at C36m_dimer_cluster_mean
start_x = 2.7; end_x = 2.1
y_position = C36m_dimer_cluster_mean
ax.annotate('',xy=(end_x,y_position), xytext=(start_x,y_position),
            arrowprops=dict(arrowstyle='->',color='purple'))

# add annotation
ax.text(2.94,C36m_dimer_cluster_mean+18,'compact\ndimer\ncluster',fontsize=13,color='purple',
        ha='center',va='center')

# add a horizontal dotted line at 59.32
plt.axhline(y=59.32,color='black',linestyle=':')
plt.ylim(bottom=55)

# add annotation
ax.text(0.5,59.32+10,'full cryo-EM structure',color='black',fontsize=14,ha='center',va='center')

# have y-ticks from 60 to 260 counting by 40
plt.yticks(np.arange(60,300,40),fontsize=14)
plt.xticks(fontsize=14)
plt.xlabel('force field                ',fontsize=14)
plt.ylabel('hydration waters per peptide',fontsize=14)

# Add a legend
hatches = {'monomer': '', 'dimer': '//'}
legend_handles = [Patch(facecolor='white', edgecolor='black', hatch=hatches['monomer'], label='monomer'),
                  Patch(facecolor='white', edgecolor='black', hatch=hatches['dimer'], label='dimer')]
ax.legend(handles=legend_handles,fontsize=12,loc='upper right')

# rename x-ticks
ax.set_xticklabels(['a03ws','a99SBdisp','     CHARMM36m'])

# save the figure as png and 300 dpi
plt.subplots_adjust(bottom=0.15)
# plt.tight_layout()
plt.savefig('violinplot.png',dpi=300)

# %%
# compare the means of the monomers and dimers in each force field
# all_data.groupby(['Forcefield','Type']).mean()
all_data.groupby(['Forcefield','Type']).median()
plt.show()
