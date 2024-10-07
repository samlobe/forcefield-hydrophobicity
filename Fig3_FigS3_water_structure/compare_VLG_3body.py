#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ffs = ['a03ws','a99SBdisp','C36m']
dfs = []
for ff in ffs:
    df = pd.read_csv(f'VLG_3body_{ff}.csv',index_col=0)
    dfs.append(df)

# put into dictionary
data = {}
for i,ff in enumerate(ffs):
    data[ff] = dfs[i]

#%%
# measure tetrahedral fraction
processed_data = {}
for ff in ffs:
    df = data[ff]
    # sum up from columns '102.5' to '117.5'
    tetra = df.loc[:,'102.5':'117.5'].sum(axis=1)*5 # 5 degree bin width
    # divide each row by the bulk value
    rel_tetra = tetra/tetra['bulk']
    # subtract each row by the bulk values
    diff = df.sub(df.loc['bulk'],axis=1)
    # sum up from columns '102.5' to '117.5'
    tetra_diff = diff.loc[:,'102.5':'117.5'].sum(axis=1)*5 # 5 degree bin width
    # save new dataframe with tetra and tetra_frac
    processed_data[ff] = pd.DataFrame({'tetra':tetra,'rel_tetra':rel_tetra,'tetra_diff':tetra_diff})

#%%
colors = ['tab:red','#039dfc','black']
fig, ax = plt.subplots(figsize=(2.75,4))
for i,ff in enumerate(ffs):
    df = processed_data[ff]
    selected_data = [df.loc['y=20-24', 'rel_tetra'], df.loc['y=5-10', 'rel_tetra']]
    percent_inc = 100*(np.array(selected_data)-1)
    print(ff,percent_inc)
    ax.plot([0,1],selected_data,marker='s',label=ff,color=colors[i])
# set y-ticks as 'y=5-10 and y=20-24'
ax.set_xticks([0,1])
# ax.set_xticklabels(['y=20-24Å\n(undocked)','y=5-10Å\n(docked)'],fontsize=12)
ax.set_xticklabels(['undocked\npeptide\n(14<y<18Å)','docked\npeptide\n(y < 4Å)'],fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("relative water tetrahedrality\naround free peptide (local)",fontsize=12)
# plt.legend(loc='upper center' ,fontsize=12)
plt.ylim(bottom=1,top=1.11)
plt.text(0.5,1.052,'a03ws',fontsize=12,ha='center',color=colors[0],rotation=12)
plt.text(0.5,1.083,'a99SBdisp',fontsize=12,ha='center',color=colors[1],rotation=16)
plt.text(0.5,1.016,'CHARMM36m',fontsize=12,ha='center',color=colors[2],rotation=8)
# shift up
# plt.subplots_adjust(bottom=0.18)
plt.tight_layout()
plt.savefig('VLG_3body_rel_tetra.png',dpi=300)


# %%
colors = ['r','b','g']
fig, ax = plt.subplots(figsize=(2.5,6))
for i,ff in enumerate(ffs):
    df = processed_data[ff]
    selected_data = [df.loc['y=20-24', 'tetra'], df.loc['y=5-10', 'tetra']]
    ax.plot([0,1],selected_data,marker='s',label=ff,color=colors[i])
# set y-ticks as 'y=5-10 and y=20-24'
ax.set_xticks([0,1])
ax.set_xticklabels(['y=20-24Å\n(undocked)','y=5-10Å\n(docked)'],fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("VLG peptide's\nrelative water tetrahedrality",fontsize=18)
plt.legend(loc='center' ,fontsize=15)
# plt.ylim(bottom=1,top=1.13)

#%%
#%%
# plot tetra_diff for each forcefield
colors = ['r','b','g']
fig, ax = plt.subplots(figsize=(2.5,6))
for i,ff in enumerate(ffs):
    df = processed_data[ff]
    selected_data = [df.loc['y=20-24', 'tetra_diff'], df.loc['y=5-10', 'tetra_diff']]
    ax.plot([0,1],selected_data,marker='s',label=ff,color=colors[i])
# set y-ticks as 'y=5-10 and y=20-24'
ax.set_xticks([0,1])
ax.set_xticklabels(['y=20-24Å\n(undocked)','y=5-10Å\n(docked)'],fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel(r"$p_{tet}$ around peptide",fontsize=18)
plt.ylim(bottom=0, top=0.03)
plt.legend(loc='upper center' ,fontsize=15)
plt.show()
