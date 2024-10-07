#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

forcefields = ['a03ws','a99SBdisp','C36m']
# from local_densities.py:
n_oxygens_bulk = [5.913931034482759, 6.008331093791991, 6.098301698301698]
n_hydrogens_bulk = [11.942344827586206, 11.995968825584521, 12.205594405594406]

avg_oxygens = []
avg_hydrogens = []
for ff in forcefields:
    df = pd.read_csv(f'local_density_{ff}.csv',index_col=0)
    avg_oxygens_bySlice = []
    avg_hydrogens_bySlice = []
    for ang_slice in np.arange(36,25,-1):
        avg_oxygens_bySlice.append(np.mean(df[f'n_oxygens_{ang_slice}']))
        avg_hydrogens_bySlice.append(np.mean(df[f'n_hydrogens_{ang_slice}']))
    avg_oxygens.append(avg_oxygens_bySlice)
    avg_hydrogens.append(avg_hydrogens_bySlice)

ys = np.arange(1,12) # angstroms above the surface
df_oxygens = pd.DataFrame(avg_oxygens,columns=ys,index=forcefields)
df_hydrogens = pd.DataFrame(avg_hydrogens,columns=ys,index=forcefields)

#%%
# loop through the forcefields and divide each row by the bulk value
for ff in forcefields:
    df_oxygens.loc[ff] = df_oxygens.loc[ff] / n_oxygens_bulk[forcefields.index(ff)]
    df_hydrogens.loc[ff] = df_hydrogens.loc[ff] / n_hydrogens_bulk[forcefields.index(ff)]

#%%

# fig, ax = plt.subplots(1,2,figsize=(10,5))
fig, ax = plt.subplots(figsize=(6,4))
colors = ['tab:red','#039dfc','black']
df_oxygens.T.plot(color=colors, ax=ax)
# df_hydrogens.T.plot(ax=ax[1],title='Hydrogen Density',xlabel='Distance from Surface (Å)',ylabel='Density (Å^-2)')

# dotted horizontal line at 1
plt.axhline(1,linestyle='--',color='gray')
# increase fontsize and ticksize
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('distance from surface (Å)',fontsize=14)
plt.ylabel('relative oxygen density',fontsize=14)
plt.legend(fontsize=14)
plt.subplots_adjust(bottom=0.15,left=0.15)
plt.savefig('oxygen_density_profile_upright_seed.png',dpi=300)
plt.show()
