#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

forcefields = ['a03ws','a99SBdisp','C36m']
types = ['monomer','dimer']
colors = {'a03ws':'red','a99SBdisp':'blue','C36m':'green'}
hatches = {'monomer':'','dimer':'//'}

all_data = pd.DataFrame()

for forcefield in forcefields:
    for type in types:
        data = pd.read_csv(f'{forcefield}_{type}.csv',index_col=0)
        if type == 'dimer':
            data = data / 2 # per monomer
        data['Forcefield'] = forcefield
        data['Type'] = type
        all_data = pd.concat([all_data, data])

fig,ax = plt.subplots()

# Create a mapping of forcefields to their position on the plot
base_positions = {'a03ws':1,'a99SBdisp':2,'C36m':3}

# Create violin plots
for i, forcefield in enumerate(forcefields):
    for j, type in enumerate(types):
        type_data = all_data[(all_data['Forcefield'] == forcefield) & (all_data['Type'] == type)]
        values = type_data.iloc[:,:-2].values.flatten()
        pos = base_positions[forcefield] + (-0.1 if type=='monomer' else 0.1)
        parts = ax.violinplot(values, [pos], points=100, widths=0.5, showmeans=False, showmedians=True)

        # for pc i