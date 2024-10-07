#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ffs = ['a03ws','a99SBdisp','C36m']
potentials = []
for ff in ffs:
    data = pd.read_csv(f'{ff}/seed_water_potential.csv')
    
    if ff == 'C36m':
        # sum over each row
        potential = data.sum(axis=1).values
        potentials.append(potential)
    else:
        # sum over each row
        potential = data.sum(axis=1).values
        potentials.append(potential)

#%%
# turn potentials into a df
potentials_df = pd.DataFrame(potentials).T
potentials_df.columns = ffs

#%%
fig, ax = plt.subplots(figsize=(6,4))
colors = ['tab:red', '#039dfc', [0.2,0.2,0.2]]
sns.violinplot(data=potentials_df, palette=colors, linecolor='gray',inner='box')
# plt.title('Potential Energy (kJ/mol)\nbetween peptides and water')
plt.xlabel('force field',fontsize=12)
plt.ylabel('Potential Energy (kJ/mol)\nbetween seed and water',fontsize=12)
plt.xticks(fontsize=12)
# rename xticks to force field names
plt.xticks(range(3),['a03ws','a99SBdisp','CHARMM36m'])

plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig('seed_water_PE.png',dpi=300)
plt.show()

