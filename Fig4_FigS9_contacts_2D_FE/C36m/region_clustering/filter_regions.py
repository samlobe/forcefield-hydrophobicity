#%%
import pandas as pd
import MDAnalysis as mda

df = pd.read_csv('../region_masks.csv')
#%%

# regions = ['i','ii','iii','iv','v']
regions = ['iv']
u = mda.Universe('../../../HP2_dimer.pdb','../../../HP2_dimer.xtc')
selection = u.select_atoms('all')

#%%
for region in regions:
    with mda.Writer(f'C36m_region{region}.xtc', selection.n_atoms) as W:
        mask = df[region].values
        for ts in u.trajectory[mask]:
            W.write(selection)

#%%