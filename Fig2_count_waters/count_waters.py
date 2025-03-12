# example script for counting hydration waters
import MDAnalysis as mda
import pandas as pd
from tqdm import tqdm
u = mda.Universe('SLS2_processed.pdb','traj.dcd')
num_waters = []

waters = 'resname SOL and name O*'
protein_residues = 'resid 2-18'
selection = f'({waters}) and around 4.25 ({protein_residues})'
print(len(u.select_atoms(protein_residues)))

for ts in tqdm(u.trajectory):
    shell_waters = u.select_atoms(selection)
    num_waters.append(len(shell_waters))

df = pd.DataFrame(num_waters,columns=['Hydration Waters'])
df.to_csv('num_waters.csv')

