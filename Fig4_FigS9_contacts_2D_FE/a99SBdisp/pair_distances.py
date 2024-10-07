#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md
from time import time

# load the trajectory
molecule = 'HP2_dimer'
traj = md.load(f'../{molecule}.xtc', top=f'../{molecule}.pdb')

#%%
# # calculate the distances between all pairs of alpha carbons

# # get the indices of the alpha carbons
# alpha_carbons = [atom.index for atom in traj.topology.atoms if atom.name == 'CA']
# # get all the pairs of alpha carbons
# pairs = np.array([[i, j] for i in alpha_carbons for j in alpha_carbons if i < j])
# # get chainID and residue number for each alpha carbon
# chainID = [traj.topology.atom(i).residue.chain.index for i in alpha_carbons]
# residue_number = [traj.topology.atom(i).residue.index for i in alpha_carbons]

# # convert chainID from 0 and 1 to A and B
# chainID = ['A' if i == 0 else 'B' for i in chainID]

# labels = [f'{chainID[i]}{residue_number[i]}' for i in range(len(alpha_carbons))]

# #%%

# # calculate the distances between all pairs of alpha carbons
# tik = time()
# distances = md.compute_distances(traj,pairs)
# tok = time()
# print(f'Done in {round(tok-tik)} seconds')

# # calculate contacts 
# tik = time()
# mindistances = md.compute_contacts(traj, scheme='closest-heavy', periodic=False)
# tok = time()
# print(f'Done in {round(tok-tik)} seconds')

#%%
# get resname of each residue
alpha_carbons = [atom.index for atom in traj.topology.atoms if atom.name == 'CA']
resname = [traj.topology.atom(i).residue.name for i in alpha_carbons]
## get one letter code for each residue
one_letter_code = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
resnames = [one_letter_code[resname[i]] for i in range(len(resname))]

chainIDs = [traj.topology.atom(i).residue.chain.index for i in alpha_carbons]
resids = [traj.topology.atom(i).residue.index for i in alpha_carbons]

#%%
chainA_resids = np.arange(0,19)
chainB_resids = np.arange(19,38)
my_residues = np.concatenate((chainA_resids,chainB_resids))

# iterate over my residue to get all pairs
my_pairs = np.array([[i, j] for i in my_residues for j in my_residues if i < j])

# calculate contacts
tik = time()
my_contacts = md.compute_contacts(traj, my_pairs, scheme='closest-heavy', periodic=False)
tok = time()
print(f'Done in {round(tok-tik)} seconds') # done in 350 seconds

#%%
# save my_contacts to a .npy file
np.save(f'mindist.npy', my_contacts[0])

# output my_contacts[1] to a .csv file
df = pd.DataFrame(my_contacts[1], columns=['residue1','residue2'])
df.to_csv('mindist_pairs.csv', index=False)


