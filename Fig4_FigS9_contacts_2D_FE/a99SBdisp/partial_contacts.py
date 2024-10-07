#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from tqdm import tqdm

# Load the distances from a .npy file
distances = np.load('mindist.npy')  # Adjust the file path if needed

# Load the CSV with residue pairs
residue_pairs = pd.read_csv('mindist_pairs.csv')

# Dimensions based on your data
n_frames = distances.shape[0]
n_pairs = distances.shape[1]

# Calculating tot_residues using the quadratic formula
a = 1
b = -1
c = -2 * n_pairs
discriminant = b**2 - 4 * a * c
tot_residues = int((-b + math.sqrt(discriminant)) / (2 * a))
n_residues = int(tot_residues / 2)

# Initialize a 3D distance matrix of zeros (residues, residues, frames)
distance_matrix = np.zeros((tot_residues, tot_residues, n_frames))

# Populate the distance matrix using mapped indices
for pair_index, (i, j) in enumerate(zip(residue_pairs['residue1'], residue_pairs['residue2'])):
    distance_matrix[i, j, :] = distances[:, pair_index]
    distance_matrix[j, i, :] = distances[:, pair_index]  # Symmetric

# Convert distances to contacts based on a sigmoid
def sigmoid(x):
    k = 20  # Steepness adjustment for the specified range
    x0 = 0.5
    return -1 / (1 + np.exp(-k * (x - x0))) + 1

contact_matrix = sigmoid(distance_matrix)

#%%
# ignore intrachain contacts between residues within two residues of each other
mask = np.ones_like(contact_matrix[:,:,0], dtype=bool)
for chain_start in [0, n_residues]:  # This iterates over the start index of each chain
    for i in range(chain_start, chain_start + n_residues):
        for offset in range(-2, 3):
            j = i + offset
            if chain_start <= j < chain_start + n_residues:  # Ensure j is within the same chain
                mask[i, j] = False
# Apply the mask to all frames
mask = np.repeat(mask[:, :, np.newaxis], n_frames, axis=2)
contact_matrix *= mask

seq = 'DNIKHVLGGGSVQIVYKPV'*2

# Function to plot contact map for a specific frame
def plot_contact_map(frame_index):
    plt.figure(figsize=(10, 8))
    plt.imshow(contact_matrix[:, :, frame_index], cmap='coolwarm_r', interpolation='none')
    plt.title(f'Contact Map for Frame {frame_index}')
    plt.xticks(np.arange(tot_residues), list(seq))
    plt.yticks(np.arange(tot_residues), list(seq))
    # draw lines to separate chains
    plt.axvline(x=n_residues-0.5, color='black', linestyle='-')
    plt.axhline(y=n_residues-0.5, color='black', linestyle='-')
    plt.colorbar()
    plt.show()

# Example of how to use the function
# plot_contact_map(1000)  # Change the frame index as needed

#%%
# create mask of inter PHF6 contacts
# Define the ranges for chain A and chain B
SPHF6a_indices = np.arange(11, 16)  # 10 to 16 inclusive
SPHF6b_indices = np.arange(30, 35)  # 29 to 35 inclusive

highlight_mask = np.zeros((tot_residues, tot_residues), dtype=bool)

def set_mask_safe(mask, i, j):
    if 0 <= i < mask.shape[0] and 0 <= j < mask.shape[1]:
        mask[i, j] = True

# Set diagonal and neighboring cells
for a, b in zip(SPHF6a_indices, SPHF6b_indices):
    # Set the diagonal element
    set_mask_safe(highlight_mask, a, b)
    set_mask_safe(highlight_mask, b, a)  # Mirror for the opposite diagonal
    # Set neighbors
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            set_mask_safe(highlight_mask, a + di, b + dj)
            set_mask_safe(highlight_mask, b + di, a + dj)

# Use this mask to visualize the interactions, for example, for frame 0
plt.figure(figsize=(10, 8))
plt.imshow(contact_matrix[:, :, 1000], cmap='coolwarm_r', interpolation='none', vmin=0, vmax=1)
# plt.colorbar()
plt.imshow(highlight_mask, cmap='Greys', alpha=0.5, interpolation='none')  # Overlay the highlight mask with transparency
# set ticks
plt.xticks(np.arange(tot_residues), list(seq))
plt.yticks(np.arange(tot_residues), list(seq))
# draw lines to separate chains
plt.axvline(x=n_residues-0.5, color='black', linestyle='-')
plt.axhline(y=n_residues-0.5, color='black', linestyle='-')
plt.show()

#%%
# count up intrachain contacts per frame
intraA_contacts = np.zeros(n_frames)
intraB_contacts = np.zeros(n_frames)

chain_a_indices = np.arange(0, n_residues)
chain_b_indices = np.arange(n_residues, tot_residues)

# create a mask for the intrachain contacts
intraA_mask = np.zeros((tot_residues, tot_residues), dtype=bool)
intraB_mask = np.zeros((tot_residues, tot_residues), dtype=bool)

for i in chain_a_indices:
    for j in chain_a_indices:
        intraA_mask[i, j] = True
        intraA_mask[j, i] = True
for i in chain_b_indices:
    for j in chain_b_indices:
        intraB_mask[i, j] = True
        intraB_mask[j, i] = True
intrachain_mask = np.logical_or(intraA_mask, intraB_mask)

# for each frame, count the number of intrachain contacts
for frame in tqdm(range(n_frames)):
    intraA_contacts[frame] = np.sum(contact_matrix[intraA_mask, frame])
    intraB_contacts[frame] = np.sum(contact_matrix[intraB_mask, frame])

# sum intra-A and intra-B contacts per frame
intra_contacts = intraA_contacts + intraB_contacts

# histogram of total contacts
plt.figure()
plt.hist(intra_contacts, bins=50)
plt.xlabel('Total Intra-chain Contacts')
plt.ylabel('Frequency')

#
# plot the time evolution of the total contacts
plt.figure()
plt.plot(intra_contacts,lw=0.1)
plt.xlabel('Frame')
plt.ylabel('Intra-chain Contacts')

#%%
# sum the contents in the highlight mask for each frame
highlight_contacts = np.zeros(n_frames)
interchain_SPHF6_contacts = np.sum(contact_matrix[highlight_mask, :], axis=0)

# histogram of interchain contacts
plt.figure()
plt.hist(interchain_SPHF6_contacts, bins=50)
plt.xlabel('Interchain SPHF6 Aligned Contacts')
plt.ylabel('Frequency')

# plot the time evolution of the interchain contacts
plt.figure()
plt.plot(interchain_SPHF6_contacts, lw=0.1)

#%%
# create a df
df = pd.DataFrame({'Intrachain': intra_contacts, 'Interchain SPHF6 Aligned': interchain_SPHF6_contacts})
df.to_csv('contacts_dimer_a99SBdisp.csv', index=False)

# %%
