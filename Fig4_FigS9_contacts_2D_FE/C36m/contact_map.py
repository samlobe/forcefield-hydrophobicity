#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

# Load the distances from a .npy file
distances = np.load('mindist.npy')  # Adjust the file path if needed

# Load the CSV with residue pairs
residue_pairs = pd.read_csv('mindist_pairs.csv')

# Define a mapping from original residue indices to 0-based indices
# Because Pritam's simulation had ACE and NME caps (C36m)
def map_residue_index(original_index):
    if original_index <= 19:
        return original_index - 1  # Maps 1-19 to 0-18
    else:
        return original_index - 3  # Maps 21-40 to 19-38

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
    mapped_i = map_residue_index(i)
    mapped_j = map_residue_index(j)
    distance_matrix[mapped_i, mapped_j, :] = distances[:, pair_index]
    distance_matrix[mapped_j, mapped_i, :] = distances[:, pair_index]  # Symmetric

# Convert distances to contacts based on the cutoff
cutoff = 0.65
contact_matrix = distance_matrix <= cutoff

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
contact_matrix = np.logical_and(contact_matrix, mask)

seq = 'DNIKHVLGGGSVQIVYKPV'*2

# Function to plot contact map for a specific frame
def plot_contact_map(frame_index):
    plt.figure(figsize=(10, 8))
    plt.imshow(contact_matrix[:, :, frame_index], cmap='Greys', interpolation='none')
    plt.title(f'Contact Map for Frame {frame_index}')
    plt.xticks(np.arange(tot_residues), list(seq))
    plt.yticks(np.arange(tot_residues), list(seq))
    # draw lines to separate chains
    plt.axvline(x=n_residues-0.5, color='red', linestyle='-')
    plt.axhline(y=n_residues-0.5, color='red', linestyle='-')
    plt.colorbar()
    plt.show()

# Example of how to use the function
plot_contact_map(1000)  # Change the frame index as needed
