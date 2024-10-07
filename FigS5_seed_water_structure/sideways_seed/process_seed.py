#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

# initalize histogramming of the water triplet distribution
bin_width = 5 # degrees
min = 40
max = 180
bins = np.arange(min,max+bin_width,bin_width)
def histo_line(data):
    histo_height, bin_edges = np.histogram(data, bins=bins, density=True)
    bin_middle = np.diff(bin_edges)/2
    bin_middle = bin_middle + bin_edges[:len(bin_middle)]
    return bin_middle, histo_height

# function to read the angles from the .txt files created by triplets.py
def read_angles(filepath):
    with open(filepath,'r') as f:
        angles = []
        for i,line in enumerate(tqdm(f)):
            if not line:
                print(f'Line {i} is empty.')
                angles.append([])
                continue
            else:
                angles.append([float(x) for x in line.split()])
                
    all_angles = np.array([item for sublist in angles for item in sublist])
    # avg_measures = len(all_angles) / len(angles) # avg number of angle measurements per frame
    return all_angles

#%% 
ffs = ['a03ws','a99SBdisp','C36m']
slices = np.arange(1,7+1) # 6 slices, each is 1 A thick, and the last one is bulk
distances = [f'y={1.0*(i-1):.0f}-{1.0*i:.0f}Å' for i in slices]

cmap = plt.get_cmap('coolwarm_r')
colors = cmap(np.linspace(0, 1, len(slices)))

# ff = ffs[2]
for ff in ffs:
    histo_heights = []

    for i,slice in enumerate(slices):
        distance = distances[i]
        angles = read_angles(f'{ff}/water_triplets/angles/slice_{slice}_angles.txt')
        # flatten the list of lists
        bin_middle, histo_height = histo_line(angles)
        histo_heights.append(histo_height)

    plt.figure()
    for i in range(len(slices[:-1])):
        plt.plot(bin_middle, histo_heights[i],label=distances[i],color=colors[i])
    # plot bulk
    plt.plot(bin_middle, histo_heights[-1],color='black',lw=2)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Probability Density')
    plt.title('3-body angle distribution\nat slices above PHF6 in seed')
    plt.ylim(bottom=0)
    plt.legend()

    # subtract the bulk water distribution from the seed water distribution
    histo_heights_diff = []
    for i in range(len(slices[:-1])):
        histo_heights_diff.append(histo_heights[i] - histo_heights[-1])

    plt.figure()
    for i in range(len(slices[:-1])):
        plt.plot(bin_middle, histo_heights_diff[i],label=distances[i],color=colors[i])
    # plot horizontal line at x=0
    plt.axhline(0, color='black',linewidth=2)
    plt.xlabel('Angle (degrees)')
    plt.ylabel('Deviation from bulk water')
    plt.title(f'{ff} 3-body angle distribution\nat slices above L301')
    plt.ylim([-.0020,.0020])
    plt.legend(loc='upper right',ncol=2,fontsize=8)
    plt.subplots_adjust(left=0.20)
    plt.savefig(f'sideways_seed_3body_{ff}.png',dpi=300)

#%%
