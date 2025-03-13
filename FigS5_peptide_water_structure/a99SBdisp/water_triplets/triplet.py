#%%
import numpy as np
import MDAnalysis as mda
import water_properties as wp
from tqdm import tqdm
import sys
import os
import argparse
import shutil
from glob import glob

# Check for the compiled fortran code in the current directory
matching_files = glob("waterlib*.so") + glob("waterlib*.pyd")
if not matching_files:
    raise FileNotFoundError("Could not find a compiled Fortran code in the current directory. (Looking for these patterns: waterlib*.so or waterlib*.pyd) \n"
                            "Please compile the Fortran code first. Example compilation using f2py (included with anaconda):`f2py -c -m waterlib waterlib.f90`\n"
                            "Please consult the tutorial (on GitHub) and leave an issue on the GitHub page if you have trouble compiling waterlib.f90's fortran code.\n")

# Setting up argparse
parser = argparse.ArgumentParser(description='Compile water triplet angles around the residue/custom group in your protein.\nExample usage: `python triplet.py myProtein_processed.gro traj.dcd --selection _____')

# Add arguments
parser.add_argument('protein',type=str,help="Processed protein structure file (e.g. 'myProtein_processed.gro')")
parser.add_argument('trajectory', type=str, help="Trajectory file name (e.g. 'traj.dcd')")
parser.add_argument('-t', '--time', type=float, default=5.0, help='Last x ns. Default is 5 ns.')
parser.add_argument('--selection', type=str, help='MDAnalysis selection string to select region of waters.')
parser.add_argument('-o', '--output',type=str, help='Outputted text file.')

args = parser.parse_args()

# Assign the arguments
protein_processed = args.protein
last_x_ns = args.time
selection = args.selection
output_file = args.output

# Ensure that the protein file exists
if not os.path.exists(protein_processed):
    print(f"Error: Can't find {protein_processed}")
    sys.exit(1)

# Ensure that the trajectory file exists
if not os.path.exists(args.trajectory):
    print(f"Error: Can't find {args.trajectory}")
    sys.exit(1)

structure_path = protein_processed  # Looking at structure file (usually in parent directory)
traj_path = args.trajectory

### LOAD THE MD TRAJECTORY
u = mda.Universe(structure_path,traj_path)
total_frames = len(u.trajectory)
timestep = 2 # ps
# convert ns to ps and then divide by timestep to get number of frames
frames_to_load = int((last_x_ns * 1e3) / timestep)

# Continue setting up analysis
waters = f'resname SOL and name O*' # looking at just water oxygens
BoxDims = u.dimensions[:3] # Å; needed for periodic boundary conditions
lowCut = 0
highCut = 3.5 # Å; cutoff distance to establish the neighbors of each hydration water

# Create a list of lists of 3-body angles.
# One sublist per configuration (i.e. per frame of trajectory)
# Each sublist contains the 3-body angle for the waters in the first shell around your group of interest.
angles_list = [[] for i in range(len(u.trajectory))]

# if using --selection, come up with an apppropriate name for the output files
if args.selection:
    my_group = args.selection
    group_name = my_group.replace(" ","_")

# Create a checkpoint file to save progress (every 1000 frames)
if args.selection:
    checkpoint_filename = f'checkpoint_{args.output}.txt'


# Load from checkpoint if exists
start_frame = 0
if os.path.exists(checkpoint_filename):
    with open(checkpoint_filename, 'r') as file:
        for line in file:
            frame, angles = line.split(":")
            angles_str = angles.strip()[1:-1]
            if angles_str:
                angles = [float(angle) for angle in angles_str.split(',')]
            else: angles = []
            angles_list[int(frame)] = angles
        start_frame = int(frame) + 1

### CALCULATE THE WATER TRIPLET ANGLES
start_frame_for_last_x_ns = max(0,total_frames - frames_to_load) # ensure it's not negative
for i,ts in enumerate(tqdm(u.trajectory[start_frame_for_last_x_ns+start_frame:])):
    # SELECT HYDRATION WATERS
    shell_waters = u.select_atoms(f'({waters}) and ({selection})') # 4.25Å is ~2 water layers from the residue
    subPos = shell_waters.positions # hydration water positions
    Pos = u.select_atoms(waters).positions # all water positions
    
    # MEAUSURE TRIPLET ANGLES
    triplet_angs, _ = wp.getTripletAngs(subPos,Pos,BoxDims,lowCut,highCut)
    #print(three_body_angs) # for debugging
    angles_list[i+start_frame] = list(np.around(triplet_angs,1))
    
    # Save checkpoint every 1000 frames
    if (i+1) % 1000 == 0:
        with open(checkpoint_filename, 'w') as txtfile:
            for j in range(i+start_frame+1):
                txtfile.write(f"{j}:{str(angles_list[j])}\n")

### SAVE FINAL RESULT
if args.selection:
    output_filename = f'{args.output}_angles.txt'

with open(output_filename,'w') as txtfile:
    # each line contains the 3-body angles for one configuration (i.e. frame of trajectory)
    # in each line, there are triplet angles for each of the hydration waters
    for line in angles_list:
        txtfile.write(str(line)[1:-1].replace(",","").replace("\n","  ")+'\n') # string formatting

# Delete checkpoint file
if os.path.exists(checkpoint_filename):
    os.remove(checkpoint_filename)

# create directory 'angles' if it doesn't exist
if not os.path.exists('angles'):
    try:
        os.makedirs('angles')
    except FileExistsError:
        pass

# move output file to angles
shutil.move(output_filename,f'angles/{output_filename}')

# %%
