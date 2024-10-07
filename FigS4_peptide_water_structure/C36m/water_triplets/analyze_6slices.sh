#!/bin/bash

# Define the common parameters
pdb_file="../SLS2_processed.pdb"
dcd_file="../traj.dcd"
selection_base="prop x > 7 and prop x < 23 and prop y > 11 and prop y < 19 and prop z >"
time="-t 10"

# Run the commands with incrementing z coordinates
for i in {0..5}
do
  z_min=$(echo "43.5 + $i" | bc)
  z_max=$(echo "44.5 + $i" | bc)
  selection="$selection_base $z_min and prop z < $z_max"
  output="slice_$((i + 1))"
  
  python triplet.py $pdb_file $dcd_file --selection "$selection" -o $output $time
done

