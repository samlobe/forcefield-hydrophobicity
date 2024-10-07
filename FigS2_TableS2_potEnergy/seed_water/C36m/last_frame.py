import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

# Load the trajectory and topology
u = mda.Universe('SLS2_processed.pdb', 'traj.dcd')

# Select all the atoms
atoms = u.select_atoms('all')

# Access the last frame
u.trajectory[-1]

# Write the last frame to a PDB file
with mda.Writer('last_frame.pdb', atoms.n_atoms) as W:
    W.write(atoms)

