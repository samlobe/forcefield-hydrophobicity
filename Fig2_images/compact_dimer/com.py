from Bio.PDB import PDBParser, NeighborSearch
import numpy as np

# Load the PDB file
pdb_file_path = "cluster1.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", pdb_file_path)

# Function to calculate the center of mass of a chain
def calculate_center_of_mass(chain):
    mass = 0.0
    center_of_mass = np.zeros(3)
    for residue in chain:
        for atom in residue:
            atom_mass = atom.mass
            center_of_mass += atom_mass * atom.coord
            mass += atom_mass
    center_of_mass /= mass
    return center_of_mass

# Extract chains A and B
chain_A = structure[0]["A"]
chain_B = structure[0]["B"]

# Calculate centers of mass
com_A = calculate_center_of_mass(chain_A)
com_B = calculate_center_of_mass(chain_B)

# Calculate the distance between the centers of mass
distance = np.linalg.norm(com_A - com_B)
print(distance)

