#%%
import numpy as np
import pandas as pd
from openmm import *
from openmm.app import *
from openmm.unit import *
import MDAnalysis as mda
from tqdm import tqdm

protein_file = 'SLS2_processed.pdb'

pdb = PDBFile(protein_file)
top = GromacsTopFile('topol.top', periodicBoxVectors=pdb.topology.getPeriodicBoxVectors())

# Store box size
box_vectors = pdb.topology.getPeriodicBoxVectors()
print("Box Vectors:", box_vectors)

# System Configuration
# nonbondedMethod = PME # not using PME! just ignoring interactions longer than 1.2 nm
nonbondedCutoff = 1.2*nanometers # I'm using a switching distance of 1.0 nm
constraints = HBonds
system = top.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=nonbondedCutoff, constraints=constraints)

#%%
seed = set([atom.index for atom in pdb.topology.atoms() if atom.residue.chain.id in ('A','B','C','D','E')])
solvent = set([atom.index for atom in pdb.topology.atoms() if atom.index not in seed])
print(f'{len(seed)} atoms in seed')
print(f'{len(solvent)} atoms in solvent')

# %%
# check what force objects there are
print('Checking the relevant force objects in the system:')
for force in system.getForces():
    if isinstance(force, NonbondedForce):
        print("\nThere is a NonbondedForce")
        print(f'{force.getNumParticles()} particles')
        print("Some NonbondedForce parameters:")
        for i in range(5):
            charge, sigma, epsilon = force.getParticleParameters(i)
            print(f"Particle {i}: charge={charge}, sigma={sigma}, epsilon={epsilon}")
# check the parameters of the CustomNonbondedForce
    if isinstance(force, CustomNonbondedForce):
        print("\nThere is a CustomNonbondedForce")
        print(f'{force.getNumParticles()} particles')
        print(f'Energy Function: {force.getEnergyFunction()}')
        print("Some CustomNonbondedForce parameters:")
        for i in range(5):
            print(f'Particle {i}: {force.getParticleParameters(i)}')

#%%
for force in system.getForces():
    if isinstance(force, NonbondedForce):
        force.setForceGroup(0)
        force.setSwitchingDistance(1.0*nanometers)
        force.addGlobalParameter("seed_scale", 1)
        force.addGlobalParameter("solvent_scale", 1)
        for i in range(force.getNumParticles()):
            charge, sigma, epsilon = force.getParticleParameters(i)
            # Set the parameters to be 0 when the corresponding parameter is 0,
            # and to have their normal values when it is 1.
            param = "seed_scale" if i in seed else "solvent_scale"
            # print(i, param)
            force.setParticleParameters(i, 0, 0, 0)
            force.addParticleParameterOffset(param, i, charge, sigma, epsilon)
        # print(force.getNumExceptions())
        # print(force.getExceptionParameters(0))
        for i in range(force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
            force.setExceptionParameters(i, p1, p2, 0, 0, 0)
    elif isinstance(force, CustomNonbondedForce):
        # print("This is a CustomNonbondedForce")
        force.setForceGroup(1)
        force.setSwitchingDistance(1.0*nanometers)
        force.addInteractionGroup(seed, solvent)
    else:
        # print the name of the force object
        # print(force)
        force.setForceGroup(2)

#%%
# create a context (prerequisite for evaluating energies)
integrator = VerletIntegrator(0.001*picosecond)
context = Context(system, integrator)

# evaluate coulomb energy between chains
def coulomb_energy(seed_scale, solvent_scale):
    context.setParameter("seed_scale", seed_scale)
    context.setParameter("solvent_scale", solvent_scale)
    return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

# load trajectory
u = mda.Universe(protein_file,'traj.dcd')
# get the number of frames in the trajectory
n_frames = len(u.trajectory)
print(f'There are {n_frames} frames in the trajectory')
coul_energies = []
lj_energies = []

for frame in tqdm(np.arange(0, n_frames, 10)):
    u.trajectory[frame]
    # set the positions of the atoms in the context
    context.setPositions(u.atoms.positions / 10 * nanometers)
    
    total_coulomb = coulomb_energy(1, 1)
    seed_coulomb = coulomb_energy(1, 0)
    solvent_coulomb = coulomb_energy(0, 1)
    # print(f"Total coulombic potential energy:\n{total_coulomb}")
    # print(f"Seed coulombic potential energy:\n{seed_coulomb}")
    # print(f"Solvent coulombic potential energy:\n{solvent_coulomb}")
    coul_energy = total_coulomb - seed_coulomb - solvent_coulomb
    coul_energies.append(coul_energy.value_in_unit(kilojoules_per_mole))
    # print(f"Coulombic potential energy:\n{coul_energy.value_in_unit(kilojoules_per_mole):.0f} kJ/mol")
    
    lj_energy = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
    lj_energies.append(lj_energy.value_in_unit(kilojoules_per_mole))
    # print(f"LJ potential energy: {lj_energy.value_in_unit(kilojoules_per_mole):.0f} kJ/mol")

#%%
# save the energies to a file
energies = pd.DataFrame({'coulomb': coul_energies, 'lj': lj_energies})
energies.to_csv('seed_water_potential.csv', index=False)

# %%

# histogram of the coulombic potential energy
import matplotlib.pyplot as plt
plt.hist(coul_energies, bins=100)
plt.xlabel('Coulombic Potential Energy (kJ/mol)')
plt.ylabel('Frequency')
plt.title('Histogram of Coulombic Potential Energy')

# histogram of the LJ potential energy
plt.hist(lj_energies, bins=100)
plt.xlabel('LJ Potential Energy (kJ/mol)')
plt.ylabel('Frequency')
plt.title('Histogram of LJ Potential Energy')
plt.show()

# add them together and histogram that
total_energies = np.array(coul_energies) + np.array(lj_energies)
plt.hist(total_energies, bins=100)
plt.xlabel('Total Potential Energy (kJ/mol)')
plt.ylabel('Frequency')
plt.title('Histogram of Total Potential Energy')
plt.show()
#%%