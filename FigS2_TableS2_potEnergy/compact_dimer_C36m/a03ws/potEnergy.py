#%%
import numpy as np
import pandas as pd
from openmm import *
from openmm.app import *
from openmm.unit import *

protein_file = 'compact_dimer_processed.pdb'

pdb = PDBFile(protein_file)
top = GromacsTopFile('topol.top', periodicBoxVectors=pdb.topology.getPeriodicBoxVectors())

# Store box size
box_vectors = pdb.topology.getPeriodicBoxVectors()
print("Box Vectors:", box_vectors)

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
constraints = HBonds
system = top.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=nonbondedCutoff, constraints=constraints)
# system = top.createSystem(nonbondedMethod=NoCutoff, constraints=constraints)
# system = top.createSystem()

#%%
chainA = set([atom.index for atom in pdb.topology.atoms() if atom.residue.chain.id == 'A'])
chainB = set([atom.index for atom in pdb.topology.atoms() if atom.residue.chain.id == 'B'])
print(f'{len(chainA)} atoms in chainA')
print(f'{len(chainB)} atoms in chainB')

# %%
# check what force objects there are
for force in system.getForces():
    print(force)
    if isinstance(force, NonbondedForce):
        print(force.getNumParticles())
        print(force.getName())
        for i in range(force.getNumParticles()):
            charge, sigma, epsilon = force.getParticleParameters(i)
            print(f"Particle {i}: charge={charge}, sigma={sigma}, epsilon={epsilon}")
#%%
# check the parameters of the CustomNonbondedForce
for force in system.getForces():
    if isinstance(force, CustomNonbondedForce):
        print(force.getNumInteractionGroups())
        # print(force.getInteractionGroupParameters(0))
        print(force.getNumParticles())
        print(force.getEnergyFunction())
        for i in range(force.getNumParticles()):
            print(force.getParticleParameters(i))

#%%
for force in system.getForces():
    if isinstance(force, NonbondedForce):
        force.setForceGroup(0)
        force.addGlobalParameter("chainA_scale", 1)
        force.addGlobalParameter("chainB_scale", 1)
        print(force.getNumParticles())
        print(force.getParticleParameters(0))
        for i in range(force.getNumParticles()):
            charge, sigma, epsilon = force.getParticleParameters(i)
            # Set the parameters to be 0 when the corresponding parameter is 0,
            # and to have their normal values when it is 1.
            param = "chainA_scale" if i in chainA else "chainB_scale"
            # print(i, param)
            force.setParticleParameters(i, 0, 0, 0)
            force.addParticleParameterOffset(param, i, charge, sigma, epsilon)
        # print(force.getNumExceptions())
        # print(force.getExceptionParameters(0))
        for i in range(force.getNumExceptions()):
            p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
            force.setExceptionParameters(i, p1, p2, 0, 0, 0)
    elif isinstance(force, CustomNonbondedForce):
        print("This is a CustomNonbondedForce")
        force.setForceGroup(1)
        force.addInteractionGroup(chainA, chainB)
    else:
        # print the name of the force object
        print(force)
        force.setForceGroup(2)

#%%
# create a context (prerequisite for evaluating energies)
integrator = VerletIntegrator(0.001*picosecond)
context = Context(system, integrator)
context.setPositions(pdb.positions)

#%%
# evaluate potential energy between chains
def potential_energy(chainA_scale, chainB_scale):
    context.setParameter("chainA_scale", chainA_scale)
    context.setParameter("chainB_scale", chainB_scale)
    return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

total_potential = potential_energy(1, 1)
chainA_potential = potential_energy(1, 0)
chainB_potential = potential_energy(0, 1)
print(f"Total  potential energy:\n{total_potential}")
print(f"Chain A potential energy:\n{chainA_potential}")
print(f"Chain B potential energy:\n{chainB_potential}")
print(f"Chain Interaction potential energy (LJ + Coulomb):\n{total_potential - chainA_potential - chainB_potential}")


#%%
print(f"LJ interaction: {context.getState(getEnergy=True, groups={1}).getPotentialEnergy()}")

# %%
