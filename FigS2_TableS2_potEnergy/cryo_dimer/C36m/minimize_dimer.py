#%%
import numpy as np
import pandas as pd
from openmm import *
from openmm.app import *
from openmm.unit import *

protein_file = 'cryo_dimer_processed.pdb'

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
# system = top.createSystem(nonbondedMethod=CutoffNonPeriodic, constraints=constraints)

#%%
platform = Platform.getPlatformByName('OpenCL')
friction = 2/picosecond
temperature = 300*kelvin
dt = 0.002*picoseconds
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

#%%
# minimize the energy
print('Minimizing...')
simulation.minimizeEnergy()
print('Done minimizing')
#%%
# output the minimized structure
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimized.pdb', 'w'))