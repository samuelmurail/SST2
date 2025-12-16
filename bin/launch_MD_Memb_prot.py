import argparse
import numpy as np
import time
import os
import sys
import logging
import pandas as pd

import openmm
from openmm import app
from openmm.app import PDBFile, PDBxFile, ForceField, Simulation
from openmm import LangevinMiddleIntegrator, unit, Platform

sys.path.insert(0, '/home/murail/Documents/Code/SST2/src/')

from SST2.st import ST, run_st
import SST2.tools as tools
import pdb_numpy


# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Add sys.sdout as handler
logger.addHandler(logging.StreamHandler(sys.stdout))



###########################
    ### READING FILES ###
###########################
logger.info(f"Read input files")
pdb_file = "../test_alaa/openmm/pth1r-pth-Gs-largeBox.rst7"
top_file =  "../test_alaa/openmm/pth1r-pth-Gs-largeBox.parm7"

inpcrd = app.AmberInpcrdFile(pdb_file)
prmtop = app.AmberPrmtopFile(top_file, periodicBoxVectors=inpcrd.boxVectors)



#################################################
### SIMULATION PARAMETERS AND SYSTEM CREATION ###
#################################################
logger.info(f"Define simulation parameters")
hmr = 3.0
eq_time_expl = 10.0

dt = 1 * unit.femtosecond
temperature = 300.0 * unit.kelvin
friction = 1.0 / unit.picoseconds
hydrogenMass = hmr * unit.amu
rigidWater = True
ewaldErrorTolerance = 0.0005
nsteps = int(np.ceil(eq_time_expl * unit.nanoseconds / dt))

logger.info(f"Set up system")
system = prmtop.createSystem(nonbondedMethod=app.PME,
                 nonbondedCutoff=0.9*unit.nanometers,
                 constraints=app.HBonds,
                 rigidWater=rigidWater,
                 hydrogenMass=hydrogenMass,
                 ewaldErrorTolerance=ewaldErrorTolerance)

logger.info(f"Set up integrator and platform")
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'single'}

logger.info(f"Add MonteCarloMembraneBarostat")
barostat = openmm.MonteCarloMembraneBarostat(
    1*unit.bar,
    0*unit.bar*unit.nanometer,
    temperature,
    openmm.MonteCarloMembraneBarostat.XYIsotropic,
    openmm.MonteCarloMembraneBarostat.ZFree,
    100)
system.addForce(barostat)

####################################
### CA ATOM POSITION RESTRAINTS ###
####################################

CA_indices = [int(i.index) for i in prmtop.topology.atoms() if i.name in ['CA']]

logger.info(f"- Add position restraints")

restraint = tools.add_pos_restr(system, CA_indices, inpcrd, k_rest=1000)

simulation = Simulation(
        prmtop.topology,
        system, 
        integrator, 
        platform, 
        platformProperties)

simulation.context.setPositions(inpcrd.positions)
simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

###############################
### ENERGY MINIMIZATION ###
###############################


OUT_PATH = 'test'
name = 'alaa'

if not os.path.exists(OUT_PATH):
    os.makedirs(OUT_PATH)

logger.info(f"- Launch energy minimization")
tools.minimize(
        simulation,
        f"{OUT_PATH}/{name}_em_water.cif",
        prmtop.topology,
        maxIterations=10000,
        overwrite=False)

###############################
### PROGRESSIVE EQUILIBRATION ###
###############################

cif = PDBxFile(f"{OUT_PATH}/{name}_em_water.cif")
simulation.context.setPositions(cif.positions)
simulation.context.setPeriodicBoxVectors(*cif.topology.getPeriodicBoxVectors())

simulation.context.setVelocitiesToTemperature(temperature)

dt = 4 * unit.femtosecond

save_step_log = 10000
save_step_dcd = 10000
eq_time_expl = 3.25 #ns
tot_steps = int(np.ceil(eq_time_expl * unit.nanoseconds / dt))

for i, k in enumerate([10.5, 4.2, 2.1, 2.1, 1]):
    print(i, k)
    simulation.context.setParameter('k', k)
    logger.info(f"Equilibration step {i} with restraint k={k}")    
    tools.simulate(
        simulation,
        prmtop.topology,
        tot_steps=tot_steps,
        dt=dt,
        generic_name=f"{OUT_PATH}/{name}_explicit_equi_{i}",
        save_step_log = save_step_log,
        save_step_dcd = save_step_dcd,
    )


###############################
### PRODUCTION ###
###############################

logger.info(f"- Launch production")
time = 200 # ns

tot_steps = int(np.ceil(time * unit.nanoseconds / dt))
simulation.context.setParameter('k', 0)

tools.simulate(
    simulation,
    prmtop.topology,
    tot_steps=tot_steps,
    dt=dt,
    generic_name=f"{OUT_PATH}/{name}_explicit_prod",
    save_step_log = save_step_log,
    save_step_dcd = save_step_dcd,
)