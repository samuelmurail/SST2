import argparse
import logging
import os
import random
import sys
import time

import numpy as np
import openmm
import pandas as pd
import pdb_numpy
from openmm import LangevinMiddleIntegrator, Platform, XmlSerializer, unit
from openmm import app
from openmm.app import ForceField, PDBFile, PDBxFile, Simulation

import SST2.tools as tools
from SST2.st import ST, run_st


# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler(sys.stdout))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crd", dest='crdfile', required=True, help='Input coordinate file : rst7, xml ...',)
    parser.add_argument("--top", dest='topfile', required=True, help='Input topology file : parm7, pdb ...',)
    parser.add_argument("--pdb", dest='pdbfile', required=True, help='Input PDB file',)

    parser.add_argument("--out", dest='folder_name', default="MD_output", help='Output folder name',)
    parser.add_argument("--name", dest='file_name', default="MD", help='Output file name prefix',)
    return parser.parse_args()

def restraints(system, coor, inpcrd, fc_bb, fc_sc, fc_lpos):
    """Add positional restraints to the system.

    Args:
        system (openmm.System): OpenMM system to which restraints will be added.
        coor (pdb_numpy.Coor): Coordinate object for atom selection.
        inpcrd (app.AmberInpcrdFile): Input coordinates containing atom positions.
        fc_bb (float): Force constant for backbone restraints.
        fc_sc (float): Force constant for sidechain restraints.
        fc_lpos (float): Force constant for lipid Z-position restraints.

    Returns:
        openmm.System: System with positional restraints added.
    """
    # PROTEIN BB + SC
    posresPROT = openmm.CustomExternalForce(
        "(k_bb*is_bb + k_sc*is_sc) * periodicdistance(x,y,z,x0,y0,z0)^2"
    )

    posresPROT.addGlobalParameter("k_bb", fc_bb)
    posresPROT.addGlobalParameter("k_sc", fc_sc)
    posresPROT.addPerParticleParameter("is_bb")
    posresPROT.addPerParticleParameter("is_sc")
    posresPROT.addPerParticleParameter("x0")
    posresPROT.addPerParticleParameter("y0")
    posresPROT.addPerParticleParameter("z0")

    bb_list = coor.select_atoms("backbone and not name H*").num -1
    sc_list = coor.select_atoms(" protein and not backbone and noh").num -1

    # Backbone restraints
    if fc_bb > 0:
        for bb_atom in bb_list:
            x, y, z = inpcrd.positions[bb_atom].value_in_unit(unit.nanometers)
            posresPROT.addParticle(bb_atom, [1, 0, x, y, z])

    # Sidechain restraints
    if fc_sc > 0:
        for sc_atom in sc_list:
            x, y, z = inpcrd.positions[sc_atom].value_in_unit(unit.nanometers)
            posresPROT.addParticle(sc_atom, [0, 1, x, y, z])

    system.addForce(posresPROT)

    # LIPID Z
    posresMEMB = openmm.CustomExternalForce(
        "k_lpos * periodicdistance(0,0,z,0,0,z0)^2"
    )
    posresMEMB.addGlobalParameter("k_lpos", fc_lpos)
    posresMEMB.addPerParticleParameter("z0")

    memb_list = coor.select_atoms("resname PC and name O31 and not name H*").num

    for memb_atom in memb_list:
        z0 = inpcrd.positions[memb_atom].value_in_unit(unit.nanometers)[2]
        posresMEMB.addParticle(memb_atom, [z0])

    system.addForce(posresMEMB)

    return system


if __name__ == "__main__":

    args = parse_args()

    ###########################
    ### READING FILES ###
    ###########################

    logger.info("- Read input files")
    crd_file = args.crdfile
    top_file = args.topfile
    pdb_file = args.pdbfile
    inpcrd = app.AmberInpcrdFile(crd_file)
    prmtop = app.AmberPrmtopFile(top_file, periodicBoxVectors=inpcrd.boxVectors)
    coor = pdb_numpy.Coor(pdb_file)

    ##############################
    ### SIMULATION PARAMETERS ###
    ##############################

    logger.info("- Define simulation parameters")
    hmr = 1.5

    dt = 4 * unit.femtosecond
    temperature = 300.0 * unit.kelvin
    friction = 1.0 / unit.picoseconds
    hydrogenMass = hmr * unit.amu
    rigidWater = True
    ewaldErrorTolerance = 0.0005

    logger.info("- Set up system")
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        rigidWater=rigidWater,
        hydrogenMass=hydrogenMass,
        ewaldErrorTolerance=ewaldErrorTolerance,
    )

    logger.info("- Set up integrator and platform")
    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    platform = Platform.getPlatformByName("CUDA")
    platformProperties = {"Precision": "single"}

    ##############################
    ### Minimisation + equi1 ###
    ##############################

    save_step_log = 10000
    save_step_dcd = 10000

    system = restraints(system, coor, inpcrd, fc_bb=4000, fc_sc=2000, fc_lpos=1000)
    simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    dt = 1 * unit.femtosecond
    tot_steps = 125000  # 0.125 ns

    OUT_PATH = args.folder_name
    name = args.file_name

    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)

    logger.info("- Initial system energy")
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    logger.info("- Launch energy minimization")
    tools.minimize(
        simulation,
        f"{OUT_PATH}/{name}_em_water.cif",
        prmtop.topology,
        maxIterations=10000,
        overwrite=False,
    )

    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    simulation.context.setVelocitiesToTemperature(temperature, random.randrange(1, 10000))

    logger.info("- Equilibration step 1")
    tools.simulate(
        simulation,
        prmtop.topology,
        tot_steps=tot_steps,
        dt=dt,
        generic_name=f"{OUT_PATH}/{name}_explicit_equi_1",
        save_step_log=save_step_log,
        save_step_dcd=save_step_dcd,
    )

    ##############################
    ### equi2 --> equi6 --> prod ###
    ##############################

    barostat = openmm.MonteCarloMembraneBarostat(
        1 * unit.bar,
        0 * unit.bar * unit.nanometer,
        temperature,
        openmm.MonteCarloMembraneBarostat.XYIsotropic,
        openmm.MonteCarloMembraneBarostat.ZFree,
        100,
    )
    barostat_added = False

    Steps = [
        dict(index="equi_2", k_bb=2000, k_sc=1000, k_lpos=400, dt=1 * unit.femtosecond, tot_steps=125000, barostat=False),
        dict(index="equi_3", k_bb=1000, k_sc=500, k_lpos=400, dt=1 * unit.femtosecond, tot_steps=125000, barostat=True),
        dict(index="equi_4", k_bb=500, k_sc=200, k_lpos=200, dt=2 * unit.femtosecond, tot_steps=250000, barostat=False),
        dict(index="equi_5", k_bb=200, k_sc=50, k_lpos=40, dt=2 * unit.femtosecond, tot_steps=250000, barostat=False),
        dict(index="equi_6", k_bb=50, k_sc=0, k_lpos=0, dt=2 * unit.femtosecond, tot_steps=500000, barostat=False),
        dict(index="prod", k_bb=0, k_sc=0, k_lpos=0, dt=4 * unit.femtosecond, tot_steps=250000, barostat=False),
    ]

    for step in Steps:
        step_name = step["index"]
        logger.info(f"- launching {step_name} step")

        if step["barostat"] and not barostat_added:
            system.addForce(barostat)
            simulation.context.reinitialize(preserveState=True)
            barostat_added = True

        simulation.context.setParameter("k_bb", step["k_bb"])
        simulation.context.setParameter("k_sc", step["k_sc"])
        simulation.context.setParameter("k_lpos", step["k_lpos"])

        tools.simulate(
            simulation,
            prmtop.topology,
            tot_steps=step["tot_steps"],
            dt=step["dt"],
            generic_name=f"{OUT_PATH}/{name}_explicit_{step_name}",
            save_step_log=save_step_log,
            save_step_dcd=save_step_dcd,
        )
    logger.info("- Simulation completed")
