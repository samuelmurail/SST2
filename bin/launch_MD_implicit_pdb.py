import argparse
import numpy as np
import time
import os
import sys
import logging
import pandas as pd

import openmm
from openmm import unit
import openmm.app as app
import pdbfixer

from openmm.app import ForceField

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src/")))

import pdb_numpy

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Add sys.sdout as handler
logger.addHandler(logging.StreamHandler(sys.stdout))


def prepare_pdb(in_pdb, out_cif, pH=7.0, overwrite=False):
    """Prepare a raw pdb file adding :
    - missing residues
    - add hydrogens at user defined pH
    - add missing atoms

    Parameters
    ----------
    in_pdb : str
        Path to the input pdb file
    out_cif : str
        Path to the output pdb file
    pH : float
        pH of the system, default is 7.0
    overwrite : bool
        Overwrite the output file, default is False
    """

    if not overwrite and os.path.isfile(out_cif):
        logger.info(f"File {out_cif} exists already, skip prepare_pdb() step")
        return

    # Fix peptide structure
    fixer = pdbfixer.PDBFixer(filename=in_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    # fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)

    app.PDBxFile.writeFile(fixer.topology, fixer.positions, open(out_cif, "w"))
    app.PDBFile.writeFile(
        fixer.topology, fixer.positions, open(out_cif[:-3] + "pdb", "w")
    )

    return


def implicit_sim(
    cif_in,
    forcefield,
    time,
    out_generic_name,
    temp=300 * unit.kelvin,
    dt=2 * unit.femtoseconds,
    kappa=None,
    min_steps=10000,
    save_coor_step=10000,
    overwrite=False,
):
    """Launch an implicit simulation with following steps:
    - Energy minimization with `maxIterations=min_steps`
    - Implicit simulation


    Parameters
    ----------
    cif_in : str
        Path to the input cif file
    forcefield : openmm ForceField
        forcefield object
    time : float
        Simulation time in ns
    out_generic_name : str
        Generic name of the output files
    temp : unit.Quantity
        Temperature, default is 300 K
    dt : unit.Quantity
        Time step, default is 2 fs
    kappa : float
        Debye-Huckel screening parameter (nm^-1), default is None
    min_steps : int
        Number of minimization steps, default is 10000
    save_coor_step : int
        Save coordinates every save_coor_step steps, default is 10000
    overwrite : bool
        Overwrite the output file, default is False

    Returns
    -------
    None
    """

    if not overwrite and os.path.isfile(f"{out_generic_name}.cif"):
        logger.info(
            f"File {out_generic_name}.cif exists already, skip implicit_sim() step"
        )
        return



    cif = app.PDBxFile(cif_in)

    system = forcefield.createSystem(
        cif.topology,
        nonbondedCutoff=3 * unit.nanometer,
        constraints=app.HBonds,
        implicitSolventKappa=kappa,
    )

    tot_steps = int(1000 * time / 0.002)

    integrator = openmm.LangevinIntegrator(temp, 1 / unit.picosecond, dt)

    simulation = app.Simulation(cif.topology, system, integrator)
    simulation.context.setPositions(cif.positions)

    # Minimize
    simulation.minimizeEnergy(maxIterations=min_steps)

    # Simulation
    simulation.reporters = []

    simulation.reporters.append(
        app.DCDReporter(f"{out_generic_name}.dcd", save_coor_step)
    )

    simulation.reporters.append(
        app.StateDataReporter(
            sys.stdout,
            save_coor_step,
            step=True,
            speed=True,
            remainingTime=True,
            totalSteps=tot_steps,
        )
    )

    simulation.reporters.append(
        app.StateDataReporter(
            f"{out_generic_name}.csv",
            save_coor_step,
            time=True,
            potentialEnergy=True,
            totalEnergy=True,
            temperature=True,
        )
    )

    logger.info(f"Launch implicit simulation of {time:.3f} ns or {tot_steps} steps")

    simulation.step(tot_steps)

    # Save position:
    positions = simulation.context.getState(
        getVelocities=False,
        getPositions=True,
        getForces=False,
        getEnergy=False,
        getParameters=False,
        groups=-1,
    ).getPositions()

    app.PDBxFile.writeFile(
        cif.topology,
        positions[: cif.topology.getNumAtoms()],
        open(f"{out_generic_name}.cif", "w"),
    )


def parser_input():
    # Parse arguments :
    parser = argparse.ArgumentParser(description="Simulate a pdb in implicit solvant.")
    parser.add_argument(
        "-pdb",
        action="store",
        dest="pdb",
        help="Input PDB file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n",
        action="store",
        dest="name",
        help="Output file name",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-dir",
        action="store",
        dest="out_dir",
        help="Output directory for intermediate files",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-pad",
        action="store",
        dest="pad",
        help="Box padding, default=1.5 nm",
        type=float,
        default=1.5,
    )
    parser.add_argument(
        "-eq_time",
        action="store",
        dest="eq_time",
        help="Implicit Solvent Equilibration time, default=10 (ns)",
        type=float,
        default=10,
    )
    parser.add_argument(
        "-hmr",
        action="store",
        dest="hmr",
        help="Hydrogen mass repartition, default=3.0 a.m.u.",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "-ref_temp",
        action="store",
        dest="ref_temp",
        help="Base temperature, default=300(K)",
        type=float,
        default=300,
    )
    return parser


if __name__ == "__main__":
    my_parser = parser_input()
    args = my_parser.parse_args()
    logger.info(args)

    OUT_PATH = args.out_dir
    name = args.name

    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)

    prepare_pdb(args.pdb, f"{OUT_PATH}/{name}_fixed.cif", pH=7.0, overwrite=False)

    # should be usable soon:
    # forcefield_files = ['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/obc2.xml']
    forcefield_files = ["amber99sbnmr.xml", "amber99_obc.xml"]
    impl_forcefield = ForceField(*forcefield_files)

    logger.info(f"- Run implicit simulation")

    implicit_sim(
        f"{OUT_PATH}/{name}_fixed.cif",
        impl_forcefield,
        args.eq_time,
        f"{OUT_PATH}/{name}_implicit_equi",
        temp=args.ref_temp,
    )
