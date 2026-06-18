import argparse
import numpy as np
import time
import os
import sys
import logging
import pandas as pd

from openmm.app import PDBFile, PDBxFile, ForceField, Simulation
from openmm import LangevinMiddleIntegrator, unit, Platform, app

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src/")))

import SST2
from SST2.st import ST, run_st
import SST2.tools as tools


# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Add sys.sdout as handler
logger.addHandler(logging.StreamHandler(sys.stdout))


def parser_input():
    # Parse arguments :
    parser = argparse.ArgumentParser(
        description="Simulate a peptide starting from a linear conformation."
    )
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
        "-eq_time",
        action="store",
        dest="eq_time",
        help="Implicit Solvent Equilibration time, default=10 (ns)",
        type=float,
        default=10,
    )
    parser.add_argument(
        "-time",
        action="store",
        dest="time",
        help="ST time, default=10.000 (ns)",
        type=float,
        default=10000,
    )
    parser.add_argument(
        "-temp_list",
        action="store",
        dest="temp_list",
        nargs="+",
        help="SST2 temperature list, default=None",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-temp_time",
        action="store",
        dest="temp_time",
        help="ST temperature time change interval, default=2.0 (ps)",
        type=float,
        default=2.0,
    )
    parser.add_argument(
        "-log_time",
        action="store",
        dest="log_time",
        help="ST log save time interval, default= temp_time=2.0 (ps)",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-min_temp",
        action="store",
        dest="min_temp",
        help="Base temperature, default=300(K)",
        type=float,
        default=300,
    )
    parser.add_argument(
        "-last_temp",
        action="store",
        dest="last_temp",
        help="Base temperature, default=500(K)",
        type=float,
        default=500,
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
        "-temp_num",
        action="store",
        dest="temp_num",
        help="Temperature rung number, default=None (computed as function of Epot)",
        type=int,
        default=None,
    )
    parser.add_argument(
        "-friction",
        action="store",
        dest="friction",
        help="Langevin Integrator friction coefficient default=1.0 (ps-1)",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "-ff",
        action="store",
        dest="ff",
        help="force field (amber99sbnmr, charmm36, GB99dms), default=amber99sbnmr",
        default="amber99sbnmr",
    )
    parser.add_argument(
        "-water_ff",
        action="store",
        dest="water_ff",
        help="force field, default=tip3p",
        default="tip3p",
    )
    parser.add_argument("-v", action="store_true", dest="verbose", help="Verbose mode")

    return parser


if __name__ == "__main__":
    my_parser = parser_input()
    args = my_parser.parse_args()
    logger.info(args)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose mode activated")
        SST2.show_log()

    OUT_PATH = args.out_dir
    name = args.name

    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)

    tools.prepare_pdb(
        args.pdb,
        f"{OUT_PATH}/{name}_fixed.cif",
        pH=7.0,
        overwrite=False
    )

    if args.ff == "amber99sbnmr":
        forcefield_files = ["amber99sbnmr.xml", "amber99_obc.xml"]
        forcefield = ForceField(*forcefield_files)
    elif args.ff == "charmm36":
        forcefield_files = ["charmm36.xml", "implicit/obc1.xml"]
        forcefield = ForceField(*forcefield_files)
    elif args.ff == "GB99dms":
        # https://github.com/greener-group/GB99dms/blob/main/GB99dms.xml
        forcefield_files = ["GB99dms.xml"]
        forcefield = ForceField(*forcefield_files)
    else:
        raise ValueError(f"Force field {args.ff} not recognized")

    ###########################
    ### BASIC EQUILIBRATION ###
    ###########################

    dt = 4 * unit.femtosecond
    temperature = args.min_temp * unit.kelvin
    friction = args.friction / unit.picosecond
    max_temp = args.last_temp * unit.kelvin
    hydrogenMass = args.hmr * unit.amu

    cif = PDBxFile(f"{OUT_PATH}/{name}_fixed.cif")
    PDBFile.writeFile(
        cif.topology, cif.positions, open(f"{OUT_PATH}/{name}_fixed.pdb", "w"), True
    )

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)

    system = tools.create_sim_system(
        cif,
        forcefield=forcefield,
        temp=temperature,
        h_mass=args.hmr,
        base_force_group=1,
        nonbondedCutoff=3.0 * unit.nanometer,
        nonbondedMethod=app.CutoffNonPeriodic,
        constraints=app.HBonds,
        pressure=None,
    )

    # Simulation Options
    platform = Platform.getPlatformByName("CUDA")
    # platform = Platform.getPlatformByName('OpenCL')
    platformProperties = {"Precision": "single"}

    simulation = Simulation(
        cif.topology, system, integrator, platform, platformProperties
    )
    simulation.context.setPositions(cif.positions)

    logger.info(f"- Minimize system")

    tools.minimize(
        simulation,
        f"{OUT_PATH}/{name}_em_water.cif",
        cif.topology,
        maxIterations=10000,
        overwrite=False,
    )

    simulation.context.setVelocitiesToTemperature(temperature)

    save_step_log = 10000
    save_step_dcd = 10000
    tot_steps = int(np.ceil(args.eq_time * unit.nanoseconds / dt))

    logger.info(f"- Launch equilibration")
    tools.simulate(
        simulation,
        cif.topology,
        tot_steps=tot_steps,
        dt=dt,
        generic_name=f"{OUT_PATH}/{name}_explicit_equi",
        save_step_log=save_step_log,
        save_step_dcd=save_step_dcd,
    )

    ##################
    # ##  ST SIM  ####
    ##################

    if args.temp_num is None and args.temp_list is None:
        ladder_num = tools.compute_ladder_num(
            f"{OUT_PATH}/{name}_explicit_equi", temperature, args.last_temp
        )
        temperatures = None
        ladder_num *= (
            2  # By experience doubling the number of rungs is better for SST2/ST
        )
    elif args.temp_list is not None:
        ladder_num = len(args.temp_list)
        temperatures = args.temp_list
    else:
        temperatures = None
        ladder_num = args.temp_num

    tot_steps = int(np.ceil(args.time * unit.nanoseconds / dt))
    save_step_dcd = 10000
    tempChangeInterval = int(args.temp_time / dt.in_units_of(unit.picosecond)._value)
    logger.info(f"Temperature change interval = {tempChangeInterval}")

    if args.log_time is not None:
        save_step_log = int(args.log_time / dt.in_units_of(unit.picosecond)._value)
    else:
        save_step_log = tempChangeInterval

    logger.info(f"Log save interval = {save_step_log}")

    temp_list = tools.compute_temperature_list(
        minTemperature=args.min_temp,
        maxTemperature=args.last_temp,
        numTemperatures=ladder_num,
    )

    logger.info(
        f"Using temperatures : {', '.join([str(round(temp.in_units_of(unit.kelvin)._value, 2)) for temp in temp_list])}"
    )
    logger.info(f"- Launch ST simulation {temp_list}")

    run_st(
        simulation,
        cif.topology,
        f"{OUT_PATH}/{name}_ST",
        tot_steps,
        dt=dt,
        temperatures=temp_list,
        save_step_dcd=save_step_dcd,
        save_step_log=save_step_log,
        tempChangeInterval=tempChangeInterval,
    )
