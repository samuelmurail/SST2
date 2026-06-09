import argparse
import numpy as np
import time
import os
import sys
import logging
import pandas as pd

from openmm.app import PDBFile, PDBxFile, ForceField, Simulation
from openmm import LangevinMiddleIntegrator, unit, Platform

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
        "-seq",
        action="store",
        dest="seq",
        help="Input Sequence",
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
        "-vec",
        action="store",
        dest="vec",
        help="Box size of one vector, provide only one value, default=None",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-eq_time_impl",
        action="store",
        dest="eq_time_impl",
        help="Implicit solvent Equilibration time, default=10 (ns)",
        type=float,
        default=10,
    )
    parser.add_argument(
        "-eq_time_expl",
        action="store",
        dest="eq_time_expl",
        help="Explicit Solvent Equilibration time, default=10 (ns)",
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
        help="Langevin Integrator friction coefficient default=10.0 (ps-1)",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "-ff",
        action="store",
        dest="ff",
        help="force field, default=amber14",
        default="amber14sb",
    )
    parser.add_argument(
        "-water_ff",
        action="store",
        dest="water_ff",
        help="force field, default=tip3p",
        default="tip3p",
    )
    parser.add_argument(
        "-ace", action="store_true", dest="ace", help="Add ACE cap to N-term"
    )
    parser.add_argument(
        "-nme", action="store_true", dest="nme", help="Add NME cap to C-term"
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

    if args.ace:
        n_term = "ACE"
    else:
        n_term = None

    if args.nme:
        c_term = "NME"
    else:
        c_term = None

    tools.create_linear_peptide(
        args.seq, f"{OUT_PATH}/{name}_linear.pdb", n_term=n_term, c_term=c_term
    )

    tools.prepare_pdb(
        f"{OUT_PATH}/{name}_linear.pdb",
        f"{OUT_PATH}/{name}_fixed.cif",
        pH=7.0,
        overwrite=False,
    )

    # should be usabble soon:
    # forcefield_files = ['amber14-all.xml', 'amber14/tip3pfb.xml', 'implicit/obc2.xml']
    if args.ff.startswith("amber"):
        forcefield_files = ["amber99sbnmr.xml", "amber99_obc.xml"]
        impl_forcefield = ForceField(*forcefield_files)
    elif args.ff == "charmm36":
        forcefield_files = ["charmm36.xml", "implicit/obc1.xml"]
        impl_forcefield = ForceField(*forcefield_files)
    else:
        raise ValueError(f"Force field {args.ff} not recognized")

    logger.info(f"- Run implicit simulation")

    tools.implicit_sim(
        f"{OUT_PATH}/{name}_fixed.cif",
        impl_forcefield,
        args.eq_time_impl,
        f"{OUT_PATH}/{name}_implicit_equi",
        temp=args.min_temp,
    )

    # forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    # forcefield = ForceField(*forcefield_files)
    forcefield = tools.get_forcefield(args.ff, args.water_ff)

    if args.water_ff in ["opc3", "tip3pfb", "tip3p", "spce"]:
        model = "tip3p"
    elif args.water_ff in ["opc", "tip4pew", "tip4pfb"]:
        model = 'tip4pew'
    else:
        logger.warning(f"Water model {args.water_ff} not recognized, using tip3pfb")
        model = "tip3p"


    if args.vec is not None:
        logger.info(f"Using box vector {args.vec} nm, pad value will be ignored.")
        box_size = args.vec
        pad = None
    else:
        box_size = None
        pad = args.pad

    tools.create_water_box(
        f"{OUT_PATH}/{name}_implicit_equi.cif",
        f"{OUT_PATH}/{name}_water.cif",
        pad=pad,
        vec=box_size,
        model=model,
        forcefield=forcefield,
        overwrite=False,
    )


    ###########################
    ### BASIC EQUILIBRATION ###
    ###########################

    dt = 4 * unit.femtosecond
    temperature = args.min_temp * unit.kelvin
    max_temp = args.last_temp * unit.kelvin
    friction = args.friction / unit.picoseconds
    hydrogenMass = args.hmr * unit.amu
    rigidWater = True
    ewaldErrorTolerance = 0.0005
    nsteps = int(np.ceil(args.eq_time_expl * unit.nanoseconds / dt))

    cif = PDBxFile(f"{OUT_PATH}/{name}_water.cif")
    PDBFile.writeFile(
        cif.topology, cif.positions, open(f"{OUT_PATH}/{name}_water.pdb", "w"), True
    )

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)

    system = tools.create_sim_system(
        cif,
        forcefield=forcefield,
        temp=temperature,
        h_mass=args.hmr,
        base_force_group=1,
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
    tot_steps = int(np.ceil(args.eq_time_expl * unit.nanoseconds / dt))

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
