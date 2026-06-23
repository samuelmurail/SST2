import argparse
import copy
import numpy as np
import math
import os
import sys
import logging
import pandas as pd
from io import StringIO


from openmm.app import PDBFile, PDBxFile, ForceField
from openmm import LangevinMiddleIntegrator, unit, app, NonbondedForce

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src/")))

import SST2
from SST2.rest2 import REST2, run_rest2
from SST2.sst2 import run_sst2
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
        help="SST2 time, default= 10 000 (ns)",
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
        help="SST2 temperature time change interval, default=2.0 (ps)",
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
        "-ref_temp",
        action="store",
        dest="ref_temp",
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
        "-exclude_Pro_omega",
        action="store_true",
        dest="exclude_Pro_omega",
        help="Exclude Proline omega dihedral scale angles",
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
        "-ace",
        action="store_true",
        dest="ace",
        help="Add ACE cap to N-term"
    )
    parser.add_argument(
        "-nme",
        action="store_true",
        dest="nme",
        help="Add NME cap to C-term"
    )
    parser.add_argument(
        "-v",
        action="store_true",
        dest="verbose",
        help="Verbose mode"
    )

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
        logger.info("Adding ACE")
        n_term = "ACE"
    else:
        n_term = None

    if args.nme:
        logger.info("Adding NME")
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
    # forcefield_files = ['amber99sbnmr.xml', 'amber99_obc.xml']
    # impl_forcefield = ForceField(*forcefield_filerfs)

    if args.ff.startswith("amber"):
        forcefield_files = ["amber99sbnmr.xml", "amber99_obc.xml"]
        impl_forcefield = ForceField(*forcefield_files)
    elif args.ff == "charmm36":
        forcefield_files = ["charmm36.xml", "implicit/obc1.xml"]
        impl_forcefield = ForceField(*forcefield_files)
    elif args.ff == "charmm36m":
        forcefield_files = ["charmm/charmm36_nowaters.xml", "implicit/obc1.xml"]
        impl_forcefield = ForceField(*forcefield_files)
    else:
        raise ValueError(f"Force field {args.ff} not recognized")

    logger.info(f"- Run implicit simulation")

    tools.implicit_sim(
        f"{OUT_PATH}/{name}_fixed.cif",
        impl_forcefield,
        args.eq_time_impl,
        f"{OUT_PATH}/{name}_implicit_equi",
        temp=args.ref_temp,
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

    #########################
    ### BASIC REST SYSTEM ###
    #########################

    dt = 4 * unit.femtosecond
    temperature = args.ref_temp * unit.kelvin
    friction = args.friction / unit.picoseconds
    hydrogenMass = args.hmr * unit.amu
    rigidWater = True
    ewaldErrorTolerance = 0.0005
    nsteps = args.eq_time_expl * unit.nanoseconds / dt

    cif = PDBxFile(f"{OUT_PATH}/{name}_water.cif")
    PDBFile.writeFile(
        cif.topology, cif.positions, open(f"{OUT_PATH}/{name}_water.pdb", "w"), True
    )

    # Get indices of the three sets of atoms.
    all_indices = [int(i.index) for i in cif.topology.atoms()]
    solute_indices = [
        int(i.index) for i in cif.topology.atoms() if i.residue.chain.id in ["A"]
    ]

    logger.info(f"- Defining solute part containing {len(solute_indices)} atoms")
    logger.info(
        f"- Defining solvent part containing {len(all_indices) - len(solute_indices)} atoms"
    )

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)

    system = tools.create_sim_system(
        cif,
        forcefield=forcefield,
        temp=temperature,
        h_mass=args.hmr,
        base_force_group=1,
    )

    charges = []
    nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
    for i in range(system.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        charges.append(charge)

    tot_charge = sum([charge._value for charge in charges])
    solute_charge = sum(
        [charges[i]._value for i in range(len(charges)) if i in solute_indices]
    )
    logger.info(f"- Total charge of the system: {tot_charge:.1f}")
    logger.info(f"- Total charge of the solute: {solute_charge:.1f}")

    if abs(tot_charge) > 0.01:
        logger.error(
            f"System is not neutral, charge = {tot_charge:.2f}. Please check the input structure."
        )

    if abs(solute_charge) > 0.01:
        logger.warning(
            f"Solute is charged, charge = {solute_charge:.2f}."
            f"Checked that your version of openmm is equal or above"
            f"8.3.1, otherwise be aware that the PME method might introduce artifacts in energy calculations.\n"
        )

    nonbondedMethod = app.PME

    sys_rest2 = REST2(
        system=system,
        pdb=cif,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator,
        dt=dt,
        temperature=temperature,
        exclude_Pro_omegas=args.exclude_Pro_omega,
        nonbondedMethod=nonbondedMethod,
    )

    logger.info(f"- Minimize system")
    tools.minimize(
        sys_rest2.simulation,
        f"{OUT_PATH}/{name}_em_water.cif",
        cif.topology,
        maxIterations=10000,
        overwrite=False,
    )

    sys_rest2.simulation.context.setVelocitiesToTemperature(temperature)

    save_step_log = 10000
    save_step_dcd = 10000
    report_rest2_Interval = 500

    logger.info(f"- Launch REST2 equilibration")

    run_rest2(
        sys_rest2,
        f"{OUT_PATH}/{name}_equi_water",
        tot_steps=nsteps,
        dt=dt,
        save_step_dcd=100000,
        save_step_log=10000,
        save_step_rest2=500,
        remove_reporters=False,
    )

    ####################
    # ##  SST2 SIM  ####
    ####################

    if args.temp_num is None and args.temp_list is None:
        ladder_num = tools.compute_ladder_num(
            f"{OUT_PATH}/{name}_equi_water_rest2",
            temperature,
            args.last_temp,
            sst2_score=True,
        )
        temperatures = None
        ladder_num *= (
            2  # By experience doubling the number of rungs is better for SST2/ST
        )
        logger.info(f"Estimated number of lambda ladder = {ladder_num}")

    elif args.temp_list is not None:
        ladder_num = len(args.temp_list)
        temperatures = args.temp_list
    else:
        temperatures = None
        ladder_num = args.temp_num

    tot_steps = args.time * unit.nanoseconds / dt
    logger.info(f"SST2 Total simulation steps = {tot_steps}")
    save_step_dcd = 10000
    # save_step_log = 100

    tempChangeInterval = int(args.temp_time / dt.in_units_of(unit.picosecond)._value)
    logger.info(f"SST2 Temp. change interval = {tempChangeInterval}")

    if args.log_time is not None:
        save_step_log = int(args.log_time / dt.in_units_of(unit.picosecond)._value)
    else:
        save_step_log = tempChangeInterval

    logger.info(f"Log save interval = {save_step_log}")

    save_check_steps = int(500.0 * unit.nanoseconds / dt)
    logger.info(f"Save checkpoint every {save_check_steps} steps")

    temp_list = tools.compute_temperature_list(
        minTemperature=args.min_temp,
        maxTemperature=args.last_temp,
        numTemperatures=ladder_num,
        refTemperature=args.ref_temp,
    )

    logger.info(
        f"Using temperatures : {', '.join([str(round(temp.in_units_of(unit.kelvin)._value, 2)) for temp in temp_list])}"
    )
    logger.info(f"- Launch SST2 simulation")

    run_sst2(
        sys_rest2,
        f"{OUT_PATH}/{name}",
        tot_steps,
        dt=dt,
        temperatures=temp_list,
        ref_temp=args.ref_temp,
        save_step_dcd=save_step_dcd,
        save_step_log=save_step_log,
        save_step_rest2=save_step_log,
        tempChangeInterval=tempChangeInterval,
        reportInterval=save_step_log,
        overwrite=False,
        save_checkpoint_steps=save_check_steps,
    )
