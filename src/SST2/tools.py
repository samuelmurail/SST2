#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import logging
import time

import numpy as np
import pandas as pd

import openmm
from openmm import unit, Platform
import openmm.app as app
import pdbfixer

# Logging
logger = logging.getLogger(__name__)


def get_forcefield(forcefield_name, water_model):
    """Create the forcefield object from the forcefield name.

    Parameters
    ----------
    forcefield_name : str
        Name of the forcefield
    water_model : str
        Name of the water forcefield

    Returns
    -------
    forcefield : openmm ForceField
        Forcefield object
    """

    forcefield_name = forcefield_name.lower().strip()

    forcefield_specs = {
        # Available forcefields in openmm, with their corresponding water models
        "amber14sb": {
            "base": ["amber14/protein.ff14SB.xml"],
            "water": {
                "tip3p": "amber14/tip3p.xml",
                "tip3pfb": "amber14/tip3pfb.xml",
                "tip4pew": "amber14/tip4pew.xml",
                "tip4pfb": "amber14/tip4pfb.xml",
                "spce": "amber14/spce.xml",
                "opc": "amber14/opc.xml",
            },
        },
        "amber99sbildn": {
            "base": ["amber99sbildn.xml"],
            "water": {
                "tip3p": "tip3p.xml",
                "tip3pfb": "tip3pfb.xml",
                "tip4pew": "tip4pew.xml",
                "tip4pfb": "tip4pfb.xml",
                "spce": "spce.xml",
                "opc": "opc.xml",
            },
        },
        "amber99sbnmr": {
            "base": ["amber99sbnm.xml"],
            "water": {
                "tip3p": "tip3p.xml",
                "tip3pfb": "tip3pfb.xml",
                "tip4pew": "tip4pew.xml",
                "tip4pfb": "tip4pfb.xml",
                "spce": "spce.xml",
                "opc": "opc.xml",
            },
        },
        # Available forcefields in openmmforcefields, with their corresponding water models
        "amber19sb": {
            "base": ["amber/protein.ff19SB.xml"],
            "water": {
                "opc3": "amber/opc3_standard.xml",
                "opc": "amber/opc_standard.xml",
                "tip3pfb": "amber/tip3pfb_standard.xml",
                "tip3p": "amber/tip3p_standard.xml",
                "tip4pew": "amber/tip4pew_standard.xml",
                "tip4pfb": "amber/tip4pfb.xml",
                "spce": "amber/spce_standard.xml",
            },
        },
        "charmm36": {
            "base": ["charmm36.xml"],
            "water": {
                "tip3": "charmm36/water.xml",
                "tip3p": "charmm36/water.xml",
                "spce": "charmm36/spce.xml",
                "tip3p-pme-b": "charmm36/tip3p-pme-b.xml",
                "tip3p-pme-f": "charmm36/tip3p-pme-f.xml",
                "tip4p2005": "charmm36/tip4p2005.xml",
                "tip4pew": "charmm36/tip4pew.xml",
                "tip5p": "charmm36/tip5p.xml",
                "tip5pew": "charmm36/tip5pew.xml",
            },
        },
        "charmm36m": {
            "base": ["charmm/charmm36_nowaters.xml"],
            "water": {
                "tip3": "charmm/waters_ions_default.xml",
                "tip3p": "charmm/waters_ions_default.xml",
                "spce": "charmm/waters_ions_spc_e.xml",
                "tip3p-pme-b": "charmm/waters_ions_tip3p_pme_b.xml",
                "tip3p-pme-f": "charmm/waters_ions_tip3p_pme_f.xml",
                "tip4p2005": "charmm/waters_ions_tip4p_2005.xml",
                "tip4pew": "charmm/waters_ions_tip4p_ew.xml",
                "tip5p": "charmm/waters_ions_tip5p.xml",
                "tip5pew": "charmm/waters_ions_tip5p_ew.xml",
            },
        },

    }

    try:
        forcefield_spec = forcefield_specs[forcefield_name]
        water_file = forcefield_spec["water"][water_model]
    except KeyError as error:
        if error.args and error.args[0] == forcefield_name:
            raise ValueError(f"Forcefield {forcefield_name} not recognized") from None
        raise ValueError(
            f"Water Forcefield {water_model} not recognized with {forcefield_name}"
        ) from None

    return app.ForceField(*forcefield_spec["base"], water_file)


def get_fastest_platform_name():
    """Get the fastest platform available on the system.

    Returns
    -------
    platform : str
        The name of the fastest platform available on the system.
    """

    platforms = {}

    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        logger.info(f"Platform {p.getName()} is available")
        platforms[p.getName()] = p

    if "CUDA" in platforms:
        return "CUDA"
    elif "HIP" in platforms:
        return "HIP"
    elif "OpenCL" in platforms:
        return "OpenCL"
    else:
        return "CPU"

def create_linear_peptide(seq, out_pdb, n_term=None, c_term=None):
    """Creates a linear peptide and save it in
    `out_pdb` file.

    Parameters
    ----------
    seq : str
        Sequence of the peptide
    out_pdb : str
        Path to the output pdb file
    n_term : str
        N-terminal residue, default is None
    c_term : str
        C-terminal residue, default is None

    Returns
    -------
    None

    """
    from pdb_numpy import abinitio

    # Create linear peptide:
    pep_coor = abinitio.make_peptide(seq, n_term=n_term, c_term=c_term)
    pep_coor.write(out_pdb)


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

    return


def implicit_sim(
    cif_in,
    forcefield,
    time,
    out_generic_name,
    temp=300 * unit.kelvin,
    dt=2 * unit.femtoseconds,
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
        cif.topology, nonbondedCutoff=3 * unit.nanometer, constraints=app.HBonds
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


def create_water_box(
    in_file,
    out_file,
    forcefield,
    pad=None,
    vec=None,
    ionicStrength=0.15 * unit.molar,
    positiveIon="Na+",
    negativeIon="Cl-",
    model='tip3p',
    overwrite=False,
):
    """Add a water box around a prepared cif file.

    Parameters
    ----------
    in_file : str
        Path to the input cif/pdb file
    out_file : str
        Path to the output cif/pdb file
    forcefield : openmm ForceField
        forcefield object
    pad : float
        Padding around the peptide in nm
    vec : float
        Vector of the box (nm), default is None
    ionicStrength : unit.Quantity
        Ionic strength of the system, default is 0.15 M
    positiveIon : str
        Positive ion, default is Na+
    negativeIon : str
        Negative ion, default is Cl-
    model : str
        Water model, default is tip3p
    overwrite : bool
        Overwrite the output file, default is False
    """

    if vec is None and pad is None:
        raise ValueError("Either pad or vec must be defined")
    if vec is not None and pad is not None:
        raise ValueError("Either pad or vec must be defined")

    if unit.is_quantity(pad):
        pad = pad.in_units_of(unit.nanometer)
    else:
        if pad is not None:
            pad = pad * unit.nanometer

    if unit.is_quantity(vec):
        vec = vec.in_units_of(unit.nanometer)
    else:
        if vec is not None:
            vec = vec * unit.nanometer

    if in_file.lower().endswith(".pdb"):
        cif = app.PDBFile(in_file)
    elif in_file.lower().endswith(".cif") or in_file.lower().endswith(".mmcif"):
        cif = app.PDBxFile(in_file)
    else:
        raise ValueError("Input file must be a pdb or cif file")

    if not overwrite and os.path.isfile(out_file):
        logger.info(f"File {out_file} exists already, skip create_water_box() step")
        if out_file.lower().endswith(".pdb"):
            cif = app.PDBFile(out_file)
        else:
            cif = app.PDBxFile(out_file)
        return cif

    # To avoid issue with clash with residues out of the box:
    x_min = min([0 * unit.nanometer] + [pos[0] for pos in cif.positions])
    y_min = min([0 * unit.nanometer] + [pos[1] for pos in cif.positions])
    z_min = min([0 * unit.nanometer] + [pos[2] for pos in cif.positions])
    min_vec = (
        openmm.Vec3(
            x_min.value_in_unit(unit.nanometer),
            y_min.value_in_unit(unit.nanometer),
            z_min.value_in_unit(unit.nanometer),
        )
        * unit.nanometer
    )
    cif.positions = [
        (pos - min_vec).value_in_unit(unit.nanometer) for pos in cif.positions
    ] * unit.nanometer

    modeller = app.Modeller(cif.topology, cif.positions)

    # Create Box

    maxSize = max(
        max((pos[i] for pos in cif.positions)) - min((pos[i] for pos in cif.positions))
        for i in range(3)
    )
    vectors = [
        openmm.Vec3(1, 0, 0),
        openmm.Vec3(1 / 3, 2 * unit.sqrt(2) / 3, 0),
        openmm.Vec3(-1 / 3, unit.sqrt(2) / 3, unit.sqrt(6) / 3),
    ]

    if vec is None:
        boxVectors = [(maxSize + pad) * v for v in vectors]
    else:
        boxVectors = [vec * v for v in vectors]
    logger.info(
        f"- Adding solvent with a {boxVectors[0][0].value_in_unit(unit.nanometer):.3} nm size box"
    )

    modeller.addSolvent(
        forcefield,
        boxVectors=boxVectors,
        ionicStrength=ionicStrength,
        positiveIon=positiveIon,
        negativeIon=negativeIon,
        model=model,
    )

    # Save
    if out_file.lower().endswith(".pdb"):
        app.PDBFile.writeFile(
            modeller.topology, modeller.positions, open(out_file, "w"), True
        )
        cif = app.PDBFile(out_file)

    elif out_file.lower().endswith(".cif") or out_file.lower().endswith(".mmcif"):
        app.PDBxFile.writeFile(
            modeller.topology, modeller.positions, open(out_file, "w"), True
        )
        cif = app.PDBxFile(out_file)

    else:
        raise ValueError("Output file must be a pdb or cif file")
    
    return cif


def create_system_simulation(
    file_io,
    forcefield,
    cif_format=True,
    dt=2 * unit.femtosecond,
    temperature=300 * unit.kelvin,
    friction=1 / unit.picoseconds,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1 * unit.nanometers,
    constraints=app.HBonds,
    platform_name="CUDA",
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
    hydrogenMass=1.0 * unit.amu,
    ignoreExternalBonds=False,
):
    """Creates a system and simulation object

    Parameters
    ----------
    file_io : str or StringIO
        Path or StringIO of the cif/pdb file
    forcefield : Openmm forcefield
        forcefield object
    cif_format : bool
        Is the file in cif format, default is True
    dt : unit.Quantity
        Time step, default is 2 fs
    temperature : unit.Quantity
        Temperature, default is 300 K
    friction : unit.Quantity
        Friction coefficient, default is 1 / ps
    nonbondedMethod : nonbonded method
        Nonbonded method, default is app.PME
    nonbondedCutoff : unit.Quantity
        Nonbonded cutoff, default is 1 nm
    constraints : constraint
        Constraints, default is app.HBonds
    platform_name : str
        Platform name, default is CUDA
    rigidWater : bool
        Rigid water, default is True
    ewaldErrorTolerance : float
        Ewald error tolerance, default is 0.0005
    hydrogenMass : unit.Quantity
        Hydrogen mass, default is 1 amu
    ignoreExternalBonds : bool
        Ignore external bonds, default is False

    Returns
    -------
    system : openmm.System
        System object
    simulation : openmm.app.Simulation
        Simulation object
    """

    if cif_format:
        pdb = app.PDBxFile(file_io)
    else:
        pdb = app.PDBFile(file_io)

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=nonbondedCutoff,
        constraints=constraints,
        rigidWater=rigidWater,
        ewaldErrorTolerance=ewaldErrorTolerance,
        hydrogenMass=hydrogenMass,
        ignoreExternalBonds=ignoreExternalBonds,
    )

    simulation = setup_simulation(
        system=system,
        position=pdb.positions,
        topology=pdb.topology,
        integrator=integrator,
        temperature=temperature,
        platform_name=platform_name
    )



    return system, simulation


def create_sim_system(
        cif,
        forcefield,
        temp=300.0,
        h_mass=1.5,
        base_force_group=1,
        rigidWater = True,
        constraints = app.HBonds,
        nonbondedMethod = app.PME,
        ewaldErrorTolerance = 0.0005,
        nonbondedCutoff = 1.0 * unit.nanometers,
        pressure = 1.0 * unit.atmospheres,
        barostatInterval = 25,
):
    """
    Create a system object from a cif file and a forcefield
    
    Parameters
    ----------
    cif : openmm.app.PDBxFile
        CIF file object
    forcefield : openmm.ForceField
        Force field object
    temp : float
        Temperature in Kelvin
    h_mass : float
        Hydrogen mass in amu
    base_force_group : int
        Base force group
    rigidWater : bool
        Rigid water flag
    constraints : openmm.app.Constraint
        Constraints
    nonbondedMethod : openmm.app.NonbondedMethod
        Nonbonded method
    ewaldErrorTolerance : float
        Ewald error tolerance
    nonbondedCutoff : float
        Nonbonded cutoff in nanometers
    pressure : float
        Pressure in atmospheres
    barostatInterval : int
        Barostat interval in steps

    Returns
    -------
    system : openmm.System
        System object
    """

    

    if unit.is_quantity(h_mass):
        hydrogenMass = h_mass.in_units_of(unit.amu)
    else:
        hydrogenMass = h_mass * unit.amu


    # Integration Options

    if unit.is_quantity(temp):
        temperature = temp.value_in_unit(unit.kelvin)
    else:
        temperature = temp * unit.kelvin

    # Prepare the Simulation

    topology = cif.topology

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=nonbondedCutoff,
        constraints=constraints,
        rigidWater=rigidWater,
        ewaldErrorTolerance=ewaldErrorTolerance,
        hydrogenMass=hydrogenMass,
    )
    if pressure is not None:
        if unit.is_quantity(pressure):
            pressure = pressure.in_units_of(unit.atmospheres)
        else:
            pressure = pressure * unit.atmospheres
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))

    for force in system.getForces():
        force.setForceGroup(base_force_group)

    return system


def minimize(simulation, out_cif, topology, maxIterations=10000, overwrite=False):
    """Minimize the energy of a system

    Parameters
    ----------
    simulation : openmm.app.Simulation
        Simulation object
    out_cif : str
        Path to the output cif file
    topology : openmm.app.Topology
        Topology object
    maxIterations : int
        Maximum number of iterations, default is 10000
    overwrite : bool
        Overwrite the output file, default is False
    """

    if not overwrite and os.path.isfile(out_cif):
        logger.info(f"File {out_cif} exists already, skip minimize() step")
        cif = app.PDBxFile(out_cif)

        # In case virtual particle are present
        # It is necessary to keep their coordinates
        positions = simulation.context.getState(
            getVelocities=False,
            getPositions=True,
            getForces=False,
            getEnergy=False,
            getParameters=False,
            groups=-1,
        ).getPositions()
        positions[: topology.getNumAtoms()] = cif.positions
        simulation.context.setPositions(positions)

        return

    simulation.minimizeEnergy(maxIterations=maxIterations)

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
        topology, positions[: topology.getNumAtoms()], open(out_cif, "w")
    )


def setup_simulation(
        system,
        position,
        topology,
        integrator,
        temperature,
        platform_name="CUDA"
):
    """Creates a simulation object

    Parameters
    ----------
    system : openmm.System
        System object
    position : unit.Quantity
        Positions
    topology : openmm.app.Topology
        Topology object
    integrator : openmm.Integrator
        Integrator object
    temperature : unit.Quantity
        Temperature
    platform_name : str
        Platform name, default is CUDA

    Returns
    -------
    simulation : openmm.app.Simulation
        Simulation object
    """

    platform = openmm.Platform.getPlatformByName(platform_name)
    prop = {}
    if platform_name != "CPU":
        prop["Precision"] = "single"

    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)

    simulation = app.Simulation(topology, system, integrator, platform, prop)
    simulation.context.setPositions(position)

    simulation.context.setVelocitiesToTemperature(temperature)
    logger.info("Created simulation")

    return simulation


def create_custom_nonbonded_force_rf(
    original_nonbonded_force,
    indexes,
    cutoff=1 * unit.nanometers,
    ONE_4PI_EPS0=138.935456,
):
    """Create a CustomNonbondedForce with Reaction Field
    to compute the solute-solvent interactions.

    Parameters
    ----------
    original_nonbonded_force : NonbondedForce
        The original NonbondedForce of the system
    indexes : list
        The list of the solute and solvent indexes
    cutoff : float * unit.nanometers
        The cutoff distance, default is 1 nm
    ONE_4PI_EPS0 : float
        The constant 1/(4*pi*epsilon0) in kJ nm / mol e2,
        default is 138.935456 kJ nm / mol e2

    Returns
    -------
    CustomNonbondedForce
        The custom nonbonded force with reaction field
    """

    eps_solvent = original_nonbonded_force.getReactionFieldDielectric()
    krf = (1 / (cutoff**3)) * (eps_solvent - 1) / (2 * eps_solvent + 1)
    crf = (1 / cutoff) * (3 * eps_solvent) / (2 * eps_solvent + 1)

    energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod*(1/r + krf*r*r - crf);"
    energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
    energy_expression += "sigma = 0.5*(sigma1+sigma2);"
    energy_expression += "krf = {:f};".format(krf.value_in_unit(unit.nanometer**-3))
    energy_expression += "crf = {:f};".format(crf.value_in_unit(unit.nanometer**-1))
    energy_expression += "chargeprod = charge1*charge2;"
    energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)

    custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)
    custom_nonbonded_force.addPerParticleParameter("charge")
    custom_nonbonded_force.addPerParticleParameter("sigma")
    custom_nonbonded_force.addPerParticleParameter("epsilon")
    custom_nonbonded_force.setNonbondedMethod(
        openmm.CustomNonbondedForce.CutoffPeriodic
    )
    custom_nonbonded_force.setCutoffDistance(cutoff)

    for index in range(original_nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = original_nonbonded_force.getParticleParameters(index)
        custom_nonbonded_force.addParticle([charge, sigma, epsilon])

    for index in range(original_nonbonded_force.getNumExceptions()):
        j, k, chargeprod, sigma, epsilon = (
            original_nonbonded_force.getExceptionParameters(index)
        )
        custom_nonbonded_force.addExclusion(j, k)

    for index in indexes:
        custom_nonbonded_force.addInteractionGroup(index[0], index[1])

    return custom_nonbonded_force


def create_custom_bonded_force_rf(
    original_nonbonded_force, atom_index, ONE_4PI_EPS0=138.935456
):
    """Create a CustomBondForce to compute the solute-solvent exceptions.
    This force is used in combination with the CustomNonbondedForce with Reaction Field.
    The CustomNonbondedForce with Reaction Field does not include exceptions,
    so we need to add them with a CustomBondForce.

    Parameters
    ----------
    original_nonbonded_force : NonbondedForce
        The original NonbondedForce of the system
    atom_index : list
        The list of the solute and solvent indexes
    ONE_4PI_EPS0 : float
        The constant 1/(4*pi*epsilon0) in kJ nm / mol e2,
        default is 138.935456 kJ nm / mol e2

    Returns
    -------
    CustomBondForce
        The custom bond force with reaction field
    """

    energy_expression = (
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r;"
    )
    energy_expression += "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)
    custom_bond_force = openmm.CustomBondForce(energy_expression)
    custom_bond_force.addPerBondParameter("chargeprod")
    custom_bond_force.addPerBondParameter("sigma")
    custom_bond_force.addPerBondParameter("epsilon")

    bond_init_param = []
    for index in range(original_nonbonded_force.getNumExceptions()):
        j, k, chargeprod, sigma, epsilon = (
            original_nonbonded_force.getExceptionParameters(index)
        )
        if j in atom_index[0] and k in atom_index[1]:
            custom_bond_force.addBond(j, k, [chargeprod, sigma, epsilon])
            bond_init_param.append((j, k, chargeprod, sigma, epsilon))

    return custom_bond_force, bond_init_param


def print_forces(system, simulation):
    """Prints the forces of the system

    Parameters
    ----------
    system : openmm.System
        System object
    simulation : openmm.app.Simulation
        Simulation object

    Returns
    -------
    None
    """

    forces_dict = get_forces(system, simulation)

    for group, force in forces_dict.items():
        print(f"{group:<3} {force['name']:<25} {force['energy']}")


def get_specific_forces(system, simulation, force_name):
    """Prints the forces of the system

    Parameters
    ----------
    system : openmm.System
        System object
    simulation : openmm.app.Simulation
        Simulation object

    Returns
    -------
    None
    """

    forces_dict = get_forces(system, simulation)

    force_tot = []

    for group, force in forces_dict.items():
        if force["name"] == force_name:
            force_tot.append(force["energy"])
    return force_tot


def get_forces(system, simulation):
    """Returns the forces of the system

    Parameters
    ----------
    system : openmm.System
        System object
    simulation : openmm.app.Simulation
        Simulation object

    Returns
    -------
    forces_dict : dict
        Dictionary with the forces
    """
    forces_dict = {}
    tot_ener = 0 * unit.kilojoules_per_mole

    for i, force in enumerate(system.getForces()):
        state = simulation.context.getState(getEnergy=True, groups={i})
        name = force.getName()
        pot_e = state.getPotentialEnergy()
        tot_ener += pot_e
        # print(f'{force.getForceGroup():<3} {name:<25} {pot_e}')

        forces_dict[force.getForceGroup()] = {"name": name, "energy": pot_e}

    forces_dict[len(forces_dict) + 1] = {"name": "Total", "energy": tot_ener}

    return forces_dict


def add_pos_restr(
    system, index_list, cif_ref, k_rest, restr_force_group=None, constant_name="k"
):
    """Add position restraints to the system

    Parameters
    ----------
    system : openmm.System
        System object
    index_list : list
        List of indices to restrain
    cif_ref : openmm.app.PDBxFile
        Reference cif file
    k_rest : float
        Force constant (KJ/mol/nm^2)
    restr_force_group : int
        Force group, default is 2
    constant_name : str
        Name of the force constant, default is k

    Returns
    -------
    restraint : openmm.CustomExternalForce
        Restraint object
    """

    restraint = openmm.CustomExternalForce(
        f"{constant_name}*periodicdistance(x, y, z, x0, y0, z0)^2"
    )
    system.addForce(restraint)
    restraint.addGlobalParameter(
        constant_name, k_rest * unit.kilojoules_per_mole / unit.nanometer**2
    )
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for index in index_list:
        restraint.addParticle(index, cif_ref.positions[index])

    if restr_force_group is not None:
        restraint.setForceGroup(restr_force_group)

    return restraint


def add_distance_restr(
    system,
    index_0_list,
    index_1_list,
    dist_min_list,
    k_rest,
    restr_force_group=None,
    constant_name="k_dist",
):
    """Add position restraints to the system

    Parameters
    ----------
    system : openmm.System
        System object
    index_0 : int
        Index of the first atom
    index_1 : int
        Index of the second atom
    dist_min : float
        Minimum distance (nm)
    k_rest : float
        Force constant (KJ/mol/nm^2)
    restr_force_group : int
        Force group, default is 2
    constant_name : str
        Name of the force constant, default is k_dist

    Returns
    -------
    restraint : openmm.CustomExternalForce
        Restraint object
    """

    assert (
        len(index_0_list) == len(index_1_list) == len(dist_min_list)
    ), "index_0, index_1 and dist_min must have the same length"

    energy_function = f"{constant_name}*(max(0, r - r0))^2"
    restraint = openmm.CustomBondForce(energy_function)
    system.addForce(restraint)
    restraint.addGlobalParameter(
        constant_name, k_rest * unit.kilojoules_per_mole / unit.nanometer**2
    )
    restraint.addPerBondParameter("r0")
    for index_0, index_1, dist_min in zip(index_0_list, index_1_list, dist_min_list):
        restraint.addBond(index_0, index_1, [dist_min * unit.nanometer])

    if restr_force_group is not None:
        restraint.setForceGroup(restr_force_group)

    return restraint


def compute_ladder_num(generic_name, min_temp, max_temp, sst2_score=False):
    """Compute the number of temperatures to simulate.

    Parameters
    ----------
    generic_name : str
        Generic name for the simulation files.
    min_temp : float
        Minimum temperature to simulate.
    max_temp : float
        Maximum temperature to simulate.
    sst2_score : bool, optional
        If True, use the SST2 score. The default is False.


    Robert Denschlag, Martin Lingenheil, Paul Tavan,
    Optimal temperature ladders in replica exchange simulations,
    Chemical Physics Letters,
    Volume 473, Issues 1–3,
    2009,

    $$ N = 1 + \frac{0.594 \sqrt{-E_{pot}}}{2 \cdot 0.534} \ln\left(\frac{T_{max}}{T_{min}}\right) $$

    Returns
    -------
    int
        Number of temperatures to simulate.
    """

    if type(min_temp) not in [int, float]:
        min_temp = min_temp._value
    if type(max_temp) not in [int, float]:
        max_temp = max_temp._value

    logger.info(f"- Extract potential energy from {generic_name}.csv")
    df_sim = pd.read_csv(generic_name + ".csv")

    # Get part number
    part = 2
    while os.path.isfile(f"{generic_name}_part_{part}.csv"):
        df_sim_part = pd.read_csv(f"{generic_name}_part_{part}.csv")
        df_sim = (
            pd.concat([df_sim, df_sim_part], axis=0, join="outer")
            .reset_index()
            .drop(["index"], axis=1)
        )
        part += 1

    # Extract potential energy
    if sst2_score:
        logger.info("- Extract SST2 potential energy")

        # New REST2 reporter format: per-fraction columns + "Solute not scaled (kJ/mole)"
        frac_cols = [col for col in df_sim.columns if col.startswith("E frac ")]
        if frac_cols:
            # New format (post-refactor)
            solute_scaled = sum(df_sim[col] for col in frac_cols)
            solute_not_scaled = df_sim["Solute not scaled (kJ/mole)"]
            solute_solvent = df_sim["Solute-Solvent (kJ/mole)"]
            df_sim["Solute(kJ/mol)"] = solute_scaled + solute_not_scaled
            df_sim["Solute-Solvent(kJ/mol)"] = solute_solvent
        elif "Solute not scaled(kJ/mol)" in df_sim.columns:
            # Legacy SST1-style format
            df_sim["Solute(kJ/mol)"] = (
                df_sim["Solute scaled(kJ/mol)"] + df_sim["Solute not scaled(kJ/mol)"]
            )
        # else: "Solute(kJ/mol)" already present in df_sim

        df_sim["new_pot"] = (
            df_sim["Solute(kJ/mol)"]
            + 0.5 * (max_temp / min_temp) ** 0.5 * df_sim["Solute-Solvent(kJ/mol)"]
        )
        E_pot = df_sim["new_pot"].mean()
    else:
        E_pot = df_sim["Potential Energy (kJ/mole)"].mean()

    if E_pot > 0:
        logger.warning(
            f"Average potential energy is positive ({E_pot:.2e} KJ.mol-1). "
            "The number of replicas cannot be estimated. "
            "Set it manually."
        )
        return 1

    logger.info(f"Average Epot = {E_pot:.2e} KJ.mol-1")
    # TO CHECK
    E_pot *= 8.314462618e-3
    logger.info(f"Average Epot = {E_pot:.2e} Kb")
    N_Nadler = 1 + 0.594 * np.sqrt(-E_pot) * np.log(max_temp / min_temp)
    logger.info(f"Nadler and Hansmann N = {N_Nadler:.2f}")
    N_Denshlag = 1 + (np.sqrt(-E_pot) / (2 * 0.534) - 0.5) * np.log(max_temp / min_temp)
    logger.info(f"Denshlag et al. N = {N_Denshlag:.2f}")
    N_Denshlag_2 = 1 + (0.594 * np.sqrt(-E_pot) - 1 / 2) * np.log(max_temp / min_temp)
    logger.info(f"Denshlag et al. 2 N = {N_Denshlag_2:.2f}")

    # print(f'\nHere N = {len(temp_list):.2f}')
    logger.info(f"Here N = {np.ceil(2*N_Denshlag):.2f}")

    return int(np.ceil(N_Denshlag))


def compute_temperature_list(
    minTemperature, maxTemperature, numTemperatures, refTemperature=None
):
    """Compute the list of temperatures to simulate.

    Parameters
    ----------
    minTemperature : float
        Minimum temperature to simulate.
    maxTemperature : float
        Maximum temperature to simulate.
    numTemperatures : int
        Number of temperatures to simulate.
    refTemperature : float, optional
        Reference temperature. The default is None.

    """

    if unit.is_quantity(minTemperature):
        minTemperature = minTemperature.in_units_of(unit.kelvin)
    else:
        minTemperature *= unit.kelvin

    if unit.is_quantity(maxTemperature):
        maxTemperature = maxTemperature.in_units_of(unit.kelvin)
    else:
        maxTemperature *= unit.kelvin

    if refTemperature is not None:
        if unit.is_quantity(refTemperature):
            refTemperature = refTemperature.in_units_of(unit.kelvin)
        else:
            refTemperature *= unit.kelvin

    # Case with refTemp is minTemp
    temperatures = [
        minTemperature
        * ((maxTemperature / minTemperature) ** (i / float(numTemperatures - 1)))
        for i in range(numTemperatures)
    ]
    if refTemperature is None or refTemperature == minTemperature:
        refTemperature = minTemperature
    else:
        # Get closest temp to ref temp
        diff_temp = [abs(temp - refTemperature) for temp in temperatures]
        ref_index = diff_temp.index(min(diff_temp))

        if ref_index > 0:
            temperatures = [
                minTemperature * ((refTemperature / minTemperature) ** (i / ref_index))
                for i in range(ref_index)
            ]
            temperatures += [
                refTemperature
                * ((maxTemperature / refTemperature))
                ** (i / (numTemperatures - ref_index - 1))
                for i in range(numTemperatures - ref_index)
            ]
        else:
            temperatures = [minTemperature] + [
                refTemperature
                * ((maxTemperature / refTemperature)) ** (i / (numTemperatures - 2))
                for i in range(numTemperatures - 1)
            ]

    return temperatures


def simulate(
    simulation,
    topology,
    tot_steps,
    dt,
    generic_name,
    additional_reporters=[],
    save_step_dcd=10000,
    save_step_log=10000,
    remove_reporters=True,
    save_checkpoint_steps=None,
    overwrite=False,
):
    """Run the simulation.

    Parameters
    ----------
    simulation : openmm.app.Simulation
        Simulation object.
    topology : openmm.app.Topology
        Topology object.
    tot_steps : int
        Total number of steps to run.
    dt : float
        Time step.
    generic_name : str
        Generic name for output files.
    additional_reporters : list, optional
        List of additional reporters. The default is [].
    save_step_dcd : int, optional
        Step to save dcd file. The default is 10000.
    save_step_log : int, optional
        Step to save log file. The default is 10000.
    save_checkpoint_steps : int, optional
        Step to save consecutive checkpoint file. The default is None.
    overwrite : bool, optional
        Overwrite previous simulation. The default is False.
    """
    tot_steps = int(tot_steps)
    final_step = tot_steps

    if not overwrite and os.path.isfile(generic_name + "_final.xml"):
        logger.info(
            f"File {generic_name}_final.xml exists already, skip simulate() step"
        )
        simulation.loadState(generic_name + "_final.xml")
        return
    elif not overwrite and os.path.isfile(generic_name + ".xml"):
        logger.info(f"File {generic_name}.xml exists, restart simulate()")
        simulation.loadState(f"{generic_name}.xml")

        # Get part number
        part = 2
        last_out_data = generic_name + ".csv"
        while os.path.isfile(f"{generic_name}_part_{part}.csv"):
            last_out_data = f"{generic_name}_part_{part}.csv"
            part += 1

        # Get last step of checkpoint:
        df_sim = pd.read_csv(last_out_data)
        chk_step = df_sim['#"Step"'][df_sim['#"Step"'] % save_step_dcd == 0].iloc[-1]

        # Bug with dcd file and step larger than 2147483647
        if chk_step >= 2147483647:
            simulation.currentStep = 0
        else:
            simulation.currentStep = int(chk_step)

        tot_steps -= chk_step
        out_name = f"{generic_name}_part_{part}"
    else:
        simulation.currentStep = 0
        out_name = generic_name

    dcd_reporter = app.DCDReporter(f"{out_name}.dcd", save_step_dcd)

    data_reporter = app.StateDataReporter(
        f"{out_name}.csv",
        save_step_log,
        totalSteps=final_step,
        step=True,
        potentialEnergy=True,
        totalEnergy=True,
        speed=True,
        temperature=True,
    )

    stdout_reporter = app.StateDataReporter(
        sys.stdout,
        save_step_dcd,
        step=True,
        temperature=True,
        speed=True,
        remainingTime=True,
        totalSteps=final_step,
    )

    check_reporter = app.CheckpointReporter(
        f"{out_name}.xml", save_step_dcd, writeState=True
    )

    # Simulation
    if remove_reporters:
        simulation.reporters = []
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(stdout_reporter)
    simulation.reporters.append(data_reporter)
    simulation.reporters.append(check_reporter)

    for reporter in additional_reporters:
        simulation.reporters.append(reporter)

    logger.info(f"Launch simulation of {tot_steps} steps")

    run_sim_check_time(
        simulation,
        tot_steps,
        dt,
        save_checkpoint_steps=save_checkpoint_steps,
        chekpoint_name=generic_name,
    )

    # simulation.step(tot_steps)

    simulation.saveState(generic_name + "_final.xml")

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
        topology, positions[: topology.getNumAtoms()], open(f"{generic_name}.cif", "w")
    )
    app.PDBFile.writeFile(
        topology, positions[: topology.getNumAtoms()], open(f"{generic_name}.pdb", "w")
    )


def run_sim_check_time(
    simulation, nsteps, dt, save_checkpoint_steps=None, chekpoint_name=None
):
    """Run a simulation and check the time

    Parameters
    ----------
    simulation : openmm.app.Simulation
        Simulation object
    nsteps : int
        Number of steps
    dt : unit.Quantity
        Time step
    save_checkpoint_steps : int
        Number of steps between each checkpoint
    chekpoint_name : str
        Name of the checkpoint file

    """

    logger.info(f"Timing {nsteps} steps of integration...")
    initial_time = time.time()
    tot_steps = nsteps

    if save_checkpoint_steps is not None:
        iter_num = int(np.ceil(nsteps / save_checkpoint_steps))
        logger.info(f"{nsteps}, {save_checkpoint_steps}, {iter_num}")
    else:
        iter_num = 1
        save_checkpoint_steps = nsteps

    for i in range(iter_num):
        if nsteps > save_checkpoint_steps:
            simulation.step(save_checkpoint_steps)
        else:
            simulation.step(nsteps)

        simulation.saveState(chekpoint_name + f"_{i:04d}.xml")

        nsteps -= save_checkpoint_steps

    final_time = time.time()
    elapsed_time = (final_time - initial_time) * unit.seconds
    elapsed_time_val = elapsed_time.value_in_unit(unit.seconds)

    dt_val = dt.value_in_unit(unit.femtoseconds)
    tot_time_val = (tot_steps * dt).value_in_unit(unit.nanoseconds)

    perfomance = (
        (tot_steps * dt).value_in_unit(unit.nanoseconds)
    ) / elapsed_time.value_in_unit(unit.day)

    logger.info(
        f"{int(tot_steps):d} steps of {dt_val:.1f} fs timestep"
        + f" ({tot_time_val:.1f} ns) took {elapsed_time_val:.1f}"
        + f" s : {perfomance:.1f} ns/day"
    )
