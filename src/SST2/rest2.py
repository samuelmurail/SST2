#!/usr/bin/env python3
# coding: utf-8


from io import StringIO
import numpy as np
import pandas as pd
import sys
import os
import logging
import pdb_numpy.format

import openmm
from openmm import unit
import openmm.app as app

from .tools import (
    setup_simulation,
    create_system_simulation,
    get_forces,
    run_sim_check_time,
)

# Logging
logger = logging.getLogger(__name__)


class Rest2Reporter(object):
    """Reporter for REST2 simulation

    Attributes
    ----------
    file : string
        The file to write to
    reportInterval : int
        The interval (in time steps) at which to write frames
    rest2 : REST2
        The REST2 object to generate the report

    Methods
    -------
    describeNextReport(simulation)
        Generate a report.


    """

    def __init__(self, file, reportInterval, rest2):
        self._out = open(file, "w", buffering=1)
        self._out.write(
            "ps,Solute scaled(kJ/mol),Solute not scaled(kJ/mol),Solvent(kJ/mol),Solute-Solvent(kJ/mol)\n"
        )
        self._reportInterval = reportInterval
        self._rest2 = rest2

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        """Generate a report.
        Compute the energies of the solute and solvent and write them to the
        file (`self._out`).

        Parameters
        ----------
        state : State
            The current state of the simulation

        Returns
        -------
        None

        """

        energies = self._rest2.compute_all_energies()

        # E_solute_scaled, E_solute_not_scaled, E_solvent, solvent_solute_nb

        time = state.getTime().value_in_unit(unit.picosecond)
        self._out.write(
            f"{time},{energies[0].value_in_unit(unit.kilojoule_per_mole)},"
            f"{energies[1].value_in_unit(unit.kilojoule_per_mole)},"
            f"{energies[2].value_in_unit(unit.kilojoule_per_mole)},"
            f"{energies[3].value_in_unit(unit.kilojoule_per_mole)}\n"
        )


class REST2:
    """REST2 class

    Attributes
    ----------
    system : System
        The system to simulate
    simulation : Simulation
        The simulation object
    positions : coordinates
        The coordinates of the system
    topology : Topology
        The topology of the system
    solute_index : list
        The list of the solute index
    solvent_index : list
        The list of the solvent index
    system_forces : dict
        The dict of the system forces
    scale : float
        The scaling factor or lambda, default is 1.0

    init_nb_param : list
        The list of the initial nonbonded parameters (charge, sigma, epsilon)
    init_nb_exept_index : list
        The list of the exception indexes
    init_nb_exept_value : list
        The list of the initial nonbonded exception parameters
        (atom1, atom2, chargeProd, sigma, epsilon)

    solute_torsion_force : CustomTorsionForce
        The torsion force of the solute
    init_torsions_index : list
        The list of the torsion indexes
    init_torsions_value : list
        The list of the initial torsion parameters

    system_solute : Solute System
        The solute system
    simulation_solute : Solute Simulation
        The solute simulation
    system_forces_solute : Solute Forces
        The solute forces

    system_solvent : Solvent System]
        The solvent system
    simulation_solvent : Solvent Simulation
        The solvent simulation
    system_forces_solvent : Solvent Forces
        The solvent forces

    init_nb_exept_solute_value : list
        The list of the initial nonbonded exception parameters of the solute
        (iatom, jatom, chargeprod, sigma, epsilon)

    Methods
    -------
    compute_all_energies()
        Compute the energies of the solute and solvent
    compute_solute_energies()
        Compute the energies of the solute
    """

    def __init__(
        self,
        system,
        pdb,
        forcefield,
        solute_index,
        integrator,
        platform_name="CUDA",
        temperature=300 * unit.kelvin,
        pressure=1.0 * unit.atmospheres,
        barostatInterval=25,
        dt=2 * unit.femtosecond,
        friction=1 / unit.picoseconds,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        hydrogenMass=1.0 * unit.amu,
    ):
        """Initialize the REST2 class

        Parameters
        ----------
        system : System
            The system to simulate
        pdb : PDBFile
            The pdb file of the system
        forcefield : ForceField
            The forcefield of the system
        solute_index : list
            The list of the solute index
        integrator : Integrator
            The integrator of the system
        platform_name : str
            The name of the platform, default is "CUDA"
        temperature : float
            The temperature of the system, default is 300 K
        pressure : float
            The pressure of the system, default is 1 atm
        barostatInterval : int
            The interval of the barostat, default is 25
        dt : float
            The timestep of the system, default is 2 fs
        friction : float
            The friction of the system, default is 1 ps-1
        nonbondedMethod : str
            The nonbonded method of the system, default is PME
        nonbondedCutoff : float
            The nonbonded cutoff of the system, default is 1 nm
        constraints : str
            The constraints of the system, default is HBonds
        rigidWater : bool
            The rigid water of the system, default is True
        ewaldErrorTolerance : float
            The Ewald error tolerance of the system, default is 0.0005
        hydrogenMass : float
            The hydrogen mass of the system, default is 1 amu


        """

        self.system = system
        self.positions = pdb.positions
        self.topology = pdb.topology
        self.solute_index = solute_index
        self.solvent_index = list(
            set(range(self.system.getNumParticles())).difference(set(self.solute_index))
        )

        assert (
            len(self.solute_index) + len(self.solvent_index)
            == self.system.getNumParticles()
        )
        assert len(self.solute_index) != 0
        assert len(self.solvent_index) != 0

        self.system_forces = {
            type(force).__name__: force for force in self.system.getForces()
        }
        self.scale = 1.0

        # Extract solute nonbonded index and values
        self.find_solute_nb_index()
        # Separate solute torsion from the solvent
        self.separate_torsion_pot()
        # Extract solute torsions index and values
        self.find_torsions()
        # Create separate solute and solvent simulation
        self.create_solute_solvent_simulation(
            forcefield=forcefield,
            platform_name=platform_name,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
            friction=friction,
            dt=dt,
        )
        # Extract solute nonbonded index and values from the solute_only system
        self.find_nb_solute_system()
        self.setup_simulation(
            integrator,
            temperature=temperature,
            pressure=pressure,
            barostatInterval=barostatInterval,
            platform_name=platform_name,
        )

    def find_solute_nb_index(self):
        """Extract initial solute nonbonded indexes and values (charge, sigma, epsilon).
        Extract also excclusion indexes and values (chargeprod, sigma, epsilon)

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        nonbonded_force = self.system_forces["NonbondedForce"]

        # Copy particles
        self.init_nb_param = []
        for particle_index in range(nonbonded_force.getNumParticles()):
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(
                particle_index
            )
            self.init_nb_param.append([charge, sigma, epsilon])

        # Copy solute-solute exclusions
        self.init_nb_exept_index = []
        self.init_nb_exept_value = []

        for exception_index in range(nonbonded_force.getNumExceptions()):
            [
                iatom,
                jatom,
                chargeprod,
                sigma,
                epsilon,
            ] = nonbonded_force.getExceptionParameters(exception_index)

            if iatom in self.solute_index and jatom in self.solute_index:
                self.init_nb_exept_index.append(exception_index)
                self.init_nb_exept_value.append(
                    [iatom, jatom, chargeprod, sigma, epsilon]
                )

    def separate_torsion_pot(self):
        """Use in the REST2 case as it avoid to modify
        twice the torsion terms in the rest2 system and
        in the solute system.

        Torsion potential is separate in two groups:
        - the solute (scaled one)
        - the solvent and not scaled solute torsion.

        As improper angles are not supposed to be scaled, here we extract only
        the proper torsion angles.

        To identify proper angles we use a trick from:
        https://github.com/maccallumlab/meld/blob/master/meld/runner/transform/rest2.py

        The original torsion potential is deleted.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        energy_expression = "k*(1+cos(period*theta-phase));"

        # Create the Solvent bond and not scaled solute torsion
        solvent_torsion_force = openmm.CustomTorsionForce(energy_expression)
        solvent_torsion_force.addPerTorsionParameter("period")
        solvent_torsion_force.addPerTorsionParameter("phase")
        solvent_torsion_force.addPerTorsionParameter("k")

        # Create the Solute bond
        solute_scaled_torsion_force = openmm.CustomTorsionForce(energy_expression)
        solute_scaled_torsion_force.addPerTorsionParameter("period")
        solute_scaled_torsion_force.addPerTorsionParameter("phase")
        solute_scaled_torsion_force.addPerTorsionParameter("k")

        # Create the not scaled Solute bond
        solute_not_scaled_torsion_force = openmm.CustomTorsionForce(energy_expression)
        solute_not_scaled_torsion_force.addPerTorsionParameter("period")
        solute_not_scaled_torsion_force.addPerTorsionParameter("phase")
        solute_not_scaled_torsion_force.addPerTorsionParameter("k")

        original_torsion_force = self.system_forces["PeriodicTorsionForce"]

        bond_idxs = [sorted([i.index, j.index]) for i, j in self.topology.bonds()]

        for i in range(original_torsion_force.getNumTorsions()):
            (
                p1,
                p2,
                p3,
                p4,
                periodicity,
                phase,
                k,
            ) = original_torsion_force.getTorsionParameters(i)

            not_improper = (
                sorted([p1, p2]) in bond_idxs
                and sorted([p2, p3]) in bond_idxs
                and sorted([p3, p4]) in bond_idxs
            )

            solute_in = (
                p1 in self.solute_index
                and p2 in self.solute_index
                and p3 in self.solute_index
                and p4 in self.solute_index
            )

            solvent_in = (
                p1 in self.solvent_index
                and p2 in self.solvent_index
                and p3 in self.solvent_index
                and p4 in self.solvent_index
            )

            if solute_in and not_improper:
                solute_scaled_torsion_force.addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            elif solute_in:
                solute_not_scaled_torsion_force.addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            elif solvent_in:
                solvent_torsion_force.addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            else:
                raise ValueError("Torsion not in solute or solvent")

        self.solute_torsion_force = solute_scaled_torsion_force

        logger.info("- Add new Torsion Forces")
        self.system.addForce(solute_scaled_torsion_force)
        self.system.addForce(solute_not_scaled_torsion_force)
        self.system.addForce(solvent_torsion_force)

        logger.info("- Delete original Torsion Forces")

        for count, force in enumerate(self.system.getForces()):
            if isinstance(force, openmm.PeriodicTorsionForce):
                self.system.removeForce(count)

    def find_torsions(self):
        """Extract the initial solute torsion index and values.
        As improper angles are not supposed to be scaled, here we extract only
        the proper torsion angles.

        To identify proper angles we use a trick from:
        https://github.com/maccallumlab/meld/blob/master/meld/runner/transform/rest2.py

        Parameters
        ----------
        None

        Returns
        -------
        None

        """

        self.init_torsions_index = []
        self.init_torsions_value = []

        torsion_force = self.solute_torsion_force

        bond_idxs = [sorted([i.index, j.index]) for i, j in self.topology.bonds()]

        for i in range(torsion_force.getNumTorsions()):
            (
                p1,
                p2,
                p3,
                p4,
                [periodicity, phase, k],
            ) = torsion_force.getTorsionParameters(i)

            # Probably useless, to check
            not_improper = (
                sorted([p1, p2]) in bond_idxs
                and sorted([p2, p3]) in bond_idxs
                and sorted([p3, p4]) in bond_idxs
            )

            # Probably useless, to check
            not_solvent = (
                p1 in self.solute_index
                and p2 in self.solute_index
                and p3 in self.solute_index
                and p4 in self.solute_index
            )

            if not_improper and not_solvent:
                self.init_torsions_index.append(i)
                self.init_torsions_value.append([p1, p2, p3, p4, periodicity, phase, k])

        logger.info(f"Solute torsion number : {len(self.init_torsions_index)}")

    def create_solute_solvent_simulation(
        self,
        forcefield,
        solute_out_pdb="solute.pdb",
        solvent_out_pdb="solvent.pdb",
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        platform_name="CUDA",
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        hydrogenMass=3.0 * unit.amu,
        friction=1 / unit.picoseconds,
        dt=2 * unit.femtosecond,
    ):
        """Extract solute only and solvent only coordinates.
        A sytem and a simulation is then created for both systems.

        Parameters
        ----------
        forcefield : str
            Forcefield name
        solute_out_pdb : str
            Output pdb file name for solute, default is "solute.pdb"
        solvent_out_pdb : str
            Output pdb file name for solvent, default is "solvent.pdb"
        nonbondedMethod : Nonbonded Method
            Nonbonded method, default is app.PME
        nonbondedCutoff : float * unit.nanometers
            Nonbonded cutoff
        constraints : Constraints
            Constraints
        platform_name : str

        """

        # Save pdb coordinates to read them with pdb_numpy

        # Redirect stdout in the variable new_stdout:
        old_stdout = sys.stdout
        stdout = new_stdout = StringIO()
        # In case of dummy atoms (position restraints, ...)
        # It has to be removed from pdb files
        top_num_atom = self.topology.getNumAtoms()

        app.PDBFile.writeFile(
            self.topology, self.positions[:top_num_atom], stdout, True
        )
        sys.stdout = old_stdout

        # Read
        solute_solvent_coor = pdb_numpy.format.pdb.parse(
            new_stdout.getvalue().split("\n")
        )

        # Separate coordinates in two pdb files:
        solute_coor = solute_solvent_coor.select_index(self.solute_index)
        solute_coor.write(solute_out_pdb, overwrite=True)

        solvent_coor = solute_solvent_coor.select_index(self.solvent_index)
        solvent_coor.write(solvent_out_pdb, overwrite=True)

        # Create system and simulations:
        self.system_solute, self.simulation_solute = create_system_simulation(
            solute_out_pdb,
            forcefield=forcefield,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            platform_name=platform_name,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
            friction=friction,
            dt=dt,
        )
        self.system_forces_solute = {
            type(force).__name__: force for force in self.system_solute.getForces()
        }

        self.system_solvent, self.simulation_solvent = create_system_simulation(
            solvent_out_pdb,
            forcefield=forcefield,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            platform_name=platform_name,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
        )

        self.system_forces_solvent = {
            type(force).__name__: force for force in self.system_solvent.getForces()
        }

    def find_nb_solute_system(self):
        """Extract in the solute only system:
        - exeption indexes and values (chargeprod, sigma, epsilon)

        Solute nonbonded values are not extracted as they are identical to
        the main system. Indexes are [0 :len(nonbonded values)]

        Exception values are stored as indexes [iatom, jatom] are different.
        """

        nonbonded_force = self.system_forces_solute["NonbondedForce"]

        # Copy particles
        self.init_nb_exept_solute_value = []
        for exception_index in range(nonbonded_force.getNumExceptions()):
            [
                iatom,
                jatom,
                chargeprod,
                sigma,
                epsilon,
            ] = nonbonded_force.getExceptionParameters(exception_index)
            self.init_nb_exept_solute_value.append(
                [iatom, jatom, chargeprod, sigma, epsilon]
            )

    def setup_simulation(
        self,
        integrator,
        temperature=300 * unit.kelvin,
        pressure=1.0 * unit.atmospheres,
        barostatInterval=25,
        platform_name="CUDA",
    ):
        """Add the simulation object.

        parameters
        ----------
        integrator : openmm.Integrator
            Integrator
        temperature : float * unit.kelvin
            Temperature, default is 300 * unit.kelvin
        pressure : float * unit.atmospheres
            Pressure, default is 1.0 * unit.atmospheres
        barostatInterval : int
            Barostat interval, default is 25
        platform_name : str
            Platform name, default is "CUDA"

        """

        # Add PT MonteCarlo barostat
        self.system.addForce(
            openmm.MonteCarloBarostat(pressure, temperature, barostatInterval)
        )

        self.simulation = setup_simulation(
            self.system,
            self.positions,
            self.topology,
            integrator=integrator,
            platform_name=platform_name,
        )

    def compute_solute_solvent_system_energy(self):
        """Update solute only and solvent only systems
        coordinates and box vector according to the solute-solvent
        system values.
        Extract then forces for each systems.

        Returns
        -------
        forces_solute : list of float * unit.kilojoules_per_mole / unit.nanometers
            Forces on solute
        forces_solvent : list of float * unit.kilojoules_per_mole / unit.nanometers
            Forces on solvent
        """

        sim_state = self.simulation.context.getState(getPositions=True, getEnergy=True)

        tot_positions = sim_state.getPositions(asNumpy=True)
        box_vector = sim_state.getPeriodicBoxVectors()

        self.simulation_solute.context.setPeriodicBoxVectors(*box_vector)
        self.simulation_solute.context.setPositions(tot_positions[self.solute_index])

        forces_solute = get_forces(self.system_solute, self.simulation_solute)

        self.simulation_solvent.context.setPeriodicBoxVectors(*box_vector)
        self.simulation_solvent.context.setPositions(tot_positions[self.solvent_index])

        forces_solvent = get_forces(self.system_solvent, self.simulation_solvent)

        return (forces_solute, forces_solvent)

    def update_torsions(self, scale):
        """Scale system solute torsion by a scale factor."""

        torsion_force = self.solute_torsion_force

        for i, index in enumerate(self.init_torsions_index):
            p1, p2, p3, p4, periodicity, phase, k = self.init_torsions_value[i]
            torsion_force.setTorsionParameters(
                index, p1, p2, p3, p4, [periodicity, phase, k * scale]
            )

        torsion_force.updateParametersInContext(self.simulation.context)

    def update_nonbonded(self, scale):
        """Scale system nonbonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        nonbonded_force = self.system_forces["NonbondedForce"]

        for i in self.solute_index:
            q, sigma, eps = self.init_nb_param[i]
            nonbonded_force.setParticleParameters(
                i, q * np.sqrt(scale), sigma, eps * scale
            )

        for i in range(len(self.init_nb_exept_index)):
            index = self.init_nb_exept_index[i]
            p1, p2, q, sigma, eps = self.init_nb_exept_value[i]
            nonbonded_force.setExceptionParameters(
                index, p1, p2, q * scale, sigma, eps * scale
            )
        # Need to fix simulation
        nonbonded_force.updateParametersInContext(self.simulation.context)

    def update_nonbonded_solute(self, scale):
        """Scale solute only system nonbonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        nonbonded_force = self.system_forces_solute["NonbondedForce"]

        for i in range(len(self.solute_index)):
            q, sigma, eps = self.init_nb_param[i]
            nonbonded_force.setParticleParameters(
                i, q * np.sqrt(scale), sigma, eps * scale
            )

        for i in range(len(self.init_nb_exept_index)):
            p1, p2, q, sigma, eps = self.init_nb_exept_solute_value[i]
            nonbonded_force.setExceptionParameters(
                i, p1, p2, q * scale, sigma, eps * scale
            )
        # Need to fix simulation
        nonbonded_force.updateParametersInContext(self.simulation_solute.context)

    def scale_nonbonded_torsion(self, scale):
        """Scale solute nonbonded potential and
        solute torsion potential
        """

        self.scale = scale
        self.update_nonbonded(scale)
        self.update_nonbonded_solute(scale)
        self.update_torsions(scale)

    def compute_all_energies(self):
        """Extract solute potential energy and solute-solvent interactions."""

        solute_force, solvent_force = self.compute_solute_solvent_system_energy()

        E_solute_not_scaled = 0 * unit.kilojoules_per_mole
        E_solute_scaled = 0 * unit.kilojoules_per_mole
        solute_not_scaled_term = ["HarmonicBondForce", "HarmonicAngleForce"]

        for i, force in solute_force.items():
            if force["name"] == "NonbondedForce":
                solute_nb = force["energy"]
                E_solute_scaled += force["energy"]
            elif force["name"] in solute_not_scaled_term:
                E_solute_not_scaled += force["energy"]

        solvent_term = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "NonbondedForce",
            "PeriodicTorsionForce",
        ]
        E_solvent = 0 * unit.kilojoules_per_mole
        for i, force in solvent_force.items():
            if force["name"] == "NonbondedForce":
                solvent_nb = force["energy"]
            if force["name"] in solvent_term:
                E_solvent += force["energy"]

        system_force = get_forces(self.system, self.simulation)

        solute_torsion_scaled_flag = True
        solute_torsion_not_scaled_flag = False
        system_term = [
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "PeriodicTorsionForce",
            "CustomTorsionForce",
        ]

        for i, force in system_force.items():
            if force["name"] == "NonbondedForce":
                all_nb = force["energy"]
            # Torsion flag to get first component of dihedral
            # forces (the solute one)
            if force["name"] == "CustomTorsionForce" and solute_torsion_scaled_flag:
                E_solute_scaled += force["energy"]
                solute_torsion_scaled_flag = False
                solute_torsion_not_scaled_flag = True
            if force["name"] == "CustomTorsionForce" and solute_torsion_not_scaled_flag:
                E_solute_not_scaled += force["energy"]
                solute_torsion_not_scaled_flag = False

        # Non scaled solvent-solute_non bonded:
        solvent_solute_nb = all_nb - solute_nb - solvent_nb
        # Scaled non bonded
        # solvent_solute_nb *= (1 / self.scale)**0.5
        # solute_nb *= 1 / self.scale

        return (
            (1 / self.scale) * E_solute_scaled,
            E_solute_not_scaled,
            E_solvent,
            (1 / self.scale) ** 0.5 * solvent_solute_nb,
        )

    def get_customPotEnergie(self):
        """Extract solute potential energy and solute-solvent interactions."""

        E_solute_scaled, _, _, solvent_solute_nb = self.compute_all_energies()

        return E_solute_scaled + 0.5 * (1 / self.scale) ** 0.5 * solvent_solute_nb

    ###################################
    ########## OLD FUNCTIONS ##########
    ###################################

    def separate_angle_pot(self):
        """Useless in the REST2 case as solute potential energy
        is obtain from the solute only system.
        """

        energy_expression = "(k/2)*(theta-theta0)^2;"

        # Create the Solvent bond
        solvent_angle_force = openmm.CustomAngleForce(energy_expression)
        solvent_angle_force.addPerAngleParameter("theta0")
        solvent_angle_force.addPerAngleParameter("k")

        # Create the Solute bond
        solute_angle_force = openmm.CustomAngleForce(energy_expression)
        solute_angle_force.addPerAngleParameter("theta0")
        solute_angle_force.addPerAngleParameter("k")

        original_angle_force = self.system_forces["HarmonicAngleForce"]

        for i in range(original_angle_force.getNumAngles()):
            p1, p2, p3, theta0, k = original_angle_force.getAngleParameters(i)
            if (
                p1 in self.solute_index
                and p2 in self.solute_index
                and p3 in self.solute_index
            ):
                solute_angle_force.addAngle(p1, p2, p3, [theta0, k])
            elif (
                p1 not in self.solute_index
                and p2 not in self.solute_index
                and p3 not in self.solute_index
            ):
                solvent_angle_force.addAngle(p1, p2, p3, [theta0, k])
            else:
                print("Wrong Angle !")
                exit()

        self.system.addForce(solute_angle_force)
        self.system.addForce(solvent_angle_force)

    def separate_bond_pot(self):
        """Useless in the REST2 case as solute potential energy
        is obtain from the solute only system.
        """

        energy_expression = "(k/2)*(r-length)^2;"

        # Create the Solvent bond
        solvent_bond_force = openmm.CustomBondForce(energy_expression)
        solvent_bond_force.addPerBondParameter("length")
        solvent_bond_force.addPerBondParameter("k")

        # Create the Solute bond
        solute_bond_force = openmm.CustomBondForce(energy_expression)
        solute_bond_force.addPerBondParameter("length")
        solute_bond_force.addPerBondParameter("k")

        original_bond_force = self.system_forces["HarmonicBondForce"]

        for i in range(original_bond_force.getNumBonds()):
            p1, p2, length, k = original_bond_force.getBondParameters(i)
            # print(p1, p2)

            if p1 in self.solute_index and p2 in self.solute_index:
                solute_bond_force.addBond(p1, p2, [length, k])
            elif p1 not in self.solute_index and p2 not in self.solute_index:
                solvent_bond_force.addBond(p1, p2, [length, k])
            else:
                print("Wrong bond !")
                exit()

        self.system.addForce(solute_bond_force)
        self.system.addForce(solvent_bond_force)

    ###########################################################
    ###  OLD FUNCTION To ensure nb calculation are correct  ###
    ###########################################################
    def add_custom_LJ_forces(self):
        """
        Create CustomNonbondedForce to capture solute-solute and solute-solvent interactions.
        Assumes PME is in use.

        Taken from:
        https://github.com/openmm/openmm/pull/2014
        """
        # nonbonded_force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)

        nonbonded_force = self.system_forces["NonbondedForce"]

        # Determine PME parameters from nonbonded_force
        cutoff_distance = nonbonded_force.getCutoffDistance()
        [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
        if (alpha_ewald / alpha_ewald.unit) == 0.0:
            # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance
            tol = nonbonded_force.getEwaldErrorTolerance()
            alpha_ewald = (1.0 / cutoff_distance) * np.sqrt(-np.log(2.0 * tol))
        print(alpha_ewald)

        # Create CustomNonbondedForce
        ONE_4PI_EPS0 = 138.935456
        energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
        # energy_expression += "epsilon = epsilon1*epsilon2;" # Why not epsilon = sqrt(epsilon1*epsilon2)
        energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"  # Why not epsilon = sqrt(epsilon1*epsilon2)
        energy_expression += "sigma = 0.5*(sigma1+sigma2);"
        custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)
        custom_nonbonded_force.addPerParticleParameter("sigma")
        custom_nonbonded_force.addPerParticleParameter("epsilon")
        # custom_nonbonded_force.addPerParticleParameter('soluteFlag')
        # custom_nonbonded_force.addInteractionGroup(solute, solvent)

        # Configure force
        custom_nonbonded_force.setNonbondedMethod(
            openmm.CustomNonbondedForce.CutoffPeriodic
        )
        custom_nonbonded_force.setCutoffDistance(cutoff_distance)
        custom_nonbonded_force.setUseLongRangeCorrection(True)  # True

        switch_flag = nonbonded_force.getUseSwitchingFunction()
        print("switch", switch_flag)
        if switch_flag:
            custom_nonbonded_force.setUseSwitchingFunction(True)
            switching_distance = nonbonded_force.getSwitchingDistance()
            custom_nonbonded_force.setSwitchingDistance(switching_distance)
        else:  # Truncated
            custom_nonbonded_force.setUseSwitchingFunction(False)

        # Create CustomBondForce for exceptions
        energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6)"
        custom_bond_force = openmm.CustomBondForce(energy_expression)
        custom_bond_force.addPerBondParameter("sigma")
        custom_bond_force.addPerBondParameter("epsilon")

        # Copy particles
        for particle_index in range(nonbonded_force.getNumParticles()):
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(
                particle_index
            )
            if sigma == 0:
                print(charge, sigma, epsilon)
            # solute_type = 1 if index in solute else 0
            # solute_type = 1
            # custom_nonbonded_force.addParticle([charge, sigma, epsilon, solute_type])
            custom_nonbonded_force.addParticle([sigma, epsilon])

        # Copy solute-solute exclusions
        for exception_index in range(nonbonded_force.getNumExceptions()):
            [
                iatom,
                jatom,
                chargeprod,
                sigma,
                epsilon,
            ] = nonbonded_force.getExceptionParameters(exception_index)
            custom_nonbonded_force.addExclusion(iatom, jatom)
            # if (iatom in solute) and (jatom in solute):
            custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])

        self.system.addForce(custom_nonbonded_force)
        self.system.addForce(custom_bond_force)

    def add_custom_Coulomb_forces(self):
        """
        Create CustomNonbondedForce to capture solute-solute and solute-solvent interactions.
        Assumes PME is in use.

        Taken from:
        https://github.com/openmm/openmm/pull/2014
        """
        # nonbonded_force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)

        nonbonded_force = self.system_forces["NonbondedForce"]

        # Determine PME parameters from nonbonded_force
        cutoff_distance = nonbonded_force.getCutoffDistance()
        [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
        if (alpha_ewald / alpha_ewald.unit) == 0.0:
            # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance
            tol = nonbonded_force.getEwaldErrorTolerance()
            alpha_ewald = (1.0 / cutoff_distance) * np.sqrt(-np.log(2.0 * tol))
        print(alpha_ewald)

        # Create CustomNonbondedForce
        ONE_4PI_EPS0 = 138.935456
        energy_expression = "ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*r)/r;"
        energy_expression += "ONE_4PI_EPS0 = {:f};".format(
            ONE_4PI_EPS0
        )  # already in OpenMM units
        energy_expression += "chargeprod = charge1*charge2;"
        energy_expression += "alpha_ewald = {:f};".format(
            alpha_ewald.value_in_unit_system(unit.md_unit_system)
        )
        custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)
        custom_nonbonded_force.addPerParticleParameter("charge")
        # custom_nonbonded_force.addPerParticleParameter('soluteFlag')
        # custom_nonbonded_force.addInteractionGroup(solute, solvent)

        # Configure force
        custom_nonbonded_force.setNonbondedMethod(
            openmm.CustomNonbondedForce.CutoffPeriodic
        )
        custom_nonbonded_force.setCutoffDistance(cutoff_distance)
        custom_nonbonded_force.setUseLongRangeCorrection(False)
        switch_flag = nonbonded_force.getUseSwitchingFunction()
        print("switch", switch_flag)
        if switch_flag:
            custom_nonbonded_force.setUseSwitchingFunction(True)
            switching_distance = nonbonded_force.getSwitchingDistance()
            custom_nonbonded_force.setSwitchingDistance(switching_distance)
        else:  # Truncated
            custom_nonbonded_force.setUseSwitchingFunction(False)

        # Create CustomBondForce for exceptions
        energy_expression = "ONE_4PI_EPS0*chargeProd_exceptions/r"
        energy_expression += (
            "- ONE_4PI_EPS0*chargeProd_product * erf(alpha_ewald * r) / r;"
        )
        energy_expression += "ONE_4PI_EPS0 = {:f};".format(
            ONE_4PI_EPS0
        )  # already in OpenMM units
        energy_expression += "alpha_ewald = {:f};".format(
            alpha_ewald.value_in_unit_system(unit.md_unit_system)
        )
        custom_bond_force = openmm.CustomBondForce(energy_expression)
        custom_bond_force.setUsesPeriodicBoundaryConditions(True)
        custom_bond_force.addPerBondParameter("chargeProd_exceptions")
        custom_bond_force.addPerBondParameter("chargeProd_product")

        # Copy particles
        for particle_index in range(nonbonded_force.getNumParticles()):
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(
                particle_index
            )
            if sigma == 0:
                print(charge, sigma, epsilon)
            # solute_type = 1 if index in solute else 0
            # solute_type = 1
            # custom_nonbonded_force.addParticle([charge, sigma, epsilon, solute_type])
            custom_nonbonded_force.addParticle([charge])

        # Copy solute-solute exclusions
        for exception_index in range(nonbonded_force.getNumExceptions()):
            [
                iatom,
                jatom,
                chargeprod,
                sigma,
                epsilon,
            ] = nonbonded_force.getExceptionParameters(exception_index)
            custom_nonbonded_force.addExclusion(iatom, jatom)
            # if (iatom in solute) and (jatom in solute):

            # Compute chargeProd_product from original particle parameters
            p1_params = custom_nonbonded_force.getParticleParameters(iatom)
            p2_params = custom_nonbonded_force.getParticleParameters(jatom)
            # print(p1_params)
            chargeProd_product = p1_params[0] * p2_params[0]
            custom_bond_force.addBond(iatom, jatom, [chargeprod, chargeProd_product])

        self.system.addForce(custom_nonbonded_force)
        self.system.addForce(custom_bond_force)

    def add_reciprocal_force(self):
        standard_nonbonded_force = openmm.NonbondedForce()

        # Set nonbonded method and related attributes
        nonbonded_force = self.system_forces["NonbondedForce"]

        standard_nonbonded_method = nonbonded_force.getNonbondedMethod()
        standard_nonbonded_force.setNonbondedMethod(standard_nonbonded_method)
        if standard_nonbonded_method in [
            openmm.NonbondedForce.CutoffPeriodic,
            openmm.NonbondedForce.CutoffNonPeriodic,
        ]:
            epsilon_solvent = nonbonded_force.getReactionFieldDielectric()
            r_cutoff = nonbonded_force.getCutoffDistance()
            standard_nonbonded_force.setReactionFieldDielectric(epsilon_solvent)
            standard_nonbonded_force.setCutoffDistance(r_cutoff)
        elif standard_nonbonded_method in [
            openmm.NonbondedForce.PME,
            openmm.NonbondedForce.Ewald,
        ]:
            print("PME")
            [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
            delta = nonbonded_force.getEwaldErrorTolerance()
            r_cutoff = nonbonded_force.getCutoffDistance()
            standard_nonbonded_force.setPMEParameters(alpha_ewald, nx, ny, nz)
            standard_nonbonded_force.setEwaldErrorTolerance(delta)
            standard_nonbonded_force.setCutoffDistance(r_cutoff)
        # Set the use of dispersion correction
        if nonbonded_force.getUseDispersionCorrection():
            print("Dispersion Correction")
            standard_nonbonded_force.setUseDispersionCorrection(True)
        else:
            standard_nonbonded_force.setUseDispersionCorrection(False)

        # Set the use of switching function
        if nonbonded_force.getUseSwitchingFunction():
            print("Switch")
            switching_distance = nonbonded_force.getSwitchingDistance()
            standard_nonbonded_force.setUseSwitchingFunction(True)
            standard_nonbonded_force.setSwitchingDistance(switching_distance)
        else:
            standard_nonbonded_force.setUseSwitchingFunction(False)

        # Disable direct space interactions
        standard_nonbonded_force.setIncludeDirectSpace(False)

        # Iterate over particles in original nonbonded force and copy to the new nonbonded force
        for particle_idx in range(nonbonded_force.getNumParticles()):
            # Get particle terms
            q, sigma, epsilon = nonbonded_force.getParticleParameters(particle_idx)

            # Add particle
            standard_nonbonded_force.addParticle(q, sigma, epsilon)

        for exception_idx in range(nonbonded_force.getNumExceptions()):
            # Get particle indices and exception terms
            p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(
                exception_idx
            )

            # Add exception
            standard_nonbonded_force.addException(p1, p2, chargeProd, sigma, epsilon)

        self.system.addForce(standard_nonbonded_force)

    def add_custom_solute_nb_forces(self, solvent, solute):
        """
        Create CustomNonbondedForce to capture solute-solute and solute-solvent interactions.
        Assumes PME is in use.

        Taken from:
        https://github.com/openmm/openmm/pull/2014
        """

        nonbonded_force = self.system_forces["NonbondedForce"]

        # Determine PME parameters from nonbonded_force
        cutoff_distance = nonbonded_force.getCutoffDistance()
        [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
        if (alpha_ewald / alpha_ewald.unit) == 0.0:
            # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance
            tol = nonbonded_force.getEwaldErrorTolerance()
            alpha_ewald = (1.0 / cutoff_distance) * np.sqrt(-np.log(2.0 * tol))
        # print(alpha_ewald)

        # Create CustomNonbondedForce
        ONE_4PI_EPS0 = 138.935456
        energy_expression = "((4*epsilon*((sigma/r)^12 - (sigma/r)^6)"
        energy_expression += "+ ONE_4PI_EPS0*chargeprod*erfc(alpha_ewald*r)/r));"
        energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
        energy_expression += "sigma = 0.5*(sigma1+sigma2);"
        energy_expression += "ONE_4PI_EPS0 = {:f};".format(
            ONE_4PI_EPS0
        )  # already in OpenMM units
        energy_expression += "chargeprod = charge1*charge2;"
        energy_expression += "alpha_ewald = {:f};".format(
            alpha_ewald.value_in_unit_system(unit.md_unit_system)
        )
        custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)
        custom_nonbonded_force.addPerParticleParameter("charge")
        custom_nonbonded_force.addPerParticleParameter("sigma")
        custom_nonbonded_force.addPerParticleParameter("epsilon")
        # custom_nonbonded_force.addPerParticleParameter('soluteFlag')
        custom_nonbonded_force.addInteractionGroup(solute, solvent)

        # Configure force
        custom_nonbonded_force.setNonbondedMethod(
            openmm.CustomNonbondedForce.CutoffPeriodic
        )
        custom_nonbonded_force.setCutoffDistance(cutoff_distance)
        custom_nonbonded_force.setUseLongRangeCorrection(False)  # Not correct for LJ

        switch_flag = nonbonded_force.getUseSwitchingFunction()
        print("switch", switch_flag)
        if switch_flag:
            custom_nonbonded_force.setUseSwitchingFunction(True)
            switching_distance = nonbonded_force.getSwitchingDistance()
            custom_nonbonded_force.setSwitchingDistance(switching_distance)
        else:  # Truncated
            custom_nonbonded_force.setUseSwitchingFunction(False)

        # Copy particles
        self.init_nb_param = []
        for particle_index in range(nonbonded_force.getNumParticles()):
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(
                particle_index
            )
            # solute_type = 1 if index in solute else 0
            # solute_type = 1
            # custom_nonbonded_force.addParticle([charge, sigma, epsilon, solute_type])
            custom_nonbonded_force.addParticle([charge, sigma, epsilon])
            self.init_nb_param.append([charge, sigma, epsilon])

        if solvent == solute:
            # Create CustomBondForce for exceptions
            energy_expression = "((4*epsilon*((sigma/r)^12 - (sigma/r)^6)"
            energy_expression += " + ONE_4PI_EPS0*chargeprod/r));"
            energy_expression += "ONE_4PI_EPS0 = {:f};".format(
                ONE_4PI_EPS0
            )  # already in OpenMM units
            custom_bond_force = openmm.CustomBondForce(energy_expression)
            custom_bond_force.addPerBondParameter("chargeprod")
            custom_bond_force.addPerBondParameter("sigma")
            custom_bond_force.addPerBondParameter("epsilon")

        # Copy solute-solute exclusions
        self.init_nb_exept_index = []
        self.init_nb_exept_value = []
        for exception_index in range(nonbonded_force.getNumExceptions()):
            [
                iatom,
                jatom,
                chargeprod,
                sigma,
                epsilon,
            ] = nonbonded_force.getExceptionParameters(exception_index)
            if iatom in solute and jatom in solute:
                self.init_nb_exept_index.append(exception_index)
                self.init_nb_exept_value.append(
                    [iatom, jatom, chargeprod, sigma, epsilon]
                )
            custom_nonbonded_force.addExclusion(iatom, jatom)
            if solvent == solute:
                if (
                    iatom in solute and jatom in solute
                ):  # or (iatom in solvent and jatom in solvent):
                    # print(iatom, jatom, chargeprod, sigma, epsilon)
                    custom_bond_force.addBond(
                        iatom, jatom, [chargeprod, sigma, epsilon]
                    )
                    if chargeprod._value != 0.0 and epsilon._value != 0.0:
                        print(iatom, jatom, chargeprod, sigma, epsilon)
                """
                if (iatom in solute and jatom in solute):
                    print("solute", iatom, jatom, chargeprod, sigma, epsilon)
                if (iatom in solute and jatom in solvent):
                    print("solute-solvent", iatom, jatom, chargeprod, sigma, epsilon)
                #if (iatom in solvent and jatom in solvent):
                #    print("Solvent", chargeprod, sigma, epsilon)
                """

        self.system.addForce(custom_nonbonded_force)
        self.system_forces["CustomNonbondedForce"] = custom_nonbonded_force
        self.add_negative_force(
            custom_nonbonded_force, name="NegativeCustomNonbondedForce"
        )

        if solvent == solute:
            self.system.addForce(custom_bond_force)
            self.add_negative_force(custom_bond_force)

    def update_custom_nonbonded(self, scale):
        custom_nonbonded_force = self.system_forces["CustomNonbondedForce"]

        for i in self.solute_index:
            q, sigma, eps = self.init_nb_param[i]
            custom_nonbonded_force.setParticleParameters(
                i, [q * np.sqrt(scale), sigma, eps * scale]
            )

        # Need to fix simulation
        custom_nonbonded_force.updateParametersInContext(self.simulation.context)

        custom_nonbonded_force = self.system_forces["NegativeCustomNonbondedForce"]

        for i in self.solute_index:
            q, sigma, eps = self.init_nb_param[i]
            custom_nonbonded_force.setParticleParameters(
                i, [q * np.sqrt(scale), sigma, eps * scale]
            )

        # Need to fix simulation
        custom_nonbonded_force.updateParametersInContext(self.simulation.context)


def run_rest2(
    sys_rest2,
    generic_name,
    tot_steps,
    dt,
    save_step_dcd=100000,
    save_step_log=500,
    save_step_rest2=500,
    rest2_reporter=True,
    overwrite=False,
    save_checkpoint_steps=None,
):
    """
    Run REST2 simulation

    Parameters
    ----------
    sys_rest2 : Rest2 object
        System to run
    generic_name : str
        Generic name for output files
    tot_steps : int
        Total number of steps to run
    dt : float
        Time step in fs
    save_step_dcd : int, optional
        Step to save dcd file, by default 100000
    save_step_log : int, optional
        Step to save log file, by default 500
    save_step_rest2 : int, optional
        Step to save rest2 file, by default 500
    rest2_reporter : bool, optional
        If True, save rest2 file, by default True
    overwrite : bool, optional
        If True, overwrite previous files, by default False
    save_checkpoint_steps : int, optional
        Step to save checkpoint file, by default None

    """

    tot_steps = np.ceil(tot_steps)

    if not overwrite and os.path.isfile(generic_name + "_final.xml"):
        logger.info(
            f"File {generic_name}_final.xml exists already, skip run_rest2() step"
        )
        sys_rest2.simulation.loadState(generic_name + "_final.xml")
        return
    elif not overwrite and os.path.isfile(generic_name + ".xml"):
        # Get part number
        part = 2
        last_out_data = generic_name + ".csv"
        while os.path.isfile(f"{generic_name}_part_{part}.csv"):
            last_out_data = f"{generic_name}_part_{part}.csv"
            part += 1

        if part != 2:
            logger.info(
                f"File {generic_name}_part_{part - 1}.xml exists, restart run_basic_rest2()"
            )
            sys_rest2.simulation.loadState(f"{generic_name}_part_{part - 1}.xml")
        else:
            logger.info(f"File {generic_name}.xml exists, restart run_basic_rest2()")
            sys_rest2.simulation.loadState(f"{generic_name}.xml")

        # Get last step of checkpoint:
        df_sim = pd.read_csv(last_out_data)

        chk_step = df_sim['#"Step"'][df_sim['#"Step"'] % save_step_dcd == 0].iloc[-1]
        # Bug with dcd file and step larger than 2147483647
        if chk_step >= 2147483647:
            sys_rest2.simulation.currentStep = 0
        else:
            sys_rest2.simulation.currentStep = int(chk_step)
        print(f"Restart at step {sys_rest2.simulation.currentStep}")

        tot_steps -= chk_step

        out_name = f"{generic_name}_part_{part}"
    else:
        sys_rest2.simulation.currentStep = 0
        out_name = generic_name

    dcd_reporter = app.DCDReporter(f"{out_name}.dcd", save_step_dcd)
    data_reporter = app.StateDataReporter(
        f"{out_name}.csv",
        save_step_log,
        totalSteps=tot_steps,
        step=True,
        potentialEnergy=True,
        totalEnergy=True,
        speed=True,
        temperature=True,
    )
    if rest2_reporter:
        rest2_reporter = Rest2Reporter(
            f"{out_name}_rest2.csv", save_step_rest2, sys_rest2
        )

    stdout_reporter = app.StateDataReporter(
        sys.stdout,
        save_step_dcd,
        step=True,
        temperature=True,
        speed=True,
        remainingTime=True,
        totalSteps=tot_steps,
    )
    check_reporter = app.CheckpointReporter(
        f"{out_name}.xml", save_step_dcd, writeState=True
    )

    sys_rest2.simulation.reporters.append(dcd_reporter)
    sys_rest2.simulation.reporters.append(stdout_reporter)
    sys_rest2.simulation.reporters.append(data_reporter)
    sys_rest2.simulation.reporters.append(check_reporter)
    if rest2_reporter:
        sys_rest2.simulation.reporters.append(rest2_reporter)

    run_sim_check_time(
        sys_rest2.simulation,
        tot_steps,
        dt,
        save_checkpoint_steps=save_checkpoint_steps,
        chekpoint_name=generic_name,
    )

    sys_rest2.simulation.saveState(generic_name + "_final.xml")

    # Save position:
    positions = sys_rest2.simulation.context.getState(
        getVelocities=False,
        getPositions=True,
        getForces=False,
        getEnergy=False,
        getParameters=False,
        groups=-1,
    ).getPositions()

    app.PDBFile.writeFile(
        sys_rest2.topology,
        positions[: sys_rest2.topology.getNumAtoms()],
        open(f"{generic_name}.pdb", "w"),
    )


if __name__ == "__main__":
    # Check energy decompostion is correct:

    # Whole system:
    OUT_PATH = "/mnt/Data_3/SST2_clean/tmp"
    name = "2HPL"

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    dt = 2 * unit.femtosecond
    temperature = 300 * unit.kelvin
    friction = 1 / unit.picoseconds

    # SYSTEM

    pdb = app.PDBFile(f"{name}_equi_water.pdb")

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=HBonds,
    )

    simulation = setup_simulation(
        system, pdb.positions, pdb.topology, integrator, "CUDA"
    )

    print("Whole system energy")
    forces_sys = print_forces(system, simulation)

    """
    nsteps=10000
    every_step = 500
    set_coor_pep = False

    scale = 1.0

    print('Timing %d steps of integration...' % nsteps)
    initial_time = time.time()
    for i in range(nsteps//every_step):
        #print('.', end='')
        simulation.step(every_step)
        scale *= 0.95
        
        #_update_custom_nonbonded(simulation, scale, solute_solvent_custom_nonbonded_force, init_nb_param)
        
        if set_coor_pep:
            tot_positions = simulation.context.getState(
                getPositions=True,
                getEnergy=True).getPositions()
            simulation_pep.context.setPositions(tot_positions[solute[0]:solute[-1]+1])
            forces = get_forces(system_pep, simulation_pep)
    print()
    final_time = time.time()
    elapsed_time = (final_time - initial_time) * unit.seconds
    integrated_time = nsteps * dt
    print(f'{nsteps} steps of {dt/unit.femtoseconds:.1f} fs timestep' +
          f'({integrated_time/unit.picoseconds:.3f} ps) took {elapsed_time/unit.seconds:.3f}'+
          f' s : {(integrated_time/elapsed_time)/(unit.nanoseconds/unit.day):.3f} ns/day')
    """

    # PEPTIDE system:

    pdb_pep = app.PDBFile(f"{name}_only_pep.pdb")

    integrator_pep = LangevinMiddleIntegrator(temperature, friction, dt)

    system_pep = forcefield.createSystem(
        pdb_pep.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=HBonds,
    )

    simulation_pep = setup_simulation(
        system_pep, pdb_pep.positions, pdb_pep.topology, integrator_pep, "CUDA"
    )

    print("Peptide forces:")
    forces_pep = print_forces(system_pep, simulation_pep)

    # NO Peptide system

    pdb_no_pep = app.PDBFile(f"{name}_no_pep.pdb")

    integrator_no_pep = LangevinMiddleIntegrator(temperature, friction, dt)

    system_no_pep = forcefield.createSystem(
        pdb_no_pep.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=HBonds,
    )

    simulation_no_pep = setup_simulation(
        system_no_pep,
        pdb_no_pep.positions,
        pdb_no_pep.topology,
        integrator_no_pep,
        "CUDA",
    )

    print("No Peptide forces:")
    forces_no_pep = print_forces(system_no_pep, simulation_no_pep)

    ####################
    # ## REST2 test ####
    ####################

    # Get indices of the three sets of atoms.
    all_indices = [int(i.index) for i in pdb.topology.atoms()]
    solvent_indices = [
        int(i.index) for i in pdb.topology.atoms() if not (i.residue.chain.id in ["B"])
    ]
    solute_indices = [
        int(i.index) for i in pdb.topology.atoms() if i.residue.chain.id in ["B"]
    ]

    print(f" {len(solute_indices)} atom in solute group")

    integrator_rest = LangevinMiddleIntegrator(temperature, friction, dt)

    test = REST2(system, pdb, forcefield, solute_indices, integrator_rest)

    print("REST2 forces 300K:")
    forces_rest2 = print_forces(test.system, test.simulation)

    (
        E_solute_scaled,
        E_solute_not_scaled,
        E_solvent,
        solvent_solute_nb,
    ) = test.compute_all_energies()
    print(f"E_solute_scaled      {E_solute_scaled}")
    print(f"E_solute_not_scaled  {E_solute_not_scaled}")
    print(f"E_solvent            {E_solvent}")
    print(f"solvent_solute_nb    {solvent_solute_nb}")

    print("\nCompare not scaled energy rest2 vs. classic:\n")
    print(
        f"HarmonicBondForce    {forces_rest2[0]['energy']/forces_sys[0]['energy']:.5e}"
    )
    print(
        f"HarmonicAngleForce   {forces_rest2[1]['energy']/forces_sys[1]['energy']:.5e}"
    )
    print("Compare scaled energy:")
    torsion_force = (
        forces_rest2[4]["energy"]
        + forces_rest2[5]["energy"]
        + forces_rest2[6]["energy"]
    )
    print(f"PeriodicTorsionForce {torsion_force/forces_sys[3]['energy']:.5e}")
    print(
        f"NonbondedForce       {forces_rest2[2]['energy']/forces_sys[2]['energy']:.5e}"
    )
    print(
        f"Total                {forces_rest2[9]['energy']/forces_sys[6]['energy']:.5e}"
    )

    print("\nCompare torsion energy rest2 vs. pep:\n")
    torsion_force = forces_rest2[4]["energy"] + forces_rest2[5]["energy"]
    print(f"PeriodicTorsionForce {torsion_force/forces_pep[3]['energy']:.5e}")

    print("\nCompare torsion energy rest2 vs. no pep:\n")
    torsion_force = forces_rest2[6]["energy"]
    print(f"PeriodicTorsionForce {torsion_force/forces_no_pep[3]['energy']:.5e}")

    print("\nCompare nonbond energy rest2 vs. no pep+pep+solvent_solute_nb:\n")
    non_bonded = (
        solvent_solute_nb + forces_pep[2]["energy"] + forces_no_pep[2]["energy"]
    )
    print(f"NonbondedForce       {non_bonded/forces_sys[2]['energy']:.5e}")

    solute_scaled_force = forces_rest2[4]["energy"] + forces_pep[2]["energy"]
    print(f"E_solute_scaled      {solute_scaled_force/E_solute_scaled:.5e}")
    solute_not_scaled_force = (
        forces_rest2[5]["energy"] + forces_pep[0]["energy"] + +forces_pep[1]["energy"]
    )
    print(f"E_solute_not_scaled  {non_bonded/forces_sys[2]['energy']:.5e}")

    print(f"E_solvent            {E_solvent/forces_no_pep[6]['energy']:.5e}")

    scale = 0.5
    test.scale_nonbonded_torsion(scale)
    print("REST2 forces 600K:")
    forces_rest2 = print_forces(test.system, test.simulation)
    (
        E_solute_scaled,
        E_solute_not_scaled,
        E_solvent,
        solvent_solute_nb,
    ) = test.compute_all_energies()
    print(f"E_solute_scaled      {E_solute_scaled}")
    print(f"E_solute_not_scaled  {E_solute_not_scaled}")
    print(f"E_solvent            {E_solvent}")
    print(f"solvent_solute_nb    {solvent_solute_nb}")

    print("\nCompare not scaled energy rest2 vs. classic:\n")
    print(
        f"HarmonicBondForce    {forces_rest2[0]['energy']/forces_sys[0]['energy']:.5e}"
    )
    print(
        f"HarmonicAngleForce   {forces_rest2[1]['energy']/forces_sys[1]['energy']:.5e}"
    )
    print("Compare scaled energy:")
    torsion_force = (
        forces_rest2[4]["energy"] * scale
        + forces_rest2[5]["energy"]
        + forces_rest2[6]["energy"]
    )
    print(f"PeriodicTorsionForce {torsion_force/forces_sys[3]['energy']:.5e}")
    print(
        f"NonbondedForce       {forces_rest2[2]['energy']/forces_sys[2]['energy']:.5e}"
    )
    print(
        f"Total                {forces_rest2[9]['energy']/forces_sys[6]['energy']:.5e}"
    )

    print("\nCompare torsion energy rest2 vs. pep:\n")
    torsion_force = forces_rest2[4]["energy"] + forces_rest2[5]["energy"]
    print(f"PeriodicTorsionForce {torsion_force/forces_pep[3]['energy']:.5e}")

    print("\nCompare torsion energy rest2 vs. no pep:\n")
    torsion_force = forces_rest2[6]["energy"]
    print(f"PeriodicTorsionForce {torsion_force/forces_no_pep[3]['energy']:.5e}")

    print("\nCompare nonbond energy rest2 vs. no pep+pep+solvent_solute_nb:\n")
    non_bonded = (
        solvent_solute_nb + forces_pep[2]["energy"] + forces_no_pep[2]["energy"]
    )
    print(f"NonbondedForce       {non_bonded/forces_sys[2]['energy']:.5e}")

    solute_scaled_force = forces_rest2[4]["energy"] + forces_pep[2]["energy"]
    print(f"E_solute_scaled      {solute_scaled_force/E_solute_scaled:.5e}")
    solute_not_scaled_force = (
        forces_rest2[5]["energy"] + forces_pep[0]["energy"] + +forces_pep[1]["energy"]
    )
    print(f"E_solute_not_scaled  {non_bonded/forces_sys[2]['energy']:.5e}")

    print(f"E_solvent            {E_solvent/forces_no_pep[6]['energy']:.5e}")

    """
    print('Timing %d steps of integration...' % nsteps)

    temp_sim = 450

    out_name = f'tmp/test_{temp_sim:d}K'
    tot_steps = 500000
    save_step_dcd = 10000
    save_step_log = 500

    scale = 300 / temp_sim

    test.simulation.reporters.append(
        DCDReporter(out_name +'_sim_temp.dcd',
                    save_step_dcd))

    test.simulation.reporters.append(
        StateDataReporter(sys.stdout, save_step_dcd,
                          step=True,
                          temperature=True,
                          speed=True,
                          remainingTime=True,
                          totalSteps=tot_steps))

    test.simulation.reporters.append(
        StateDataReporter(out_name +'_sim_temp.csv', 
                          save_step_log,
                          step=True,
                          potentialEnergy=True,
                          totalEnergy=True,
                          speed=True,
                          temperature=True))

    test.simulation.reporters.append(
        CheckpointReporter(
            out_name +'_sim_temp.chk',
            100000))

    class Rest2Reporter(object):
        def __init__(self, file, reportInterval, rest2):
            self._out = open(file, 'w')
            self._out.write("ps,Solute(kJ/mol),Solvent(kJ/mol),Solute-Solvent(kJ/mol)\n")
            self._reportInterval = reportInterval
            self._rest2 = rest2
            
        def __del__(self):
            self._out.close()

        def describeNextReport(self, simulation):
            steps = self._reportInterval - simulation.currentStep%self._reportInterval
            return (steps, False, False, False, False, None)

        def report(self, simulation, state):

            energies = self._rest2.compute_all_energies()

            time = state.getTime().value_in_unit(unit.picosecond)
            self._out.write(f'{time},{energies[0]._value},{energies[1]._value},{energies[2]._value}\n') 

    test.simulation.reporters.append(
        Rest2Reporter(
            out_name +'_energie_sim_temp.csv',
            500,
            test))

    # At 500K (βm/β0)
    # T0/Tm


    set_coor_pep = True
    nsteps = tot_steps
    every_step = 1000

    test.scale_nonbonded_torsion(scale)

    initial_time = time.time()
    for i in range(nsteps//every_step):
        #print('.', end='')
        test.simulation.step(every_step)
        
        #_update_custom_nonbonded(simulation, scale, solute_solvent_custom_nonbonded_force, init_nb_param)
        
        if set_coor_pep:
            #solute_force = test.compute_solute_energy()
            #solvent_force = test.compute_solvent_energy()

            #for i, force in solute_force.items():
            #    if force['name'] == 'NonbondedForce':
            #        solute_nb = force['energy']._value

            #for i, force in solvent_force.items():
            #    if force['name'] == 'NonbondedForce':
            #        solvent_nb = force['energy']._value

            #print(f"solute nonbonde = {solute_nb:>10.2f} solvent nonbonde = {solvent_nb:>10.2f}")

            #test.scale_nonbonded_torsion(scale)
            test.compute_all_energies()
            #test.scale_nonbonded_torsion(scale)
            #test.update_nonbonded(scale)
            #test.update_custom_nonbonded(scale)
            #test.update_torsions(scale)

            #test.scale_nonbonded_torsion(scale)
            #tot_positions = test.simulation.context.getState(
            #    getPositions=True,
            #    getEnergy=True).getPositions()
            #simulation_pep.context.setPositions(tot_positions[solute_indices[0]:solute_indices[-1]+1])
            #forces = get_forces(system_pep, simulation_pep)

    print()
    print(scale)
    final_time = time.time()
    elapsed_time = (final_time - initial_time) * unit.seconds
    integrated_time = nsteps * dt
    print(f'{nsteps} steps of {dt/unit.femtoseconds:.1f} fs timestep' +
          f'({integrated_time/unit.picoseconds:.3f} ps) took {elapsed_time/unit.seconds:.3f}'+
          f' s : {(integrated_time/elapsed_time)/(unit.nanoseconds/unit.day):.3f} ns/day')

    forces_rest2 = get_forces(test.system, test.simulation)

    for key, val in forces_rest2.items():
        print(f"{key:3} {val['name']:20} {val['energy']}")


    # Classic
    # 10000 steps of 2.0 fs timestep(20.000 ps) took 4.925 s : 350.893 ns/day
    # With separate Bonds, Angles, Torsions and 
    # specific Non bonded Solvent-Solute
    # 10000 steps of 2.0 fs timestep(20.000 ps) took 5.258 s : 328.624 ns/day
    # With separate Bonds, Angles, Torsions and 
    # specific Non bonded Solute-Solute and Solvent-Solute
    # 10000 steps of 2.0 fs timestep(20.000 ps) took 5.698 s : 303.262 ns/day

    """

    """
    vmd test_2HPL/2HPL_em_water.pdb test_2HPL/2HPL_equi_water.dcd -m 2HPL.pdb
    pbc wrap -molid 0 -first 0 -last last -compound fragment -center com -centersel "chain A and protein" -orthorhombic
    """
