#!/usr/bin/env python3
# coding: utf-8


from io import StringIO
import numpy as np
import pandas as pd
import sys
import os
import logging
import pdb_numpy.format
import time

import openmm
from openmm import unit
import openmm.app as app

# Test to launch directly the script
try:
    from .tools import (
        setup_simulation,
        create_system_simulation,
        get_forces,
        simulate,
        print_forces,
        create_custom_nonbonded_force_rf,
        create_custom_bonded_force_rf,
    )
    from .topology import get_subset
except ImportError:
    from tools import (
        setup_simulation,
        create_system_simulation,
        get_forces,
        simulate,
        print_forces,
        create_custom_nonbonded_force_rf,
        create_custom_bonded_force_rf,
    )
    from topology import get_subset

# Logging
logger = logging.getLogger(__name__)


class Rest2Reporter(object):
    """Reporter for REST2 simulation.

    Attributes
    ----------
    file : string
        The file to write to
    reportInterval : int
        The interval (in time steps) at which to write frames
    rest2 : REST2
        The REST2 object to generate the report
    """

    def __init__(self, file, reportInterval, rest2):
        self._reportInterval = reportInterval
        self._rest2 = rest2

        # Build header dynamically from fractions available in this system
        self._fracs = self._rest2.fractional_terms
        
        frac_headers = ",".join(
            f"E frac {frac} (kJ/mole)" for frac in self._fracs
        )

        self._out = open(file, "w", buffering=1)
        self._out.write(
            f"Step,Lambda,"
            f"{frac_headers},"
            f"Solute not scaled (kJ/mole),"
            f"Solvent (kJ/mole),"
            f"Solute-Solvent (kJ/mole)\n"
        )

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = (
            self._reportInterval
            - simulation.currentStep % self._reportInterval
        )
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        """Generate a report.

        Compute the energies of:
        - per-fraction unscaled solute energies (one column per fraction)
        - the solute not scaled
        - the solvent
        - the solute-solvent (unscaled)

        Parameters
        ----------
        simulation : Simulation
            The current simulation
        state : State
            The current state of the simulation
        """
        energies = self._rest2.compute_all_energies()

        E_frac_dict         = energies[0]   # {frac: unscaled_energy}
        E_solute_not_scaled = energies[1]
        E_solvent           = energies[2]
        E_solute_solvent    = energies[3]

        step = state.getStepCount()

        # Per-fraction energy columns (0.0 if fraction not present in system)
        frac_values = ",".join(
            f"{E_frac_dict.get(frac, 0 * unit.kilojoules_per_mole).value_in_unit(unit.kilojoule_per_mole):.2f}"
            for frac in self._fracs
        )

        self._out.write(
            f"{step},{self._rest2.scale:.3f},"
            f"{frac_values},"
            f"{E_solute_not_scaled.value_in_unit(unit.kilojoule_per_mole):.2f},"
            f"{E_solvent.value_in_unit(unit.kilojoule_per_mole):.2f},"
            f"{E_solute_solvent.value_in_unit(unit.kilojoule_per_mole):.2f}\n"
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

    system_solvent : Solvent System
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
        exclude_Pro_omegas=False,
    ):
        """Initialize the REST2 class

        To initialize the REST2 class, the following steps are performed:
        - Extract solute nonbonded index and values
        - Separate solute's torsions from the solvent
        - Extract solute's torsions index and values
        - Create separate solute and solvent simulations
        - Extract solute's nonbonded index and values from the solute_only system
        - Setup simulation

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
        exclude_Pro_omegas : bool
            The exclusion of the proline omegas scaling, default is False
        """

        self.system = system
        self.positions = pdb.positions
        self.topology = pdb.topology
        self.solute_index = solute_index
        self.fractional_terms = []
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

        solute_charge = 0 * unit.elementary_charge
        solvent_charge = 0 * unit.elementary_charge
        nonbonded = self.system_forces["NonbondedForce"]
        for i in self.solute_index:
            charge, _, _ = nonbonded.getParticleParameters(i)
            solute_charge += charge
        for i in self.solvent_index:
            charge, _, _ = nonbonded.getParticleParameters(i)
            solvent_charge += charge

        logger.info(f"- Solute total charge : {solute_charge.value_in_unit(unit.elementary_charge):.2f}")
        logger.info(f"- Solvent total charge: {solvent_charge.value_in_unit(unit.elementary_charge):.2f}")

        if abs(solute_charge.value_in_unit(unit.elementary_charge) + solvent_charge.value_in_unit(unit.elementary_charge)) > 0.01:
            logger.error(f"System is not neutral, charge = {solute_charge.value_in_unit(unit.elementary_charge) + solvent_charge.value_in_unit(unit.elementary_charge):.2f}. Please check the input structure.")

        self.solute_charge = solute_charge
        self.solvent_charge = solvent_charge
        self.scale = 1.0

        self.CMAP_flag = False
        self.NBFIX_flag = False
        self.lj14_flag = False

        # Charmm36 and Amber19SB forcefields have CMAP potential, so we need to separate it from the solvent if it is the case
        if "CMAPTorsionForce" in self.system_forces:
            self.CMAP_flag = True
            logger.info("CMAP Founded")
            self.separate_cmap_pot()
        
        # In charmm forcefield, NBFIX is implemented as a CustomNonbondedForce, so we need to separate it from the solvent if it is the case and if the user want to scale NBFIX
        if "CustomNonbondedForce" in self.system_forces:
            self.NBFIX_flag = True
            EnergyFunction = self.system_forces["CustomNonbondedForce"].getEnergyFunction()
            logger.info(f"CustomNonbondedForce Founded, energy function: {EnergyFunction}")
            
            NBFIX_force = self.system_forces["CustomNonbondedForce"]

            logger.info(f"Add solute flag to CustomNonbondedForce and update energy function to scale NBFIX terms")
            self.add_scale_NBFIX(NBFIX_force, self.solute_index)
        
        # In charmm36 forcefield, the Lennard-Jones 1-4 interactions are implemented as a CustomBondForce, so we need to separate it from the solvent if it is the case and if the user want to scale LJ14
        if "CustomBondForce" in self.system_forces:
            self.lj14_flag = True
            logger.info("CustomBondForce Founded")
            self.add_scale_LJ14(self.system, self.solute_index)

        # Charmm36 forcefield has a specific dihedral potential improper dihedrals
        if "CustomTorsionForce" in self.system_forces:
            logger.info("CustomTorsionForce Founded, not supposed to be scaled, check if it is the case")

        # Extract solute nonbonded index and values
        self.find_solute_nb_index()
        # Separate solute torsion from the solvent
        self.separate_torsion_pot(exclude_Pro_omegas=exclude_Pro_omegas)
        # Create separate solute and solvent simulation

        if nonbondedMethod == app.PME:
            logger.info("Create systems with PME")
            self.reaction_field = False

            # Extract alpha parameters form PME NonbondedForce
            self.alpha_ewald = self.system_forces["NonbondedForce"].getPMEParameters()[0]
            if self.alpha_ewald == 0.0 * unit.nanometers**-1:
                logger.warning(
                    "Warning: alpha parameter of PME is 0.0, "
                    "alpha chosen based on the Ewald error tolerance"
                )
                tolerance = self.system_forces["NonbondedForce"].getEwaldErrorTolerance()
                cutoff = self.system_forces["NonbondedForce"].getCutoffDistance()
                logger.info(f"- Ewald error tolerance: {tolerance}")
                logger.info(f"- PME cutoff: {cutoff}")

                self.alpha_ewald = (np.sqrt(-np.log(tolerance)) / cutoff).in_units_of(unit.nanometers**-1)

            logger.info(f"- PME alpha: {self.alpha_ewald}")

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
            self.find_nb_solute_system()


        elif nonbondedMethod == app.CutoffPeriodic:

            logger.info("Create systems with Reaction Field")
            self.reaction_field = True

            self.create_rf_simulation(
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
            logger.info("- Reaction Field object forces:")
            print_forces(self.system_rf, self.simulation_rf)


        else:
            raise ValueError("nonbondedMethod not supported")


        # Extract solute nonbonded index and values from the solute_only system
        self.setup_simulation(
            integrator,
            temperature=temperature,
            pressure=pressure,
            barostatInterval=barostatInterval,
            platform_name=platform_name,
        )


        self.fractional_terms = sorted(set(self.fractional_terms))
        logger.info(f"- Fractional terms: {self.fractional_terms}")

        print("- REST2 object forces:")
        print_forces(self.system, self.simulation)
        print("- REST2 solute object forces:")
        print_forces(self.system_solute, self.simulation_solute)
        print(f"\ncompute_all_energies() output (unscaled physical energies):")
        (
            E_frac_dict,
            E_solute_not_scaled,
            E_solvent,
            E_cross_nb,
        ) = self.compute_all_energies()

        # Sum of all fractional unscaled energies (for display)
        E_solute_scaled_total = sum(
            E_frac_dict.values(), 0 * unit.kilojoules_per_mole
        )

        print(f"  E_frac_dict:")
        for frac, E_t in sorted(E_frac_dict.items()):
            print(f"    frac={frac:.4f}  E_t = {E_t}")
        print(f"  E_solute_scaled_total = {E_solute_scaled_total}  (sum of fractions)")
        print(f"  E_solute_not_scaled   = {E_solute_not_scaled}")
        print(f"  E_solvent             = {E_solvent}")
        print(f"  E_cross_nb            = {E_cross_nb}")


    def add_scale_NBFIX(self, NBFIX_force, solute_index):
        """Extract initial solute Lennard-Jones indexes and values (sigma, epsilon).

        Parameters
        ----------
        None

        Returns
        -------
        None
        """


        if 'acoef' not in NBFIX_force.getEnergyFunction() and 'bcoef' not in NBFIX_force.getEnergyFunction():
             logger.error("CustomNonbondedForce does not seem to be a NBFIX force, no acoef and bcoef found in the energy function.")
             return

        NBFIX_force.addGlobalParameter('lambda', 1.0)
            
        # Add per-particle solute flag
        NBFIX_force.addPerParticleParameter('is_solute')
        
        # Update energy expression
        old_expr = NBFIX_force.getEnergyFunction()
        new_expr = f"""
            scale * ({old_expr[:-1]});
            scale = is_solute1*is_solute2*lambda
                    + (is_solute1*(1-is_solute2) 
                    + is_solute2*(1-is_solute1))*sqrt(lambda)
                    + (1-is_solute1)*(1-is_solute2)*1.0
        """ 
        NBFIX_force.setEnergyFunction(new_expr)

        logger.info(f"- Extract solute Lennard-Jones parameters for {NBFIX_force.getNumParticles()} particles")
        self.init_LJ_param = []
        for particle_index in range(NBFIX_force.getNumParticles()):
            params = NBFIX_force.getParticleParameters(particle_index)
            NBFIX_force.setParticleParameters(particle_index, [params[0], particle_index in solute_index])
            params = NBFIX_force.getParticleParameters(particle_index)

    def add_scale_LJ14(self, system, solute_index):
        """Separate CHARMM36 LennardJones14 CustomBondForce into:
        - solvent (unscaled)
        - solute boundary (epsilon scaled by sqrt(lambda))
        - solute-solute (epsilon scaled by lambda)

        Parameters
        ----------
        system : System
            The system to modify
        solute_index : list
            The list of the solute index

        Returns
        -------
        None        
        """
        
        
        # Find the LJ14 CustomBondForce
        for count, force in enumerate(system.getForces()):
            if isinstance(force, openmm.CustomBondForce):
                if 'epsilon' in force.getEnergyFunction():
                    lj14_force = force
                    lj14_index = count
                    break
        old_expr = lj14_force.getEnergyFunction()
        logger.info(f"Old LJ14 energy function: {old_expr}")
  
        if 'epsilon' not in lj14_force.getEnergyFunction():
            logger.error("CustomBondForce does not seem to be a Lennard-Jones 1-4 force, no epsilon found in the energy function.")
            return

        if lj14_force is None:
            logger.info("No LennardJones14 CustomBondForce found, skipping.")
            return

        energy_expression = "4*epsilon*((sigma/r)^12-(sigma/r)^6);"

        # Solvent (unscaled)
        lj14_solvent = openmm.CustomBondForce(energy_expression)
        lj14_solvent.setName("LJ14_solvent")
        lj14_solvent.addPerBondParameter("sigma")
        lj14_solvent.addPerBondParameter("epsilon")

        # Boundary: 1 solute atom, epsilon scaled by sqrt(lambda)
        lj14_boundary = openmm.CustomBondForce(
            "4*lambda_lj14_12*epsilon*((sigma/r)^12-(sigma/r)^6);"
        )
        lj14_boundary.setName("LJ14_solute_boundary")
        lj14_boundary.addGlobalParameter("lambda_lj14_12", 1.0)
        lj14_boundary.addPerBondParameter("sigma")
        lj14_boundary.addPerBondParameter("epsilon")

        # Solute-solute: both atoms in solute, epsilon scaled by lambda
        lj14_solute = openmm.CustomBondForce(
            "4*lambda_lj14_22*epsilon*((sigma/r)^12-(sigma/r)^6);"
        )
        lj14_solute.setName("LJ14_solute_solute")
        lj14_solute.addGlobalParameter("lambda_lj14_22", 1.0)
        lj14_solute.addPerBondParameter("sigma")
        lj14_solute.addPerBondParameter("epsilon")

        # Distribute bonds
        for i in range(lj14_force.getNumBonds()):
            p1, p2, (sigma, epsilon) = lj14_force.getBondParameters(i)

            solute_num_atom = (p1 in solute_index) + (p2 in solute_index)

            # logger.info(f"Bond {i}: p1={p1}, p2={p2}, epsilon={epsilon:.3f}, sigma={sigma:.3f} -> solute_num_atom={solute_num_atom}")
            # logger.info(f"  [{i}] {lj14_force.getPerBondParameterName(i)}")

            if solute_num_atom == 2:
                lj14_solute.addBond(p1, p2, [sigma, epsilon])
            elif solute_num_atom == 1:
                lj14_boundary.addBond(p1, p2, [sigma, epsilon])
            else:
                lj14_solvent.addBond(p1, p2, [sigma, epsilon])

        # # Verify parameter order in original force
        # logger.info("LJ14 per-bond parameter order:")
        # for i in range(lj14_force.getNumPerBondParameters()):
        #     logger.info(f"  [{i}] {lj14_force.getPerBondParameterName(i)}")
        # # Print first bond as sanity check
        # if lj14_force.getNumBonds() > 0:
        #     p1, p2, params = lj14_force.getBondParameters(0)
        #     logger.info(f"First bond: p1={p1}, p2={p2}, params={params}")
        # if lj14_solvent.getNumBonds() > 0:
        #     p1, p2, params = lj14_solvent.getBondParameters(0)
        #     logger.info(f"First bond: p1={p1}, p2={p2}, params={params}")
        # if lj14_boundary.getNumBonds() > 0:
        #     p1, p2, params = lj14_boundary.getBondParameters(0)
        #     logger.info(f"First bond: p1={p1}, p2={p2}, params={params}")
        # if lj14_solute.getNumBonds() > 0:
        #     p1, p2, params = lj14_solute.getBondParameters(0)
        #     logger.info(f"First bond: p1={p1}, p2={p2}, params={params}")

        # Store references
        self.lj14_solvent = lj14_solvent
        self.lj14_boundary = lj14_boundary
        self.lj14_solute = lj14_solute
        self.lj14_boundary_flag = False
        self.lj14_solute_flag = False

        logger.info("- Add new LJ14 Forces")
        if lj14_solvent.getNumBonds() > 0:
            logger.info(
                f"  Adding LJ14_solvent"
                f"({lj14_solvent.getNumBonds()} bonds)"
            )
            system.addForce(lj14_solvent)
        if lj14_boundary.getNumBonds() > 0:
            logger.info(
                f"  Adding LJ14_solute_boundary "
                f"({lj14_boundary.getNumBonds()} bonds)"
            )
            self.lj14_boundary_flag = True
            system.addForce(lj14_boundary)
        if lj14_solute.getNumBonds() > 0:
            logger.info(
                f"  Adding LJ14_solute_solute "
                f"({lj14_solute.getNumBonds()} bonds)"
            )
            self.lj14_solute_flag = True
            system.addForce(lj14_solute)

        logger.info("- Delete original LJ14 Force")
        system.removeForce(lj14_index)

    def find_solute_nb_index(self):
        """Extract initial solute nonbonded indexes and values (charge, sigma, epsilon).
        Extract also exclusion indexes and values (chargeprod, sigma, epsilon)

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

    def separate_cmap_pot(self):
        """
        CMAP potential is separate in two groups:
        - the solute (scaled one)
        - the solvent

        The original cmap potential is deleted.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # extract original cmap torsion parameters
        logger.info("- Separate CMAP torsions in solvent and solute CMAP forces")

        original_cmap_force = self.system_forces["CMAPTorsionForce"]
        logger.info("- Create CMAP torsions force for solvent")
        solvent_cmap_force = openmm.CMAPTorsionForce()
        solvent_cmap_force.setName("CMAP_solvent")
        solute_cmap_dict = {}

        for i in range(original_cmap_force.getNumTorsions()):
            cmap_indexes = original_cmap_force.getTorsionParameters(i)
            solute_in = sum(
                [cmap_indexes[j + 1] in self.solute_index for j in range(8)]
            )
            solvent_in = all(
                [cmap_indexes[j + 1] in self.solvent_index for j in range(8)]
            )

            if solvent_in:
                # logger.info(f"Add CMap torsion {i} in solvent")
                solvent_cmap_force.addTorsion(*cmap_indexes)
            elif solute_in > 0:
                # Need to create a CMAP for all atomnumber cases (1 to 8 solute atoms) as the scaling factor is different for each case
                if solute_in not in solute_cmap_dict:
                    logger.info(f"- Create CMAP torsions force for solute with {solute_in}/8 solute atoms /8")
                    solute_cmap_dict[solute_in] = openmm.CMAPTorsionForce()
                    solute_cmap_dict[solute_in].setName(f"CMAP_solute_{solute_in}/8")
                    self.fractional_terms.append(solute_in/8.0)
                solute_cmap_dict[solute_in].addTorsion(*cmap_indexes)
            else:
                raise ValueError("CMap not in solute or solvent")

        logger.info("Extract original CMAP maps parameters and add them to the new CMAP forces")
        self.cmap_force_map = []
        for i in range(original_cmap_force.getNumMaps()):
            map_param = original_cmap_force.getMapParameters(i)
            logger.info(f"Add CMAP Map index {i}")
            for solute_cmap_force in solute_cmap_dict.values():
                solute_cmap_force.addMap(map_param[0], map_param[1])
            solvent_cmap_force.addMap(map_param[0], map_param[1])
            self.cmap_force_map.append([map_param[0], map_param[1]])

        logger.info("Add new CMAP Forces")
        self.system.addForce(solvent_cmap_force)
        for solute_cmap_force in solute_cmap_dict.values():
            self.system.addForce(solute_cmap_force)
        self.solute_cmap_force_dict = solute_cmap_dict

        logger.info("Delete original Torsion Forces")

        # Remove the first CMAP force found in the system forces, it should be the original one as the new ones are added at the end of the forces list
        for count, force in enumerate(self.system.getForces()):
            if isinstance(force, openmm.CMAPTorsionForce):
                logger.info(f"Remove CMAP Force {count}")
                self.system.removeForce(count)
                break

    def separate_torsion_pot(self, exclude_Pro_omegas=False):
        """Use in the REST2 case as it avoids to modify
        twice the torsion terms in the rest2 system and
        in the solute system.

        Torsion potential is separated in two groups:
        - the solute (scaled one), itself split by fraction (gREST k/l scaling)
        - the solvent and not scaled solute torsion.

        As improper angles are not supposed to be scaled, here we extract only
        the proper torsion angles.

        To identify proper angles we use a trick from:
        https://github.com/maccallumlab/meld/blob/master/meld/runner/transform/rest2.py

        The original torsion potential is deleted.

        Parameters
        ----------
        exclude_Pro_omegas : bool
            The exclusion of the proline omegas scaling, default is False

        Returns
        -------
        None
        """

        energy_expression = "k*(1+cos(period*theta-phase));"

        # Create the Solvent and not-scaled solute torsion
        solvent_torsion_force = openmm.CustomTorsionForce(energy_expression)
        solvent_torsion_force.setName("Torsion_solvent")
        solvent_torsion_force.addPerTorsionParameter("period")
        solvent_torsion_force.addPerTorsionParameter("phase")
        solvent_torsion_force.addPerTorsionParameter("k")

        # Create the not scaled solute torsion
        solute_not_scaled_torsion_force = openmm.CustomTorsionForce(energy_expression)
        solute_not_scaled_torsion_force.setName("Torsion_solute_not_scaled")
        solute_not_scaled_torsion_force.addPerTorsionParameter("period")
        solute_not_scaled_torsion_force.addPerTorsionParameter("phase")
        solute_not_scaled_torsion_force.addPerTorsionParameter("k")

        # Create 4 scaled solute torsion forces, one per fraction (gREST k/l)
        # fraction 1/4, 2/4, 3/4, 4/4
        scaled_expressions = {
            frac: f"lambda_{frac}_4 * k*(1+cos(period*theta-phase));"
            for frac in [1, 2, 3, 4]
        }

        solute_scaled_torsion_forces = {}
        for frac in [1, 2, 3, 4]:
            logger.info(f"- Create scaled solute torsion force for {frac}/4 scaling: {scaled_expressions[frac]}")
            force = openmm.CustomTorsionForce(scaled_expressions[frac])
            force.setName(f"Torsion_solute_scaled_{frac}/4")
            force.addGlobalParameter(f"lambda_{frac}_4", 1.0)
            force.addPerTorsionParameter("period")
            force.addPerTorsionParameter("phase")
            force.addPerTorsionParameter("k")
            solute_scaled_torsion_forces[frac] = force

        original_torsion_force = self.system_forces["PeriodicTorsionForce"]

        bond_idxs = [sorted([i.index, j.index]) for i, j in self.topology.bonds()]

        # Identify proline backbone atoms
        if exclude_Pro_omegas:
            bond_idxs_pro = []
            for i, j in self.topology.bonds():
                if i.residue.name == "PRO" and i.name == "N" and j.name == "C":
                    bond_idxs_pro.append(sorted([i.index, j.index]))
                elif j.residue.name == "PRO" and j.name == "N" and i.name == "C":
                    bond_idxs_pro.append(sorted([i.index, j.index]))
            logger.info(f"bond_idxs_pro {bond_idxs_pro}")

        for i in range(original_torsion_force.getNumTorsions()):
            (
                p1, p2, p3, p4,
                periodicity, phase, k,
            ) = original_torsion_force.getTorsionParameters(i)

            not_improper = (
                sorted([p1, p2]) in bond_idxs
                and sorted([p2, p3]) in bond_idxs
                and sorted([p3, p4]) in bond_idxs
            )

            solute_in = (
                p1 in self.solute_index
                or p2 in self.solute_index
                or p3 in self.solute_index
                or p4 in self.solute_index
            )

            solvent_in = (
                p1 in self.solvent_index
                and p2 in self.solvent_index
                and p3 in self.solvent_index
                and p4 in self.solvent_index
            )

            not_pro_omega = True
            if solute_in and not_improper and exclude_Pro_omegas:
                if sorted([p2, p3]) in bond_idxs_pro:
                    logger.info(f"Proline omega torsion detected {p1}-{p2}-{p3}-{p4}")
                    not_pro_omega = False

            if solute_in and not_improper and not_pro_omega:
                # Count how many of the 4 atoms are in the solute (gREST k/l)
                solute_num_atom = (
                    (p1 in self.solute_index)
                    + (p2 in self.solute_index)
                    + (p3 in self.solute_index)
                    + (p4 in self.solute_index)
                )
                solute_scaled_torsion_forces[solute_num_atom].addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            elif solute_in:
                # Improper or proline omega — not scaled
                solute_not_scaled_torsion_force.addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            elif solvent_in:
                solvent_torsion_force.addTorsion(
                    p1, p2, p3, p4, [periodicity, phase, k]
                )
            else:
                raise ValueError(f"Torsion {p1}-{p2}-{p3}-{p4} not in solute or solvent")

        # Store reference to scaled forces for update
        self.solute_scaled_torsion_forces = solute_scaled_torsion_forces

        logger.info("- Add new Solvent Torsion Forces")
        if solvent_torsion_force.getNumTorsions() > 0:
            logger.info(
                f"  Adding Torsion_solvent "
                f"({solvent_torsion_force.getNumTorsions()} torsions)"
            )
            self.system.addForce(solvent_torsion_force)
        logger.info("- Add new Solute not scaled Torsion Forces")
        if solute_not_scaled_torsion_force.getNumTorsions() > 0:
            logger.info(
                f"  Adding Torsion_solute_not_scaled "
                f"({solute_not_scaled_torsion_force.getNumTorsions()} torsions)"
            )
            self.system.addForce(solute_not_scaled_torsion_force)
        self.lambda_torsion_frac = []
        for frac in [1, 2, 3, 4]:
            if solute_scaled_torsion_forces[frac].getNumTorsions() > 0:
                self.lambda_torsion_frac.append(frac)
                logger.info(
                    f"  Adding Torsion_solute_scaled_{frac}4 "
                    f"({solute_scaled_torsion_forces[frac].getNumTorsions()} torsions)"
                )
                self.system.addForce(solute_scaled_torsion_forces[frac])
                self.fractional_terms.append(frac/4.0)


        logger.info("- Delete original Torsion Forces")
        for count, force in enumerate(self.system.getForces()):
            if isinstance(force, openmm.PeriodicTorsionForce):
                self.system.removeForce(count)

    def create_solute_solvent_simulation(
        self,
        forcefield,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        platform_name="CUDA",
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        hydrogenMass=1.0 * unit.amu,
        friction=1 / unit.picoseconds,
        dt=2 * unit.femtosecond,
        temperature=300 * unit.kelvin,
    ):
        """Extract solute only and solvent only coordinates.
        A sytem and a simulation is then created for both systems.

        Parameters
        ----------
        forcefield : str
            Forcefield name
        nonbondedMethod : Nonbonded Method
            Nonbonded method, default is app.PME
        nonbondedCutoff : float * unit.nanometers
            Nonbonded cutoff
        constraints : Constraints
            Constraints
        platform_name : str
            Platform name, default is CUDA
        rigidWater : bool
            Rigid water, default is True
        ewaldErrorTolerance : float
            Ewald error tolerance, default is 0.0005
        hydrogenMass : float * unit.amu
            Hydrogen mass, default is 1.0 * unit.amu
        friction : float / unit.picoseconds
            Friction, default is 1 / unit.picoseconds
        dt : float * unit.femtosecond
            Time step, default is 2 * unit.femtosecond

        """

        # Save pdb coordinates to read them with pdb_numpy

        # Redirect stdout in the variable new_stdout:
        old_stdout = sys.stdout
        solvent_stdout = StringIO()
        solute_stdout = StringIO()
        # In case of dummy atoms (position restraints, ...)
        # It has to be removed from pdb files
        # top_num_atom = self.topology.getNumAtoms()


        solvent_top, solvent_pos = get_subset(
            self.topology, self.positions, keep=self.solvent_index, types="atom"
        )
        # print(len(solvent_pos), len([self.positions[i] for i in self.solvent_index]))
        app.PDBFile.writeFile(solvent_top, solvent_pos, solvent_stdout, True)

        # Need to use the get_subset function because of small molecule issue related
        solute_top, solute_pos = get_subset(
            self.topology, self.positions, keep=self.solute_index, types="atom"
        )

        app.PDBFile.writeFile(solute_top, solute_pos, solute_stdout, True)

        sys.stdout = old_stdout
        # Create system and simulations:
        # I don't understand why I need to pass solute_stdout (a StringIO) to StringIO(solute_stdout.getvalue())
        # It's a dirty fix
        # self.system_solute, self.simulation_solute = create_system_simulation(
        #     StringIO(solute_stdout.getvalue()),
        #     cif_format=False,
        #     forcefield=forcefield,
        #     nonbondedMethod=nonbondedMethod,
        #     nonbondedCutoff=nonbondedCutoff,
        #     constraints=constraints,
        #     platform_name=platform_name,
        #     rigidWater=rigidWater,
        #     ewaldErrorTolerance=ewaldErrorTolerance,
        #     hydrogenMass=hydrogenMass,
        #     friction=friction,
        #     dt=dt,
        #     ignoreExternalBonds=True,
        # )

        # In the case of Charmm, forcefield I need to modify the NBIX force in the system
        # before creating the simulation

        solute_integrator = openmm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            friction,
            dt
        )

        pdb_solute = app.PDBFile(StringIO(solute_stdout.getvalue()))
        print(f"Number of atoms in solute pdb: {pdb_solute.topology}")

        self.system_solute = forcefield.createSystem(
            pdb_solute.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
            ignoreExternalBonds=True,
        )
        self.system_forces_solute = {
            type(force).__name__: force for force in self.system_solute.getForces()
        }

        # In charmm forcefield, NBFIX is implemented as a CustomNonbondedForce, so we need to separate it from the solvent if it is the case and if the user want to scale NBFIX
        if self.NBFIX_flag:
            NBFIX_force = self.system_forces_solute["CustomNonbondedForce"]
            logger.info(f"Add solute flag to CustomNonbondedForce and update energy function to scale NBFIX terms")
            self.add_scale_NBFIX(NBFIX_force, range(len(self.solute_index)))
        # This is useless as the value is accessible in the original system
        # if self.lj14_flag:
        #     self.add_scale_LJ14(self.system_solute, range(len(self.solute_index)))

        self.simulation_solute = setup_simulation(
            system=self.system_solute,
            position=pdb_solute.positions,
            topology=pdb_solute.topology,
            integrator=solute_integrator,
            temperature=temperature,
            platform_name=platform_name
        )

        self.system_solvent, self.simulation_solvent = create_system_simulation(
            file_io=StringIO(solvent_stdout.getvalue()),
            cif_format=False,
            forcefield=forcefield,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            temperature=temperature,
            constraints=constraints,
            platform_name=platform_name,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
            ignoreExternalBonds=True,
        )

        self.system_forces_solvent = {
            type(force).__name__: force for force in self.system_solvent.getForces()
        }

    def create_rf_simulation(
        self,
        forcefield,
        nonbondedMethod=app.CutoffPeriodic,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        platform_name="CUDA",
        temperature=300 * unit.kelvin,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        hydrogenMass=1.0 * unit.amu,
        friction=1 / unit.picoseconds,
        dt=2 * unit.femtosecond,
    ):
        """Extract solute only and solvent only coordinates.
        A sytem and a simulation is then created for both systems.

        Parameters
        ----------
        forcefield : str
            Forcefield name
        nonbondedMethod : Nonbonded Method
            Nonbonded method, default is app.PME
        nonbondedCutoff : float * unit.nanometers
            Nonbonded cutoff
        constraints : Constraints
            Constraints
        platform_name : str
            Platform name, default is CUDA
        rigidWater : bool
            Rigid water, default is True
        ewaldErrorTolerance : float
            Ewald error tolerance, default is 0.0005
        hydrogenMass : float * unit.amu
            Hydrogen mass, default is 1.0 * unit.amu
        friction : float / unit.picoseconds
            Friction, default is 1 / unit.picoseconds
        dt : float * unit.femtosecond
            Time step, default is 2 * unit.femtosecond

        """

        # Save pdb coordinates to read them with pdb_numpy

        # Redirect stdout in the variable new_stdout:
        old_stdout = sys.stdout
        all_stdout = StringIO()
        # In case of dummy atoms (position restraints, ...)
        # It has to be removed from pdb files
        # top_num_atom = self.topology.getNumAtoms()

        all_top, all_pos = get_subset(self.topology, self.positions, types="atom")
        app.PDBFile.writeFile(all_top, all_pos, all_stdout, True)


        logger.info("- Create System with Reaction Field for non bonded electrostatic.")


        pdb = app.PDBFile(StringIO(all_stdout.getvalue()))

        integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

        self.system_rf = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=nonbondedMethod,
            nonbondedCutoff=nonbondedCutoff,
            constraints=constraints,
            rigidWater=rigidWater,
            ewaldErrorTolerance=ewaldErrorTolerance,
            hydrogenMass=hydrogenMass,
        )

        for force in self.system_rf.getForces():
            if isinstance(force, openmm.NonbondedForce):
                original_nonbonded_force = force
                break

        # Solute Solute

        custom_nonbonded_force_pp = create_custom_nonbonded_force_rf(
            original_nonbonded_force, [[self.solute_index, self.solute_index]]
        )
        custom_nonbonded_force_pp.setForceGroup(7)
        custom_nonbonded_force_pp.setName("Nonbonded_pp")
        self.nonbonded_pp_rf_force = custom_nonbonded_force_pp

        custom_bonded_force_pp, self.bond_rf_param_pp = create_custom_bonded_force_rf(
            original_nonbonded_force, [self.solute_index, self.solute_index]
        )
        logger.info(f"- Bonded rf parameter pp num: {len(self.bond_rf_param_pp)}")
        custom_bonded_force_pp.setForceGroup(8)
        custom_bonded_force_pp.setName("Bonded_pp")
        self.bonded_pp_rf_force = custom_bonded_force_pp

        self.system_rf.addForce(custom_nonbonded_force_pp)
        self.system_rf.addForce(custom_bonded_force_pp)

        # Solvent Solute

        custom_nonbonded_force_wp = create_custom_nonbonded_force_rf(
            original_nonbonded_force, [[self.solvent_index, self.solute_index]]
        )
        custom_nonbonded_force_wp.setForceGroup(9)
        custom_nonbonded_force_wp.setName("Nonbonded_wp")
        self.nonbonded_wp_rf_force = custom_nonbonded_force_wp

        custom_bonded_force_wp, self.bond_rf_param_wp = create_custom_bonded_force_rf(
            original_nonbonded_force, [self.solvent_index, self.solute_index]
        )
        logger.info(f"- Bonded rf parameter wp num: {len(self.bond_rf_param_wp)}")
        custom_bonded_force_wp.setForceGroup(10)
        custom_bonded_force_wp.setName("Bonded_wp")
        self.bonded_wp_rf_force = custom_bonded_force_wp

        self.system_rf.addForce(custom_nonbonded_force_wp)
        self.system_rf.addForce(custom_bonded_force_wp)

        
        logger.info("- Delete the original NonbondedForce.")
        for count, force in enumerate(self.system_rf.getForces()):
            if isinstance(force, openmm.NonbondedForce):
                self.system_rf.removeForce(count)
                break

        logger.info("- Delete PeriodicTorsionForce.")
        for count, force in enumerate(self.system_rf.getForces()):
            if isinstance(force, openmm.PeriodicTorsionForce):
                self.system_rf.removeForce(count)
                break

        logger.info("- Delete CMMotionRemover.")
        for count, force in enumerate(self.system_rf.getForces()):
            if isinstance(force, openmm.CMMotionRemover):
                self.system_rf.removeForce(count)
                break

        logger.info(" - Remove HarmonicBondForce element which does not concern solute")
        for count, force in enumerate(self.system_rf.getForces()):
            if isinstance(force, openmm.HarmonicBondForce):
                harmonic_bond_force = force
                break

        remove_num = 0
        for bond_index in range(harmonic_bond_force.getNumBonds()):
            p1, p2, l, k = harmonic_bond_force.getBondParameters(
                bond_index
            )
            if p1 not in self.solute_index and p2 not in self.solute_index:
                harmonic_bond_force.setBondParameters(bond_index, p1, p2, l, 0.0 * k.unit)
                remove_num += 1
                # print(f"Remove {p1}-{p2}")
            elif p1 in self.solute_index and p2 not in self.solute_index:
                logger.error(f"Bond between solute and solvent detected {p1}-{p2}")
            elif p1 not in self.solute_index and p2 in self.solute_index:
                logger.error(f"Bond between solute and solvent detected {p1}-{p2}")
        logger.info(f"   - Remove {remove_num} bonds from HarmonicBondForce")


        logger.info(" - Remove HarmonicAngleForce element which does not concern solute")
        for count, force in enumerate(self.system_rf.getForces()):
            if isinstance(force, openmm.HarmonicAngleForce):
                harmonic_angle_force = force
                break
        
        remove_num = 0
        for angle_index in range(harmonic_angle_force.getNumAngles()):
            p1, p2, p3, angle, k = harmonic_angle_force.getAngleParameters(
                angle_index
            )
            if p1 not in self.solute_index and p2 not in self.solute_index and p3 not in self.solute_index:
                harmonic_angle_force.setAngleParameters(angle_index, p1, p2, p3, angle, 0.0 * k.unit)
                remove_num += 1
                # print(f"Remove {p1}-{p2}-{p3}")
            elif p1 in self.solute_index and (p2 not in self.solute_index or p3 not in self.solute_index):
                logger.error(f"Bond between solute and solvent detected {p1}-{p2}-{p3}")
            elif p2 in self.solute_index and (p1 not in self.solute_index or p3 not in self.solute_index):
                logger.error(f"Bond between solute and solvent detected {p1}-{p2}-{p3}")
            elif p3 in self.solute_index and (p1 not in self.solute_index or p2 not in self.solute_index):
                logger.error(f"Bond between solute and solvent detected {p1}-{p2}-{p3}")
        logger.info(f"   - Remove {remove_num} angles from HarmonicAngleForce")

        # Create the simulation
        self.simulation_rf = setup_simulation(
            self.system_rf, pdb.positions, pdb.topology, integrator, platform_name
        )  

        self.system_forces_rf = {
            type(force).__name__: force for force in self.system_rf.getForces()
        }
        
        sys.stdout = old_stdout



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
            system=self.system,
            position=self.positions,
            topology=self.topology,
            temperature=temperature,
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

        sim_state = self.simulation.context.getState(
            getPositions=True,
            getEnergy=False,
            getVelocities=False,
            getForces=False,
            getParameters=False,)

        tot_positions = sim_state.getPositions(asNumpy=True)
        box_vector = sim_state.getPeriodicBoxVectors()

        if self.reaction_field:
            self.simulation_rf.context.setPeriodicBoxVectors(*box_vector)
            self.simulation_rf.context.setPositions(tot_positions)
            forces_rf = get_forces(self.system_rf, self.simulation_rf)
            forces_all = get_forces(self.system, self.simulation)

            return (None, forces_rf, forces_all)
        
        else:

            self.simulation_solute.context.setPeriodicBoxVectors(*box_vector)
            self.simulation_solute.context.setPositions(tot_positions[self.solute_index])
            forces_solute = get_forces(self.system_solute, self.simulation_solute)
            self.simulation_solvent.context.setPeriodicBoxVectors(*box_vector)
            self.simulation_solvent.context.setPositions(tot_positions[self.solvent_index])
            forces_solvent = get_forces(self.system_solvent, self.simulation_solvent)
            forces_all = get_forces(self.system, self.simulation)

            return (forces_solute, forces_solvent, forces_all)

    def update_torsion(self, lam):
        """Update torsion scaling via global lambda parameters.
        Each force group uses λ^(frac) precomputed on CPU,
        sent as a single float per group.

        Parameters
        ----------
        lam : float
            REST2 scaling parameter (T0/T_solute), between 0 and 1
        """
        context = self.simulation.context
        for frac in self.lambda_torsion_frac:
            if frac in self.solute_scaled_torsion_forces:
                context.setParameter(f"lambda_{frac}_4", lam ** (frac / 4.0))

    def update_cmap(self, scale):
        """Scale system solute cmap by a scale factor."""

        for j, map in enumerate(self.cmap_force_map):
            for atom_num, solute_cmap_force in self.solute_cmap_force_dict.items():
                scale_local = scale ** (atom_num / 8)
                solute_cmap_force.setMapParameters(j, map[0], map[1] * scale_local)

        for solute_cmap_force in self.solute_cmap_force_dict.values():
            solute_cmap_force.updateParametersInContext(self.simulation.context)


    def update_nonbonded(self, scale):
        """Scale system nonbonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        nonbonded_force = self.system_forces["NonbondedForce"]
        sqrt_scale = np.sqrt(scale)

        for i in self.solute_index:
            q, sigma, eps = self.init_nb_param[i]
            nonbonded_force.setParticleParameters(
                i, q * sqrt_scale, sigma, eps * scale
            )

        for i in range(len(self.init_nb_exept_index)):
            index = self.init_nb_exept_index[i]
            p1, p2, q, sigma, eps = self.init_nb_exept_value[i]
            nonbonded_force.setExceptionParameters(
                index, p1, p2, q * scale, sigma, eps * scale
            )
        # Need to fix simulation
        nonbonded_force.updateParametersInContext(self.simulation.context)

    def update_NBFIX(self, scale):
        """Scale system NBFIX interaction:
        - LJ epsilon by `scale`
        """
        self.simulation.context.setParameter('lambda', scale)
        self.simulation_solute.context.setParameter('lambda', scale)

    def update_lj14(self, scale):
        """Scale system 1-4 LJ interaction:"""
        if self.lj14_boundary_flag:
            self.simulation.context.setParameter("lambda_lj14_12", scale ** 0.5)
        if self.lj14_solute_flag:
            self.simulation.context.setParameter("lambda_lj14_22", scale)
        # self.simulation_solute.context.setParameter("lambda_lj14_22", scale)

    def update_nonbonded_reaction_field(self, scale):
        """Scale system nonbonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        for i in self.solute_index:
            q, sigma, eps = self.init_nb_param[i]
            self.nonbonded_pp_rf_force.setParticleParameters(
                i, [q * np.sqrt(scale), sigma, eps * scale]
            )
            self.nonbonded_wp_rf_force.setParticleParameters(
                i, [q * np.sqrt(scale), sigma, eps * scale]
            )
        self.nonbonded_pp_rf_force.updateParametersInContext(self.simulation_rf.context)
        self.nonbonded_wp_rf_force.updateParametersInContext(self.simulation_rf.context)

    def update_bonded_reaction_field(self, scale):
        """Scale system bonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        for i in range(len(self.bond_rf_param_pp)):
            p1, p2, q, sigma, eps = self.bond_rf_param_pp[i]
            self.bonded_pp_rf_force.setBondParameters(
                 i, p1, p2, [q * scale, sigma, eps * scale])
        
        self.bonded_pp_rf_force.updateParametersInContext(self.simulation_rf.context)
        
        if len(self.bond_rf_param_wp) > 0:
            for i in range(len(self.bond_rf_param_wp)):
                p1, p2, q, sigma, eps = self.bond_rf_param_wp[i]
                self.bonded_wp_rf_force.setBondParameters(
                    i, p1, p2, [q * np.sqrt(scale), sigma, eps * scale])
            
            self.bonded_wp_rf_force.updateParametersInContext(self.simulation_rf.context)

    def update_nonbonded_solute(self, scale):
        """Scale solute only system nonbonded interaction:
        - LJ epsilon by `scale`
        - Coulomb charges by `sqrt(scale)`
        - charge product is scaled by `scale`
        """

        nonbonded_force = self.system_forces_solute["NonbondedForce"]
        scale_sqrt = np.sqrt(scale)

        # assert len(self.init_nb_param) == nonbonded_force.getNumParticles()

        for i, index in enumerate(self.solute_index):
            q, sigma, eps = self.init_nb_param[index]
            nonbonded_force.setParticleParameters(
                i, q * scale_sqrt, sigma, eps * scale
            )

        for i in range(nonbonded_force.getNumExceptions()):
            p1, p2, q, sigma, eps = self.init_nb_exept_solute_value[i]
            nonbonded_force.setExceptionParameters(
                i, p1, p2, q * scale, sigma, eps * scale
            )

        nonbonded_force.updateParametersInContext(self.simulation_solute.context)

    def scale_nonbonded_torsion(self, scale):
        """Scale solute nonbonded potential and
        solute torsion potential
        """

        self.scale = scale
        if self.CMAP_flag:
            self.update_cmap(scale)
        self.update_nonbonded(scale)
        if self.reaction_field:
            self.update_nonbonded_reaction_field(scale)
            self.update_bonded_reaction_field(scale)
        else:
            self.update_nonbonded_solute(scale)
        if self.NBFIX_flag:
            self.update_NBFIX(scale)
        if self.lj14_flag:
            self.update_lj14(scale)
        self.update_torsion(scale)

    def compute_all_energies(self):
        """Extract solute potential energy and solute-solvent interactions.

        Returns
        -------
        E_frac_dict : dict {frac: unit.Quantity}
            Unscaled energy per fraction group, e.g. {0.25: ..., 0.5: ..., 1.0: ...}
        E_solute_not_scaled : unit.Quantity
            Solute energy terms not scaled in REST2 (bonds, angles, impropers)
        E_solvent : unit.Quantity
            Solvent energy
        E_solute_solvent_nb : unit.Quantity
            Unscaled solute-solvent nonbonded energy
        """

        if self.reaction_field:

            E_solute_not_scaled = 0 * unit.kilojoules_per_mole

            _, rf_force, system_force = (
                self.compute_solute_solvent_system_energy()
            )
            for i, force in rf_force.items():
                if force["name"] == "Nonbonded_wp":
                    nonbonded_rf_wp = force["energy"]
                elif force["name"] == "Bonded_wp":
                    bonded_rf_wp = force["energy"]
                elif force["name"] == "Nonbonded_pp":
                    nonbonded_rf_pp = force["energy"]
                elif force["name"] == "Bonded_pp":
                    bonded_rf_pp = force["energy"]
                elif force["name"] in ["HarmonicBondForce", "HarmonicAngleForce"]:
                    E_solute_not_scaled += force["energy"]

            E_solute_scaled = nonbonded_rf_pp + bonded_rf_pp

            solute_torsion_scaled_flag    = True
            solute_torsion_not_scaled_flag = True

            for i, force in system_force.items():
                if force["name"] == "CustomTorsionForce" and solute_torsion_scaled_flag:
                    E_solute_scaled += force["energy"]
                    solute_torsion_scaled_flag = False
                elif force["name"] == "CustomTorsionForce" and solute_torsion_not_scaled_flag:
                    E_solute_not_scaled += force["energy"]
                    solute_torsion_not_scaled_flag = False
                    break

            solvent_solute_nb = nonbonded_rf_wp + bonded_rf_wp

            # Reaction field: single fraction=1.0, return as dict
            return (
                {1.0: (1 / self.scale) * E_solute_scaled},
                E_solute_not_scaled,
                0 * unit.kilojoules_per_mole,
                (1 / self.scale) ** 0.5 * solvent_solute_nb,
            )

        solute_force, solvent_force, system_force = (
            self.compute_solute_solvent_system_energy()
        )

        def unscale(energy, frac=1.0):
            return energy / (self.scale ** frac)

        # ------------------------------------------------------------------ #
        # 1. Solute system                                                    #
        # ------------------------------------------------------------------ #

        SOLUTE_NB_SCALED = {
            "NonbondedForce",
            "LennardJones"
        }
        SOLUTE_NB_RAW    = {"LennardJones14"}
        SOLUTE_NOT_SCALED = {
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "CustomTorsionForce",   # CHARMM impropers
        }
        NOT_USED = {
            "PeriodicTorsionForce",
            "CMAPTorsionForce",
            "CMMotionRemover",
            "Total",
        }

        # Per-fraction unscaled energy accumulator
        E_frac_dict         = {}   # {frac: unscaled_energy}
        E_solute_not_scaled = 0 * unit.kilojoules_per_mole
        solute_nb           = 0 * unit.kilojoules_per_mole

        def add_frac(energy, frac):
            """Accumulate unscaled energy into the correct fraction bucket."""
            unscaled = unscale(energy, frac)
            if frac not in E_frac_dict:
                E_frac_dict[frac] = 0 * unit.kilojoules_per_mole
            E_frac_dict[frac] += unscaled

        for force in solute_force.values():
            name   = force["name"]
            energy = force["energy"]

            if name in SOLUTE_NB_SCALED:
                solute_nb += energy
                add_frac(energy, frac=1.0)
            elif name in SOLUTE_NB_RAW:
                solute_nb += energy     # subtraction only, no unscaling
            elif name in SOLUTE_NOT_SCALED:
                E_solute_not_scaled += energy
            elif name in NOT_USED:
                pass
            else:
                logger.warning(
                    f"compute_all_energies: unhandled solute force '{name}'"
                )

        # ------------------------------------------------------------------ #
        # 2. Solvent system                                                   #
        # ------------------------------------------------------------------ #

        SOLVENT_TERMS = {
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "NonbondedForce",
            "PeriodicTorsionForce",
        }
        SOLVENT_NB = {
            "NonbondedForce",
            "LennardJones"}

        E_solvent  = 0 * unit.kilojoules_per_mole
        solvent_nb = 0 * unit.kilojoules_per_mole

        for force in solvent_force.values():
            name   = force["name"]
            energy = force["energy"]

            if name in SOLVENT_NB:
                solvent_nb += energy
            if name in SOLVENT_TERMS:
                E_solvent += energy
            elif name == "CMMotionRemover":
                pass

        # ------------------------------------------------------------------ #
        # 3. Full system                                                      #
        # ------------------------------------------------------------------ #

        TORSION_SCALED_PREFIX = "Torsion_solute_scaled_"
        CMAP_SOLUTE_PREFIX    = "CMAP_solute_"

        SYSTEM_NB = {
            "NonbondedForce",
            "LennardJones"
        }

        SYSTEM_IGNORED = {
            "CMMotionRemover",
            "MonteCarloBarostat",
            "HarmonicBondForce",
            "HarmonicAngleForce",
            "CustomTorsionForce",
            "Torsion_solvent",
            "Torsion_solute_not_scaled",
            "CMAP_solvent",
            "LJ14_solvent",
            "Total",
        }

        all_nb = 0 * unit.kilojoules_per_mole

        for force in system_force.values():
            name   = force["name"]
            energy = force["energy"]

            if name in SYSTEM_NB:
                all_nb += energy

            elif name == "LJ14_solute_solute":
                add_frac(energy, frac=1.0)

            elif name == "LJ14_solute_boundary":
                add_frac(energy, frac=0.5)

            elif name.startswith(TORSION_SCALED_PREFIX):
                # Parse "Torsion_solute_scaled_3/4" → frac = 0.75
                num, den = name[len(TORSION_SCALED_PREFIX):].split("/")
                add_frac(energy, frac=int(num) / int(den))

            elif name.startswith(CMAP_SOLUTE_PREFIX):
                # Parse "CMAP_solute_7/8" → frac = 0.875
                num, den = name[len(CMAP_SOLUTE_PREFIX):].split("/")
                add_frac(energy, frac=int(num) / int(den))

            elif name == "Torsion_solute_not_scaled":
                E_solute_not_scaled += energy

            elif name in SYSTEM_IGNORED:
                pass

            else:
                logger.warning(
                    f"compute_all_energies: unhandled system force '{name}'"
                )

        # ------------------------------------------------------------------ #
        # 4. Solute-solvent nonbonded (by subtraction)                       #
        # ------------------------------------------------------------------ #

        solute_solvent_nb   = all_nb - solute_nb - solvent_nb
        E_solute_solvent_nb = unscale(solute_solvent_nb, frac=0.5)

        return (
            E_frac_dict,
            E_solute_not_scaled,
            E_solvent,
            E_solute_solvent_nb,
        )

    def get_customPotEnergie(self):
        """Extract solute potential energy and solute-solvent interactions."""

        E_solute_scaled, _, _, solvent_solute_nb = self.compute_all_energies()

        return E_solute_scaled + 0.5 * (1 / self.scale) ** 0.5 * solvent_solute_nb


def run_rest2(
    sys_rest2,
    generic_name,
    tot_steps,
    dt,
    save_step_dcd=100000,
    save_step_log=500,
    save_step_rest2=500,
    overwrite=False,
    remove_reporters=True,
    add_REST2_reporter=True,
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
    overwrite : bool, optional
        If True, overwrite previous files, by default False
    save_checkpoint_steps : int, optional
        Step to save checkpoint file, by default None

    """

    if not overwrite and os.path.isfile(generic_name + "_final.xml"):
        logger.info(
            f"File {generic_name}_final.xml exists already, skip simulate() step"
        )
        sys_rest2.simulation.loadState(generic_name + "_final.xml")
        return

    new_reporter = []
    if add_REST2_reporter:
        new_reporter = [
            Rest2Reporter(f"{generic_name}_rest2.csv", save_step_rest2, sys_rest2)
        ]

    simulate(
        sys_rest2.simulation,
        sys_rest2.topology,
        tot_steps,
        dt,
        generic_name,
        additional_reporters=new_reporter,
        save_step_log=save_step_log,
        save_step_dcd=save_step_dcd,
        remove_reporters=remove_reporters,
        save_checkpoint_steps=save_checkpoint_steps,
    )

    """

    tot_steps = np.ceil(tot_steps)
    final_step = tot_steps

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
        totalSteps=final_step,
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
        totalSteps=final_step,
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

    """


if __name__ == "__main__":
    # Check energy decomposition is correct:
    # Validates that compute_all_energies() correctly decomposes the system
    # at both 300 K (lambda=1) and a target temperature (lambda<1).
    # Uses force names (not indices) so it works across forcefields and platforms.
    import tools

    logger.setLevel(logging.INFO)
    if not logger.hasHandlers():
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # ------------------------------------------------------------------ #
    # Configuration — change these to test different setups               #
    # ------------------------------------------------------------------ #
    name = "2HPL"
    selection = "chain B"
    # selection = "(chain A and resid > 10 and resid < 21)"
    charmm_use = True
    platform_name = "OpenCL"   # OpenCL / CUDA / CPU
    nonbondedMethod = app.PME  # app.PME or app.CutoffPeriodic
    target_temperature = 600   # K — the REST2 "hot" replica temperature

    # ------------------------------------------------------------------ #
    # Force field                                                         #
    # ------------------------------------------------------------------ #
    if not charmm_use:
        forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    else:
        forcefield = app.ForceField("charmm36.xml", "charmm36/tip3p-pme-b.xml")

    dt = 2 * unit.femtosecond
    temperature = 300 * unit.kelvin
    friction = 1 / unit.picoseconds

    # ------------------------------------------------------------------ #
    # Helper: sum energies from a force-dict by name                      #
    # ------------------------------------------------------------------ #
    def sum_by_name(force_dict, names):
        total = 0.0 * unit.kilojoules_per_mole
        for f in force_dict.values():
            if f["name"] in names:
                total += f["energy"]
        return total

    def ratio(a, b):
        av = a.value_in_unit(unit.kilojoules_per_mole)
        bv = b.value_in_unit(unit.kilojoules_per_mole)
        return av / bv if bv != 0 else float("nan")

    # ------------------------------------------------------------------ #
    # Build solute / solvent PDB files                                    #
    # ------------------------------------------------------------------ #
    equi_coor = pdb_numpy.Coor(f"src/SST2/tests/inputs/{name}_equi_water.pdb")
    solute_indices = equi_coor.get_index_select(selection)
    equi_coor.select_atoms(selection).write(f"tmp_{name}_only_pep.pdb", overwrite=True)
    equi_coor.select_atoms(f"not ({selection})").write(f"tmp_{name}_no_pep.pdb", overwrite=True)

    logger.info(f"Selected {len(solute_indices)} atoms in solute group")
    logger.info(f"Selected {equi_coor.len - len(solute_indices)} atoms in solvent group")
    logger.info(f"Selected {equi_coor.len} atoms in total system")

    pdb = app.PDBFile(f"src/SST2/tests/inputs/{name}_equi_water.pdb")

    # ------------------------------------------------------------------ #
    # Classic whole-system simulation (reference)                         #
    # ------------------------------------------------------------------ #
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
    )
    simulation = setup_simulation(
        system, pdb.positions, pdb.topology,
        openmm.LangevinMiddleIntegrator(temperature, friction, dt),
        platform_name,
    )
    print("\nWhole system energy")
    tools.print_forces(system, simulation)
    forces_sys = tools.get_forces(system, simulation)

    # ------------------------------------------------------------------ #
    # Solute-only simulation (reference)                                  #
    # ------------------------------------------------------------------ #
    pdb_pep = app.PDBFile(f"tmp_{name}_only_pep.pdb")
    system_pep = forcefield.createSystem(
        pdb_pep.topology,
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        ignoreExternalBonds=True,
    )
    simulation_pep = setup_simulation(
        system_pep, pdb_pep.positions, pdb_pep.topology,
        openmm.LangevinMiddleIntegrator(temperature, friction, dt),
        platform_name,
    )
    print("\nPeptide (solute) forces:")
    tools.print_forces(system_pep, simulation_pep)
    forces_pep = tools.get_forces(system_pep, simulation_pep)

    # ------------------------------------------------------------------ #
    # Solvent-only simulation (reference)                                 #
    # ------------------------------------------------------------------ #
    pdb_no_pep = app.PDBFile(f"tmp_{name}_no_pep.pdb")
    system_no_pep = forcefield.createSystem(
        pdb_no_pep.topology,
        nonbondedMethod=nonbondedMethod,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
        ignoreExternalBonds=True,
    )
    simulation_no_pep = setup_simulation(
        system_no_pep, pdb_no_pep.positions, pdb_no_pep.topology,
        openmm.LangevinMiddleIntegrator(temperature, friction, dt),
        platform_name,
    )
    print("\nSolvent forces:")
    tools.print_forces(system_no_pep, simulation_no_pep)
    forces_no_pep = tools.get_forces(system_no_pep, simulation_no_pep)

    # Reference totals from separate simulations (force-name based)
    NB_NAMES = {"NonbondedForce"}
    BOND_NAMES = {"HarmonicBondForce"}
    ANGLE_NAMES = {"HarmonicAngleForce"}
    TORSION_NAMES = {"PeriodicTorsionForce", "CustomTorsionForce", "CMAPTorsionForce"}

    ref_pep_nb  = sum_by_name(forces_pep, NB_NAMES)
    ref_sol_nb  = sum_by_name(forces_no_pep, NB_NAMES)
    ref_sys_nb  = sum_by_name(forces_sys, NB_NAMES)
    ref_sys_tot = sum_by_name(forces_sys, BOND_NAMES | ANGLE_NAMES | TORSION_NAMES | NB_NAMES)

    ####################
    # ## REST2 setup ###
    ####################

    print(f"\nSolute indices: {len(solute_indices)} atoms")

    test = REST2(
        system,
        pdb,
        forcefield,
        solute_indices,
        openmm.LangevinMiddleIntegrator(temperature, friction, dt),
        nonbondedMethod=nonbondedMethod,
        platform_name=platform_name,
    )

    def check_decomposition(label, scale):
        print(f"\n{'='*60}")
        print(f"  {label}  (lambda = {scale:.4f})")
        print(f"{'='*60}")

        tools.print_forces(test.system, test.simulation)
        forces_rest2 = tools.get_forces(test.system, test.simulation)

        # ------------------------------------------------------------------ #
        # compute_all_energies() now returns:                                 #
        #   [0] E_frac_dict         : {frac: unscaled_energy}                #
        #   [1] E_solute_not_scaled : unit.Quantity                           #
        #   [2] E_solvent           : unit.Quantity                           #
        #   [3] E_cross_nb          : unscaled solute-solvent NB              #
        # ------------------------------------------------------------------ #
        (
            E_frac_dict,
            E_solute_not_scaled,
            E_solvent,
            E_cross_nb,
        ) = test.compute_all_energies()

        # Sum of all fractional unscaled energies (for display)
        E_solute_scaled_total = sum(
            E_frac_dict.values(), 0 * unit.kilojoules_per_mole
        )

        print(f"\ncompute_all_energies() output (unscaled physical energies):")
        print(f"  E_frac_dict:")
        for frac, E_t in sorted(E_frac_dict.items()):
            print(f"    frac={frac:.4f}  E_t = {E_t}")
        print(f"  E_solute_scaled_total = {E_solute_scaled_total}  (sum of fractions)")
        print(f"  E_solute_not_scaled   = {E_solute_not_scaled}")
        print(f"  E_solvent             = {E_solvent}")
        print(f"  E_cross_nb            = {E_cross_nb}")

        # ------------------------------------------------------------------ #
        # Reconstruct total REST2 potential energy from decomposition:        #
        # U_REST2(lambda) = sum_t lambda^f_t * E_t_unscaled                  #
        #                 + E_solute_not_scaled                               #
        #                 + E_solvent                                         #
        #                 + sqrt(lambda) * E_cross_nb                        #
        # ------------------------------------------------------------------ #
        E_rebuilt = E_solute_not_scaled + E_solvent
        for frac, E_t in E_frac_dict.items():
            E_rebuilt += (scale ** frac) * E_t
        E_rebuilt += np.sqrt(scale) * E_cross_nb

        # Potential energy from the simulation context
        E_state = test.simulation.context.getState(
            getEnergy=True
        ).getPotentialEnergy()

        print(f"\nClosure check (rebuilt vs context):")
        print(f"  E_rebuilt (from decomp) = {E_rebuilt}")
        print(f"  E_state   (from context)= {E_state}")
        print(f"  ratio rebuilt/state     = {ratio(E_rebuilt, E_state):.6f}  (expect 1.0)")

        # ------------------------------------------------------------------ #
        # At lambda=1: decomposed energies must match separate simulations   #
        # ------------------------------------------------------------------ #
        if abs(scale - 1.0) < 1e-6:
            print(f"\nCompare at lambda=1 with separate simulations:")

            # Nonbonded cross term recovered by subtraction
            rest2_nb        = sum_by_name(forces_rest2, NB_NAMES)
            cross_nb_from_sub = rest2_nb - ref_pep_nb - ref_sol_nb
            print(f"  NB cross (subtraction)   = {cross_nb_from_sub}")
            print(f"  NB cross (decomposition) = {E_cross_nb}")
            print(
                f"  ratio                    = "
                f"{ratio(E_cross_nb, cross_nb_from_sub):.6f}  (expect 1.0)"
            )

            # Total NB should match whole system
            rest2_total_nb = sum_by_name(forces_rest2, NB_NAMES)
            print(f"\n  REST2 total NB = {rest2_total_nb}")
            print(f"  System NB      = {ref_sys_nb}")
            print(
                f"  ratio          = "
                f"{ratio(rest2_total_nb, ref_sys_nb):.6f}  (expect 1.0)"
            )

            # Per-fraction check at lambda=1: scaled == unscaled
            print(f"\n  Per-fraction energies at lambda=1 (scaled == unscaled):")
            for frac, E_t in sorted(E_frac_dict.items()):
                print(f"    frac={frac:.4f}  E_t_unscaled = {E_t}")

        print(f"\nREST2 solute forces:")
        tools.print_forces(test.system_solute, test.simulation_solute)
        
    # ------------------------------------------------------------------ #
    # Check at 300 K (lambda = 1)                                         #
    # ------------------------------------------------------------------ #
    check_decomposition("300 K  (lambda=1, unscaled reference)", scale=1.0)

    # ------------------------------------------------------------------ #
    # Time individual update calls                                         #
    # ------------------------------------------------------------------ #
    scale = temperature.value_in_unit(unit.kelvin) / target_temperature  # 300/600 = 0.5

    t0 = time.time()
    test.scale_nonbonded_torsion(scale)
    print(f"\nscale_nonbonded_torsion:     {1000*(time.time()-t0):.4f}ms")

    # Reset and time each sub-update individually
    test.scale_nonbonded_torsion(1.0)

    t0 = time.time()
    test.update_nonbonded(scale)
    print(f"update_nonbonded:            {1000*(time.time()-t0):.4f}ms")

    t0 = time.time()
    test.update_torsion(scale)
    print(f"update_torsion:              {1000*(time.time()-t0):.4f}ms")

    if test.CMAP_flag:
        t0 = time.time()
        test.update_cmap(scale)
        print(f"update_cmap:                 {1000*(time.time()-t0):.4f}ms")

    if not test.reaction_field:
        t0 = time.time()
        test.update_nonbonded_solute(scale)
        print(f"update_nonbonded_solute:     {1000*(time.time()-t0):.4f}ms")

    if test.NBFIX_flag:
        t0 = time.time()
        test.update_NBFIX(scale)
        print(f"update_NBFIX:                {1000*(time.time()-t0):.4f}ms")

    if test.lj14_flag:
        t0 = time.time()
        test.update_lj14(scale)
        print(f"update_lj14:                 {1000*(time.time()-t0):.4f}ms")

    test.scale = scale  # keep internal state consistent after manual sub-updates

    # ------------------------------------------------------------------ #
    # Check at target temperature (lambda < 1)                            #
    # ------------------------------------------------------------------ #
    check_decomposition(
        f"{target_temperature} K  (lambda={scale:.4f})", scale=scale
    )

    # ------------------------------------------------------------------ #
    # Round-trip: restore lambda=1, verify energies are unchanged         #
    # ------------------------------------------------------------------ #
    print(f"\n{'='*60}")
    print("  Round-trip check: restore lambda=1")
    print(f"{'='*60}")
    test.scale_nonbonded_torsion(1.0)
    forces_back = tools.get_forces(test.system, test.simulation)
    nb_back  = sum_by_name(forces_back, NB_NAMES)
    nb_orig  = sum_by_name(tools.get_forces(system, simulation), NB_NAMES)
    print(f"  NB after round-trip = {nb_back}")
    print(f"  NB reference (orig) = {nb_orig}")
    print(f"  ratio               = {ratio(nb_back, nb_orig):.6f}  (expect 1.0)")

    """
    vmd test_2HPL/2HPL_em_water.pdb test_2HPL/2HPL_equi_water.dcd -m 2HPL.pdb
    pbc wrap -molid 0 -first 0 -last last -compound fragment -center com -centersel "chain A and protein" -orthorhombic
    """
