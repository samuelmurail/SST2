#!/usr/bin/env python3
# coding: utf-8

"""
Tests for rest2 functions
"""

import os
import pytest
import copy
import numpy as np

import openmm
from openmm import unit
import openmm.app as app

import pdb_numpy

import SST2.tools as tools
import SST2.rest2 as rest2


from .datafiles import PDB_PROT_PEP_SOL, PDB_5AWL, PDB_2HPL


def test_peptide_protein_complex_scratch(tmp_path):

    OUT_PATH = tmp_path
    name = 'tmp'

    tools.prepare_pdb(
        PDB_2HPL,
        f"{OUT_PATH}/{name}_fixed.cif",
        pH=7.0,
        overwrite=False)

    forcefield_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    forcefield = app.ForceField(*forcefield_files)

    tools.create_water_box(
        f"{OUT_PATH}/{name}_fixed.cif",
        f"{OUT_PATH}/{name}_water.cif",
        pad=1.5,
        forcefield=forcefield,
        overwrite=False)

    cif = app.PDBxFile(f"{OUT_PATH}/{name}_water.cif")

    dt = 4 * unit.femtosecond
    temperature = 310.0 * unit.kelvin
    hydrogenMass = 1.5 * unit.amu
    friction = 1.0 / unit.picoseconds

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    # Get indices of the three sets of atoms.
    all_indices = [int(i.index) for i in cif.topology.atoms()]
    solute_indices = [int(i.index) for i in cif.topology.atoms() if i.residue.chain.id in ["B"]]

    system = tools.create_sim_system(cif,
        temp=temperature,
        forcefield=forcefield,
        h_mass=hydrogenMass,
        base_force_group=1)

    sys_rest2 = rest2.REST2(
        system=system,
        pdb=cif,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator,
        temperature=temperature,
        dt=dt)


    assert sys_rest2.system.getNumForces() == 9
    assert sys_rest2.solute_index == solute_indices
    assert len(sys_rest2.solute_index) == 74
    assert len(sys_rest2.solvent_index) == pytest.approx(
        14305, 0.1
    )
    assert sys_rest2.scale == 1.0
    assert len(sys_rest2.init_nb_param) == len(sys_rest2.solute_index) + len(sys_rest2.solvent_index)
    assert len(sys_rest2.init_nb_exept_index) == 387
    assert len(sys_rest2.init_nb_exept_value) == 387
    assert len(sys_rest2.init_nb_exept_solute_value) == 387
    assert len(sys_rest2.init_torsions_index) == 206
    assert len(sys_rest2.init_torsions_value) == 206
    assert sys_rest2.solute_torsion_force.getNumTorsions() == 206


def test_peptide_protein_complex(tmp_path):
    """Test peptide protein complex"""

    tolerance = 0.0001

    # Prepare system coordinates
    prot_pep_coor = pdb_numpy.Coor(PDB_PROT_PEP_SOL)

    # Write solute and solvent coordinates
    solute_coor = prot_pep_coor.select_atoms("chain B")
    solute_coor.write(os.path.join(tmp_path, "solute.pdb"))

    solvent_coor = prot_pep_coor.select_atoms("not chain B")
    solvent_coor.write(os.path.join(tmp_path, "solvent.pdb"))

    # Set system parameters
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    dt = 2 * unit.femtosecond
    temperature = 300 * unit.kelvin
    friction = 1 / unit.picoseconds

    # Set system
    pdb = app.PDBFile(PDB_PROT_PEP_SOL)
    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
    )
    simulation = tools.setup_simulation(system, pdb.positions, pdb.topology, integrator)

    print("System Forces:")
    tools.print_forces(system, simulation)
    forces_sys = tools.get_forces(system, simulation)

    HarmonicBondForce_sys = tools.get_specific_forces(system, simulation, "HarmonicBondForce")[0]
    assert 1339.07153 == pytest.approx(
        HarmonicBondForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_sys  = tools.get_specific_forces(system, simulation, "HarmonicAngleForce")[0]
    assert 3484.75097 == pytest.approx(
        HarmonicAngleForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_sys  = tools.get_specific_forces(system, simulation, "NonbondedForce")[0]
    assert -192565.3092 == pytest.approx(
        NonbondedForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_sys  = tools.get_specific_forces(system, simulation, "PeriodicTorsionForce")[0]
    assert 5415.10546 == pytest.approx(
        PeriodicTorsionForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_sys = tools.get_specific_forces(system, simulation, "CMMotionRemover")[0]
    assert 0.0 == CMMotionRemover_sys.value_in_unit(unit.kilojoule_per_mole)

    assert -182326.38123 == pytest.approx(
        forces_sys[6]["energy"].value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    # Solute system:

    pdb_pep = app.PDBFile(os.path.join(tmp_path, "solute.pdb"))

    integrator_pep = openmm.LangevinMiddleIntegrator(temperature, friction, dt)
    system_pep = forcefield.createSystem(
        pdb_pep.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
    )
    simulation_pep = tools.setup_simulation(
        system_pep, pdb_pep.positions, pdb_pep.topology, integrator_pep
    )

    print("Solute Forces:")
    tools.print_forces(system_pep, simulation_pep)
    forces_solute = tools.get_forces(system_pep, simulation_pep)

    HarmonicBondForce_solute = tools.get_specific_forces(system_pep, simulation_pep, "HarmonicBondForce")[0]
    assert 55.969520 == pytest.approx(
        HarmonicBondForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_solute  = tools.get_specific_forces(system_pep, simulation_pep, "HarmonicAngleForce")[0]
    assert 123.60030 == pytest.approx(
        HarmonicAngleForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_solute  = tools.get_specific_forces(system_pep, simulation_pep, "NonbondedForce")[0]
    assert -451.44809 == pytest.approx(
        NonbondedForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_solute  = tools.get_specific_forces(system_pep, simulation_pep, "PeriodicTorsionForce")[0]
    assert 200.288879 == pytest.approx(
        PeriodicTorsionForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_solute = tools.get_specific_forces(system_pep, simulation_pep, "CMMotionRemover")[0]
    assert 0.0 == CMMotionRemover_solute.value_in_unit(unit.kilojoule_per_mole)

    assert -71.589387 == pytest.approx(
        forces_solute[6]["energy"].value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    pdb_no_pep = app.PDBFile(os.path.join(tmp_path, "solvent.pdb"))

    integrator_no_pep = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    system_no_pep = forcefield.createSystem(
        pdb_no_pep.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometers,
        constraints=app.HBonds,
    )

    simulation_no_pep = tools.setup_simulation(
        system_no_pep, pdb_no_pep.positions, pdb_no_pep.topology, integrator_no_pep
    )

    print("Solvent Forces:")
    tools.print_forces(system_no_pep, simulation_no_pep)
    forces_solvent = tools.get_forces(system_no_pep, simulation_no_pep)

    HarmonicBondForce_solvent = tools.get_specific_forces(system_no_pep, simulation_no_pep, "HarmonicBondForce")[0]
    assert 1283.102050 == pytest.approx(
        HarmonicBondForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_solvent = tools.get_specific_forces(system_no_pep, simulation_no_pep, "HarmonicAngleForce")[0]
    assert 3361.15087 == pytest.approx(
        HarmonicAngleForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_solvent = tools.get_specific_forces(system_no_pep, simulation_no_pep, "NonbondedForce")[0]
    assert -189514.66834 == pytest.approx(
        NonbondedForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_solvent = tools.get_specific_forces(system_no_pep, simulation_no_pep, "PeriodicTorsionForce")[0]
    assert 5214.81640 == pytest.approx(
        PeriodicTorsionForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_solvent = tools.get_specific_forces(system_no_pep, simulation_no_pep, "CMMotionRemover")[0]
    assert 0.0 == CMMotionRemover_solvent.value_in_unit(unit.kilojoule_per_mole)

    assert -179655.59901 == pytest.approx(
        forces_solvent[6]["energy"].value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    ####################
    # ## REST2 test ####
    ####################

    # Get indices of the solute atoms.
    solute_indices = [
        int(i.index) for i in pdb.topology.atoms() if i.residue.chain.id in ["B"]
    ]

    integrator_rest = openmm.LangevinMiddleIntegrator(temperature, friction, dt)
    test = rest2.REST2(system, pdb, forcefield, solute_indices, integrator_rest)

    assert test.system.getNumForces() == 8
    assert test.solute_index == solute_indices
    assert len(test.solute_index) == 74
    assert len(test.solvent_index) == 11790
    assert test.scale == 1.0
    assert len(test.init_nb_param) == len(test.solute_index) + len(test.solvent_index)
    assert len(test.init_nb_exept_index) == 387
    assert len(test.init_nb_exept_value) == 387
    assert len(test.init_nb_exept_solute_value) == 387
    assert len(test.init_torsions_index) == 206
    assert len(test.init_torsions_value) == 206
    assert test.solute_torsion_force.getNumTorsions() == 206

    # solute_scaled_torsion_force, solute_not_scaled_torsion_force, solvent_torsion_force
    torsion_len = [206, 17, 5435]
    force_i = 0

    for force in test.system.getForces():
        if isinstance(force, openmm.CustomTorsionForce):
            assert force.getNumTorsions() == torsion_len[force_i]
            force_i += 1

    print("REST2 forces 300K:")
    tools.print_forces(test.system, test.simulation)
    forces_rest2 = tools.get_forces(test.system, test.simulation)

    HarmonicBondForce_rest2 = tools.get_specific_forces(
        test.system, test.simulation, "HarmonicBondForce")[0]
    HarmonicAngleForce_rest2  = tools.get_specific_forces(
        test.system, test.simulation, "HarmonicAngleForce")[0]
    NonbondedForce_rest2  = tools.get_specific_forces(
        test.system, test.simulation, "NonbondedForce")[0]
    CustomTorsionForce_rest2  = tools.get_specific_forces(
        test.system, test.simulation, "CustomTorsionForce")
    CustomTorsionForce_rest2_sum  = np.sum(CustomTorsionForce_rest2)

    print("Compare not scaled energy rest2 vs. classic:\n")
    # HarmonicBondForce
    assert (
        pytest.approx(HarmonicBondForce_rest2 / HarmonicBondForce_sys, tolerance)
        == 1.0
    )
    # HarmonicAngleForce
    assert (
        pytest.approx(HarmonicAngleForce_rest2 / HarmonicAngleForce_sys, tolerance)
        == 1.0
    )
    # NonbondedForce
    assert (
        pytest.approx(NonbondedForce_rest2 / NonbondedForce_sys, tolerance)
        == 1.0
    )
    # CustomTorsionForce
    assert (
        pytest.approx(CustomTorsionForce_rest2_sum / PeriodicTorsionForce_sys,
            tolerance,
        )
        == 1.0
    )
    # Total
    assert (
        pytest.approx(forces_rest2[9]["energy"] / forces_sys[6]["energy"], tolerance)
        == 1.0
    )

    (
        E_solute_scaled,
        E_solute_not_scaled,
        E_solvent,
        solvent_solute_nb,
    ) = test.compute_all_energies()

    print(f"E_solute_scaled      {E_solute_scaled}")
    print(f"E_solute_not_scaled  {E_solute_not_scaled}")
    print(f"E_solute             {E_solute_scaled + E_solute_not_scaled}")
    print(f"E_solvent            {E_solvent}")
    print(f"E Solvent Solute     {solvent_solute_nb}")

    tolerance = 0.0001

    # check E_solute_scaled
    assert (
        pytest.approx(
            (CustomTorsionForce_rest2[0] + NonbondedForce_solute) / E_solute_scaled,
            tolerance,
        )
        == 1.0
    )
    # check E_solute_not_scaled
    assert (
        pytest.approx(
            (
                CustomTorsionForce_rest2[1]
                + HarmonicBondForce_solute
                + HarmonicAngleForce_solute
            )
            / E_solute_not_scaled,
            tolerance,
        )
        == 1.0
    )
    # check E_solute
    assert (
        pytest.approx(
            (E_solute_scaled + E_solute_not_scaled) / forces_solute[6]["energy"],
            tolerance,
        )
        == 1.0
    )
    # check E_solvent
    assert pytest.approx((E_solvent) / forces_solvent[6]["energy"], tolerance) == 1.0
    # check E_solvent_solute_nb
    assert (
        pytest.approx(
            (
                NonbondedForce_rest2
                - NonbondedForce_solvent
                - NonbondedForce_solute
            )
            / solvent_solute_nb,
            tolerance,
        )
        == 1.0
    )

    # Compare REST2 solute and solvent simulation with previous classic simulation
    solute_force, solvent_force = test.compute_solute_solvent_system_energy()
    print(f"Solute Force")
    for i, force in solute_force.items():
        print(
            f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol"
        )
        if force["name"] != "CMMotionRemover":
            assert (
                pytest.approx(force["energy"] / forces_solute[i]["energy"], tolerance)
                == 1.0
            )

    print(f"Solvent Force")
    for i, force in solvent_force.items():
        print(
            f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol"
        )
        if force["name"] != "CMMotionRemover":
            assert (
                pytest.approx(force["energy"] / forces_solvent[i]["energy"], tolerance)
                == 1.0
            )

    tolerance = 0.0001

    for scale in [0.5, 1.0, 1.5]:
        test.scale_nonbonded_torsion(scale)
        assert test.scale == scale
        print(f"\nREST2 forces lambda = {scale:.1f}  Temp  = {300/scale:.1f} K\n")
        forces_rest2_new = tools.get_forces(test.system, test.simulation)

        # compare scaled with non scaled REST2:
        for i, force in forces_rest2_new.items():
            print(
                f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol"
            )
            if scale != 1.0:
                if force['name'] in ['NonbondedForce', 'Total']:  # NonbondedForce, Total
                    assert (
                        pytest.approx(
                            force["energy"] / forces_rest2[i]["energy"], tolerance
                        )
                        != 1.0
                    )
                elif force['name'] == "CustomTorsionForce" and i in [4]:  # CustomTorsionForce scaled
                    assert (
                        pytest.approx(
                            (force["energy"] / scale) / forces_rest2[i]["energy"],
                            tolerance,
                        )
                        == 1.0
                    )
                elif force['name'] in ['HarmonicBondForce', 'HarmonicAngleForce'] or i in [
                    5,
                    6,
                ]:  # HarmonicBondForce, HarmonicAngleForce, CustomTorsionForces not scaled
                    assert (
                        pytest.approx(
                            force["energy"] / forces_rest2[i]["energy"], tolerance
                        )
                        == 1.0
                    )
            elif force['name'] in ['HarmonicBondForce', 'HarmonicAngleForce', 'NonbondedForce', "CustomTorsionForce", 'Total']:
                assert (
                    pytest.approx(
                        (force["energy"] / scale) / forces_rest2[i]["energy"], tolerance
                    )
                    == 1.0
                )

        (
            E_solute_scaled_new,
            E_solute_not_scaled_new,
            E_solvent_new,
            solvent_solute_nb_new,
        ) = test.compute_all_energies()

        print(f"E_solute_scaled      {E_solute_scaled_new}")
        print(f"E_solute_not_scaled  {E_solute_not_scaled_new}")
        print(f"E_solvent            {E_solvent_new}")
        print(f"solvent_solute_nb    {solvent_solute_nb_new}")

        assert pytest.approx(E_solute_scaled_new / E_solute_scaled, tolerance) == 1.0
        assert (
            pytest.approx(E_solute_not_scaled_new / E_solute_not_scaled, tolerance)
            == 1.0
        )
        assert pytest.approx(E_solvent_new / E_solvent_new, tolerance) == 1.0
        assert (
            pytest.approx(solvent_solute_nb_new / solvent_solute_nb, tolerance) == 1.0
        )

        # Compare REST2 solute and solvent simulation with previous classic simulation
        (
            solute_force_new,
            solvent_force_new,
        ) = test.compute_solute_solvent_system_energy()
        print(f"Solute Force")
        for i, force in solute_force_new.items():
            print(
                f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol {forces_solute[i]['energy']}"
            )
            if force['name'] in ['HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce']:  # HarmonicBondForce, HarmonicAngleForce, CustomTorsionForces not scaled
                assert (
                    pytest.approx(
                        force["energy"] / forces_solute[i]["energy"], tolerance
                    )
                    == 1.0
                )
            if force['name'] in ['NonbondedForce']:  # NonbondedForce
                assert (
                    pytest.approx(
                        (force["energy"] / scale) / forces_solute[i]["energy"],
                        tolerance,
                    )
                    == 1.0
                )

        print(f"Solvent Force")
        for i, force in solvent_force_new.items():
            print(
                f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol"
            )
            if force["name"] != "CMMotionRemover":
                assert (
                    pytest.approx(
                        force["energy"] / forces_solvent[i]["energy"], tolerance
                    )
                    == 1.0
                )

        # check E_solute_scaled
        print(solute_force_new)
        for i in solute_force_new:
            force = solute_force_new[i]
            if force["name"] == "NonbondedForce":
                solute_nonbonded_new = force["energy"]
                break
        assert (
            pytest.approx(
                (
                    (forces_rest2_new[4]["energy"] + solute_nonbonded_new)
                    / scale
                )
                / E_solute_scaled_new,
                tolerance,
            )
            == 1.0
        )
        # check E_solute_not_scaled
        solute_nonbonded_new = 0 * unit.kilojoule_per_mole
        for i in solute_force_new:
            force = solute_force_new[i]
            if force["name"] in ["HarmonicBondForce", "HarmonicAngleForce"]:
                solute_nonbonded_new += force["energy"]
        assert (
            pytest.approx(
                (
                    forces_rest2_new[5]["energy"]
                    + solute_nonbonded_new
                )
                / E_solute_not_scaled_new,
                tolerance,
            )
            == 1.0
        )
        # check E_solute
        Correct_E_solute = 0 * unit.kilojoule_per_mole
        for i in forces_rest2_new:
            force = forces_rest2_new[i]
            if i in [4,5]:
                Correct_E_solute += force["energy"]

        for i in solute_force_new:
            force = solute_force_new[i]
            if force["name"] in ["HarmonicBondForce", "HarmonicAngleForce", "NonbondedForce"]:
                Correct_E_solute += force["energy"]

        assert (
            pytest.approx(
                (E_solute_scaled_new * scale + E_solute_not_scaled_new)
                / Correct_E_solute,
                tolerance,
            )
            == 1.0
        )
        # check E_solvent
        assert (
            pytest.approx((E_solvent_new) / solvent_force_new[6]["energy"], tolerance)
            == 1.0
        )
        # check E_solvent_solute_nb
        for i in solvent_force_new:
            force = solvent_force_new[i]
            if force["name"] == "NonbondedForce":
                solvent_nonbonded_new = force["energy"]
                break
        for i in solute_force_new:
            force = solute_force_new[i]
            if force["name"] == "NonbondedForce":
                solute_nonbonded_new = force["energy"]
                break
        for i in forces_rest2_new:
            force = forces_rest2_new[i]
            if force["name"] == "NonbondedForce":
                rest2_nonbonded_new = force["energy"]
                break
        assert (
            pytest.approx(
                (
                    rest2_nonbonded_new
                    - solute_nonbonded_new
                    - solvent_nonbonded_new
                )
                / (scale**0.5 * solvent_solute_nb_new),
                tolerance,
            )
            == 1.0
        )



def test_5awl_omega_PRO(tmp_path):
    """Test peptide protein complex"""

    tolerance = 0.00001

 
    name = "5awl"
    OUT_PATH = tmp_path

    if not os.path.exists(OUT_PATH):
        os.makedirs(OUT_PATH)

    tools.prepare_pdb(
        PDB_5AWL, os.path.join(tmp_path, f"{name}_fixed.cif"), pH=7.0, overwrite=False
    )

    forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    forcefield = app.ForceField(*forcefield_files)

    tools.create_water_box(
        os.path.join(tmp_path, f"{name}_fixed.cif"),
        os.path.join(tmp_path, f"{name}_water.cif"),
        pad=1.5,
        forcefield=forcefield,
        overwrite=False,
    )

    dt = 4 * unit.femtosecond
    temperature = 300 * unit.kelvin
    friction = 1 / unit.picoseconds
    hydrogenMass = 3 * unit.amu
    rigidWater = True
    ewaldErrorTolerance = 0.0005
    nsteps = 0.01 * unit.nanoseconds / dt

    pdb = app.PDBxFile(os.path.join(tmp_path, f"{name}_water.cif"))

    # Get indices of the three sets of atoms.
    all_indices = [int(i.index) for i in pdb.topology.atoms()]
    solute_indices = [
        int(i.index) for i in pdb.topology.atoms() if i.residue.chain.id in ["A"]
    ]

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    system = tools.create_sim_system(
        pdb,
        temp=temperature,
        forcefield=forcefield,
        h_mass=hydrogenMass,
        base_force_group=1,
    )

    system_2 = copy.deepcopy(system)

    test = rest2.REST2(
        system=system,
        pdb=pdb,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator,
        dt=dt,
    )

    assert test.system.getNumForces() == 9
    assert test.solute_index == solute_indices
    assert len(test.solute_index) == 166
    assert len(test.solvent_index) == pytest.approx(4000, abs=2000)
    assert test.scale == 1.0
    assert len(test.init_nb_param) == len(test.solute_index) + len(test.solvent_index)
    assert len(test.init_nb_exept_index) == 899
    assert len(test.init_nb_exept_value) == 899
    assert len(test.init_nb_exept_solute_value) == 899
    assert len(test.init_torsions_index) == 521
    assert len(test.init_torsions_value) == 521
    assert test.solute_torsion_force.getNumTorsions() == 521

    torsion_len = [521, 46, 0]
    force_i = 0

    for force in test.system.getForces():
        if isinstance(force, openmm.CustomTorsionForce):
            assert force.getNumTorsions() == torsion_len[force_i]
            force_i += 1


    integrator_2 = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    test_2 = rest2.REST2(
        system=system_2,
        pdb=pdb,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator_2,
        dt=dt,
        exclude_Pro_omegas=True,
    )

    assert test_2.system.getNumForces() == 9
    assert test_2.solute_index == solute_indices
    assert len(test_2.solute_index) == 166
    assert len(test_2.solvent_index) ==  pytest.approx(4000, abs=2000)
    assert test_2.scale == 1.0
    assert len(test_2.init_nb_param) == len(test_2.solute_index) + len(test_2.solvent_index)
    assert len(test_2.init_nb_exept_index) == 899
    assert len(test_2.init_nb_exept_value) == 899
    assert len(test_2.init_nb_exept_solute_value) == 899
    assert len(test_2.init_torsions_index) == 521 - 4
    assert len(test_2.init_torsions_value) == 521 - 4
    assert test_2.solute_torsion_force.getNumTorsions() == 521 - 4

    torsion_len = [521-4, 46+4, 0]
    force_i = 0

    for force in test_2.system.getForces():
        if isinstance(force, openmm.CustomTorsionForce):
            assert force.getNumTorsions() == torsion_len[force_i]
            force_i += 1

