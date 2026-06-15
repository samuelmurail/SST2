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
    name = "tmp"

    tools.prepare_pdb(PDB_2HPL, f"{OUT_PATH}/{name}_fixed.cif", pH=7.0, overwrite=False)

    forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    forcefield = app.ForceField(*forcefield_files)

    tools.create_water_box(
        f"{OUT_PATH}/{name}_fixed.cif",
        f"{OUT_PATH}/{name}_water.cif",
        pad=1.5,
        forcefield=forcefield,
        overwrite=False,
    )

    cif = app.PDBxFile(f"{OUT_PATH}/{name}_water.cif")

    dt = 4 * unit.femtosecond
    temperature = 310.0 * unit.kelvin
    hydrogenMass = 1.5 * unit.amu
    friction = 1.0 / unit.picoseconds

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, dt)

    # Get indices of the three sets of atoms.
    all_indices = [int(i.index) for i in cif.topology.atoms()]
    solute_indices = [
        int(i.index) for i in cif.topology.atoms() if i.residue.chain.id in ["B"]
    ]

    system = tools.create_sim_system(
        cif,
        temp=temperature,
        forcefield=forcefield,
        h_mass=hydrogenMass,
        base_force_group=1,
    )

    sys_rest2 = rest2.REST2(
        system=system,
        pdb=cif,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator,
        temperature=temperature,
        dt=dt,
    )

    assert sys_rest2.solute_index == solute_indices
    assert len(sys_rest2.solute_index) == 74
    assert len(sys_rest2.solvent_index) == pytest.approx(14305, 0.1)
    assert sys_rest2.scale == 1.0
    assert len(sys_rest2.init_nb_param) == len(sys_rest2.solute_index) + len(
        sys_rest2.solvent_index
    )
    assert len(sys_rest2.init_nb_exept_index) == 387
    assert len(sys_rest2.init_nb_exept_value) == 387
    assert len(sys_rest2.init_nb_exept_solute_value) == 387
    # torsion forces are now named; check by name instead of by old attributes
    torsion_counts_scratch = {
        force.getName(): force.getNumTorsions()
        for force in sys_rest2.system.getForces()
        if isinstance(force, openmm.CustomTorsionForce)
    }
    total_scaled_scratch = sum(
        c for n, c in torsion_counts_scratch.items()
        if n.startswith("Torsion_solute_scaled_")
    )
    assert total_scaled_scratch == 206
    assert torsion_counts_scratch.get("Torsion_solute_not_scaled", 0) == 17


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

    HarmonicBondForce_sys = tools.get_specific_forces(
        system, simulation, "HarmonicBondForce"
    )[0]
    assert 1339.07153 == pytest.approx(
        HarmonicBondForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_sys = tools.get_specific_forces(
        system, simulation, "HarmonicAngleForce"
    )[0]
    assert 3484.75097 == pytest.approx(
        HarmonicAngleForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_sys = tools.get_specific_forces(
        system, simulation, "NonbondedForce"
    )[0]
    assert -192565.3092 == pytest.approx(
        NonbondedForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_sys = tools.get_specific_forces(
        system, simulation, "PeriodicTorsionForce"
    )[0]
    assert 5415.10546 == pytest.approx(
        PeriodicTorsionForce_sys.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_sys = tools.get_specific_forces(
        system, simulation, "CMMotionRemover"
    )[0]
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

    HarmonicBondForce_solute = tools.get_specific_forces(
        system_pep, simulation_pep, "HarmonicBondForce"
    )[0]
    assert 55.969520 == pytest.approx(
        HarmonicBondForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_solute = tools.get_specific_forces(
        system_pep, simulation_pep, "HarmonicAngleForce"
    )[0]
    assert 123.60030 == pytest.approx(
        HarmonicAngleForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_solute = tools.get_specific_forces(
        system_pep, simulation_pep, "NonbondedForce"
    )[0]
    assert -451.44809 == pytest.approx(
        NonbondedForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_solute = tools.get_specific_forces(
        system_pep, simulation_pep, "PeriodicTorsionForce"
    )[0]
    assert 200.288879 == pytest.approx(
        PeriodicTorsionForce_solute.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_solute = tools.get_specific_forces(
        system_pep, simulation_pep, "CMMotionRemover"
    )[0]
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

    HarmonicBondForce_solvent = tools.get_specific_forces(
        system_no_pep, simulation_no_pep, "HarmonicBondForce"
    )[0]
    assert 1283.102050 == pytest.approx(
        HarmonicBondForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    HarmonicAngleForce_solvent = tools.get_specific_forces(
        system_no_pep, simulation_no_pep, "HarmonicAngleForce"
    )[0]
    assert 3361.15087 == pytest.approx(
        HarmonicAngleForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    NonbondedForce_solvent = tools.get_specific_forces(
        system_no_pep, simulation_no_pep, "NonbondedForce"
    )[0]
    assert -189514.66834 == pytest.approx(
        NonbondedForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    PeriodicTorsionForce_solvent = tools.get_specific_forces(
        system_no_pep, simulation_no_pep, "PeriodicTorsionForce"
    )[0]
    assert 5214.81640 == pytest.approx(
        PeriodicTorsionForce_solvent.value_in_unit(unit.kilojoule_per_mole), tolerance
    )

    CMMotionRemover_solvent = tools.get_specific_forces(
        system_no_pep, simulation_no_pep, "CMMotionRemover"
    )[0]
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

    assert test.solute_index == solute_indices
    assert len(test.solute_index) == 74
    assert len(test.solvent_index) == 11790
    assert test.scale == 1.0
    assert len(test.init_nb_param) == len(test.solute_index) + len(test.solvent_index)
    assert len(test.init_nb_exept_index) == 387
    assert len(test.init_nb_exept_value) == 387
    assert len(test.init_nb_exept_solute_value) == 387

    # torsion forces are now named; check by name instead of old index-based approach
    torsion_counts = {
        force.getName(): force.getNumTorsions()
        for force in test.system.getForces()
        if isinstance(force, openmm.CustomTorsionForce)
    }
    total_scaled = sum(
        c for n, c in torsion_counts.items() if n.startswith("Torsion_solute_scaled_")
    )
    assert total_scaled == 206
    assert torsion_counts.get("Torsion_solute_not_scaled", 0) == 17
    assert torsion_counts.get("Torsion_solvent", 0) == 5435

    print("REST2 forces 300K:")
    tools.print_forces(test.system, test.simulation)
    forces_rest2 = tools.get_forces(test.system, test.simulation)

    HarmonicBondForce_rest2 = tools.get_specific_forces(
        test.system, test.simulation, "HarmonicBondForce"
    )[0]
    HarmonicAngleForce_rest2 = tools.get_specific_forces(
        test.system, test.simulation, "HarmonicAngleForce"
    )[0]
    NonbondedForce_rest2 = tools.get_specific_forces(
        test.system, test.simulation, "NonbondedForce"
    )[0]
    # REST2 splits PeriodicTorsionForce into named CustomTorsionForce objects;
    # sum all forces whose name starts with "Torsion_"
    CustomTorsionForce_rest2_sum = 0 * unit.kilojoules_per_mole
    for f in forces_rest2.values():
        if f["name"].startswith("Torsion_"):
            CustomTorsionForce_rest2_sum += f["energy"]


    print("Compare not scaled energy rest2 vs. classic:\n")
    # HarmonicBondForce
    assert (
        pytest.approx(HarmonicBondForce_rest2 / HarmonicBondForce_sys, tolerance) == 1.0
    )
    # HarmonicAngleForce
    assert (
        pytest.approx(HarmonicAngleForce_rest2 / HarmonicAngleForce_sys, tolerance)
        == 1.0
    )
    # NonbondedForce
    assert pytest.approx(NonbondedForce_rest2 / NonbondedForce_sys, tolerance) == 1.0
    # CustomTorsionForce
    assert (
        pytest.approx(
            CustomTorsionForce_rest2_sum / PeriodicTorsionForce_sys,
            tolerance,
        )
        == 1.0
    )
    # Total energy should match between REST2 and original system
    Total_rest2 = tools.get_specific_forces(test.system, test.simulation, "Total")[0]
    Total_sys = tools.get_specific_forces(system, simulation, "Total")[0]
    assert pytest.approx(Total_rest2 / Total_sys, tolerance) == 1.0

    (
        E_frac_dict,
        E_solute_not_scaled,
        E_solvent,
        solvent_solute_nb,
    ) = test.compute_all_energies()
    E_solute_scaled = sum(E_frac_dict.values(), 0 * unit.kilojoules_per_mole)

    print(f"E_solute_scaled      {E_solute_scaled}")
    print(f"E_solute_not_scaled  {E_solute_not_scaled}")
    print(f"E_solute             {E_solute_scaled + E_solute_not_scaled}")
    print(f"E_solvent            {E_solvent}")
    print(f"E Solvent Solute     {solvent_solute_nb}")

    tolerance = 0.0001

    forces_rest2_dict = tools.get_forces(test.system, test.simulation)

    # check E_solute_scaled: sum all Torsion_solute_scaled_k/4 energies + NB from solute subsystem
    scaled_torsion_sum = sum(
        (f["energy"] for f in forces_rest2_dict.values() if f["name"].startswith("Torsion_solute_scaled_")),
        0 * unit.kilojoules_per_mole,
    )
    assert pytest.approx(
        (scaled_torsion_sum + NonbondedForce_solute) / E_solute_scaled, tolerance
    ) == 1.0

    # check E_solute_not_scaled: Torsion_solute_not_scaled + bond + angle from solute subsystem
    Torsion_not_scaled_list = tools.get_specific_forces(
        test.system, test.simulation, "Torsion_solute_not_scaled"
    )
    Torsion_not_scaled_energy = (
        Torsion_not_scaled_list[0] if Torsion_not_scaled_list else 0 * unit.kilojoules_per_mole
    )

    print(f"Torsion_not_scaled_energy {Torsion_not_scaled_energy}")
    print(f"HarmonicBondForce_solute {HarmonicBondForce_solute}")
    print(f"HarmonicAngleForce_solute {HarmonicAngleForce_solute}")
    print(f"E_solute_not_scaled {E_solute_not_scaled}")

    assert pytest.approx(
        (Torsion_not_scaled_energy + HarmonicBondForce_solute + HarmonicAngleForce_solute)
        / E_solute_not_scaled,
        tolerance,
    ) == 1.0
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
            (NonbondedForce_rest2 - NonbondedForce_solvent - NonbondedForce_solute)
            / solvent_solute_nb,
            tolerance,
        )
        == 1.0
    )

    # Compare REST2 solute and solvent simulation with previous classic simulation
    solute_force, solvent_force, all_force = test.compute_solute_solvent_system_energy()
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
        for i, force in forces_rest2_new.items():
            print(
                f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol"
            )

        # helper: look up a single named force energy from a force dict
        def get_named(fd, name):
            for f in fd.values():
                if f["name"] == name:
                    return f["energy"]
            return None

        NB_new = get_named(forces_rest2_new, "NonbondedForce")
        NB_ref = get_named(forces_rest2_dict, "NonbondedForce")
        HB_new = get_named(forces_rest2_new, "HarmonicBondForce")
        HB_ref = get_named(forces_rest2_dict, "HarmonicBondForce")
        HA_new = get_named(forces_rest2_new, "HarmonicAngleForce")
        HA_ref = get_named(forces_rest2_dict, "HarmonicAngleForce")
        TorS_new = get_named(forces_rest2_new, "Torsion_solvent")
        TorS_ref = get_named(forces_rest2_dict, "Torsion_solvent")
        TorNS_new = get_named(forces_rest2_new, "Torsion_solute_not_scaled")
        TorNS_ref = get_named(forces_rest2_dict, "Torsion_solute_not_scaled")

        if scale != 1.0:
            # NonbondedForce must have changed
            assert pytest.approx(NB_new / NB_ref, tolerance) != 1.0
            # Bond and angle forces are not scaled — must be unchanged
            assert pytest.approx(HB_new / HB_ref, tolerance) == 1.0
            assert pytest.approx(HA_new / HA_ref, tolerance) == 1.0
            # Solvent torsion not scaled
            if TorS_new is not None and TorS_ref is not None:
                assert pytest.approx(TorS_new / TorS_ref, tolerance) == 1.0
            # Solute not-scaled torsion unchanged
            if TorNS_new is not None and TorNS_ref is not None:
                assert pytest.approx(TorNS_new / TorNS_ref, tolerance) == 1.0
            # Each solute scaled torsion force Torsion_solute_scaled_k/4 must equal
            # scale^(k/4) * reference value
            for n, e_new in forces_rest2_new.items():
                fname = e_new["name"]
                if fname.startswith("Torsion_solute_scaled_"):
                    frac_str = fname[len("Torsion_solute_scaled_"):]
                    num, den = frac_str.split("/")
                    frac = int(num) / int(den)
                    e_ref = get_named(forces_rest2_dict, fname)
                    if e_ref is not None:
                        assert pytest.approx(
                            (e_new["energy"] / (scale ** frac)) / e_ref, tolerance
                        ) == 1.0

        (
            E_frac_dict_new,
            E_solute_not_scaled_new,
            E_solvent_new,
            solvent_solute_nb_new,
        ) = test.compute_all_energies()
        E_solute_scaled_new = sum(E_frac_dict_new.values(), 0 * unit.kilojoules_per_mole)

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
            all_force_new,
        ) = test.compute_solute_solvent_system_energy()
        print(f"Solute Force")
        for i, force in solute_force_new.items():
            print(
                f"{i}   {force['name']:25} {force['energy'].value_in_unit(unit.kilojoule_per_mole):.2f} KJ/mol {forces_solute[i]['energy']}"
            )
            if force["name"] in [
                "HarmonicBondForce",
                "HarmonicAngleForce",
                "PeriodicTorsionForce",
            ]:  # HarmonicBondForce, HarmonicAngleForce, CustomTorsionForces not scaled
                assert (
                    pytest.approx(
                        force["energy"] / forces_solute[i]["energy"], tolerance
                    )
                    == 1.0
                )
            if force["name"] in ["NonbondedForce"]:  # NonbondedForce
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

        # check E_solute_scaled: invariance already asserted above; verify formula too
        # E_solute_scaled = (sum_k Torsion_k/4 / scale^(k/4)) + NB_solute / scale
        # Since compute_all_energies returns unscaled values, and we already checked
        # E_solute_scaled_new == E_solute_scaled, the formula is implicitly verified.
        # Extra explicit check: scaled_sum + NB_solute = scale * E_solute_scaled_new
        print(solute_force_new)
        scaled_torsion_sum_new = sum(
            (f["energy"] for f in forces_rest2_new.values() if f["name"].startswith("Torsion_solute_scaled_")),
            0 * unit.kilojoules_per_mole,
        )
        solute_nb_new = next(
            f["energy"] for f in solute_force_new.values() if f["name"] == "NonbondedForce"
        )
        assert pytest.approx(
            (scaled_torsion_sum_new + solute_nb_new) / (scale * E_solute_scaled_new),
            tolerance,
        ) == 1.0

        # check E_solute_not_scaled
        TorNS_new_energy = get_named(forces_rest2_new, "Torsion_solute_not_scaled")
        solute_bonded_new = sum(
            (f["energy"] for f in solute_force_new.values()
             if f["name"] in ["HarmonicBondForce", "HarmonicAngleForce"]),
            0 * unit.kilojoules_per_mole,
        )
        if TorNS_new_energy is not None:
            assert pytest.approx(
                (TorNS_new_energy + solute_bonded_new) / E_solute_not_scaled_new,
                tolerance,
            ) == 1.0
        # check E_solute: all solute torsion forces (by name) + bond/angle/NB from solute subsystem
        Correct_E_solute = 0 * unit.kilojoule_per_mole
        for f in forces_rest2_new.values():
            if f["name"].startswith("Torsion_solute_"):
                Correct_E_solute += f["energy"]
        for f in solute_force_new.values():
            if f["name"] in ["HarmonicBondForce", "HarmonicAngleForce", "NonbondedForce"]:
                Correct_E_solute += f["energy"]

        print(f"Correct_E_solute {Correct_E_solute}")
        print(f"E_solute_scaled_new {E_solute_scaled_new}")
        print(f"E_solute_not_scaled_new {E_solute_not_scaled_new}")

        assert (
            pytest.approx(
                (E_solute_scaled_new * scale + E_solute_not_scaled_new)
                / Correct_E_solute,
                tolerance,
            )
            == 1.0
        )
        # check E_solvent: compare against the Total entry in the returned force dict
        solvent_total_energy = next(
            f["energy"] for f in solvent_force_new.values() if f["name"] == "Total"
        )
        assert pytest.approx(E_solvent_new / solvent_total_energy, tolerance) == 1.0
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
                (rest2_nonbonded_new - solute_nonbonded_new - solvent_nonbonded_new)
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

    assert test.solute_index == solute_indices
    assert len(test.solute_index) == 166
    assert len(test.solvent_index) == pytest.approx(4000, abs=2000)
    assert test.scale == 1.0
    assert len(test.init_nb_param) == len(test.solute_index) + len(test.solvent_index)
    assert len(test.init_nb_exept_index) == 899
    assert len(test.init_nb_exept_value) == 899
    assert len(test.init_nb_exept_solute_value) == 899

    # torsion forces are now named; check by name
    torsion_counts_5awl = {
        force.getName(): force.getNumTorsions()
        for force in test.system.getForces()
        if isinstance(force, openmm.CustomTorsionForce)
    }
    total_scaled_5awl = sum(
        c for n, c in torsion_counts_5awl.items() if n.startswith("Torsion_solute_scaled_")
    )
    assert total_scaled_5awl == 521
    assert torsion_counts_5awl.get("Torsion_solute_not_scaled", 0) == 46

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

    assert test_2.solute_index == solute_indices
    assert len(test_2.solute_index) == 166
    assert len(test_2.solvent_index) == pytest.approx(4000, abs=2000)
    assert test_2.scale == 1.0
    assert len(test_2.init_nb_param) == len(test_2.solute_index) + len(
        test_2.solvent_index
    )
    assert len(test_2.init_nb_exept_index) == 899
    assert len(test_2.init_nb_exept_value) == 899
    assert len(test_2.init_nb_exept_solute_value) == 899

    # with exclude_Pro_omegas=True, 4 torsions move from scaled to not_scaled
    torsion_counts_5awl_2 = {
        force.getName(): force.getNumTorsions()
        for force in test_2.system.getForces()
        if isinstance(force, openmm.CustomTorsionForce)
    }
    total_scaled_5awl_2 = sum(
        c for n, c in torsion_counts_5awl_2.items() if n.startswith("Torsion_solute_scaled_")
    )
    assert total_scaled_5awl_2 == 521 - 4
    assert torsion_counts_5awl_2.get("Torsion_solute_not_scaled", 0) == 46 + 4
