#!/usr/bin/env python3
# coding: utf-8

"""
Tests for rest2 functions
"""

import os
import pandas as pd
import pytest

import openmm
from openmm import unit
import openmm.app as app


import SST2.tools as tools
import SST2.rest2 as rest2
import SST2.sst2 as sst2

from .datafiles import PDB_5AWL


def test_peptide_protein_complex(tmp_path):
    """Test SST2 protocol"""

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

    sys_rest2 = rest2.REST2(
        system=system,
        pdb=pdb,
        forcefield=forcefield,
        solute_index=solute_indices,
        integrator=integrator,
        dt=dt,
        platform_name=tools.get_fastest_platform_name(),
    )

    tools.minimize(
        sys_rest2.simulation,
        os.path.join(tmp_path, f"{name}_em_water.cif"),
        pdb.topology,
        maxIterations=1000,
        overwrite=False,
    )

    sys_rest2.simulation.context.setVelocitiesToTemperature(temperature)

    rest2.run_rest2(
        sys_rest2,
        os.path.join(tmp_path, f"{name}_equi_water"),
        tot_steps=nsteps,
        dt=dt,
        save_step_dcd=100000,
        save_step_log=10000,
        save_step_rest2=500,
    )

    ladder_num = tools.compute_ladder_num(
        os.path.join(tmp_path, f"{name}_equi_water_rest2"), 300, 500, sst2_score=True
    )

    assert ladder_num == 3

    temp_list = tools.compute_temperature_list(
        minTemperature=300,
        maxTemperature=500,
        numTemperatures=ladder_num,
        refTemperature=300,
    )

    temp_expected = [
        300.0 * unit.kelvin,
        387.2983346207417 * unit.kelvin,
        500.0 * unit.kelvin,
    ]

    assert len(temp_expected) == len(temp_list) == ladder_num

    for i in range(len(temp_list)):
        assert temp_list[i] == temp_expected[i]

    tot_steps = 10000
    save_step_dcd = 1000
    save_step_log = 100
    tempChangeInterval = int(2.0 * unit.picosecond / dt.in_units_of(unit.picosecond))
    save_check_steps = 5000

    sst2.run_sst2(
        sys_rest2,
        os.path.join(tmp_path, f"{name}"),
        tot_steps,
        dt=dt,
        temperatures=temp_list,
        ref_temp=temp_list[1],
        save_step_dcd=save_step_dcd,
        save_step_log=save_step_log,
        save_step_rest2=save_step_log,
        tempChangeInterval=tempChangeInterval,
        reportInterval=save_step_log,
        overwrite=False,
        save_checkpoint_steps=save_check_steps,
    )

    print(tmp_path)

    # Check if the output files are created
    assert os.path.exists(os.path.join(tmp_path, f"{name}_sst2.dcd"))
    assert os.path.exists(os.path.join(tmp_path, f"{name}_sst2.pdb"))
    assert os.path.exists(os.path.join(tmp_path, f"{name}_sst2.cif"))
    assert os.path.exists(os.path.join(tmp_path, f"{name}_sst2.xml"))

    file_path = os.path.join(tmp_path, f"{name}_sst2.csv")
    assert os.path.exists(file_path)
    num_lines = sum(1 for _ in open(file_path))
    assert num_lines == tot_steps / save_step_log + 1

    file_path = os.path.join(tmp_path, f"{name}_sst2_full.csv")
    assert os.path.exists(file_path)
    num_lines = sum(1 for _ in open(file_path))
    assert num_lines == tot_steps / save_step_log + 1

    ener_df = pd.read_csv(file_path)
    assert ener_df.shape[0] == tot_steps / save_step_log
    assert ener_df["Step"].iloc[-1] == tot_steps
    assert ener_df["Aim Temp (K)"].max() <= 500.0
    assert ener_df["Aim Temp (K)"].min() == 300.0
    frac_cols = [col for col in ener_df.columns if col.startswith("E frac ")]
    assert sum(ener_df[col].mean() for col in frac_cols) == pytest.approx(-709, abs=500)
    assert ener_df["E solute not scaled (kJ/mole)"].mean() == pytest.approx(
        448.0, abs=50
    )
    assert ener_df["E solvent (kJ/mole)"].mean() == pytest.approx(-60000, abs=20000)
    assert ener_df["E solvent-solute (kJ/mole)"].mean() == pytest.approx(-2272, abs=500)

    file_path = os.path.join(tmp_path, f"{name}_sst2_final.xml")
    assert os.path.exists(file_path)
    num_lines = sum(1 for _ in open(file_path))
    assert num_lines == pytest.approx(9000, abs=2000)

    file_path = os.path.join(tmp_path, f"{name}_sst2.pdb")
    assert os.path.exists(file_path)
    num_lines = sum(1 for _ in open(file_path))
    assert num_lines == pytest.approx(4000, abs=2000)
