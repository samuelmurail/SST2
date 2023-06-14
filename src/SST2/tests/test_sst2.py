#!/usr/bin/env python3
# coding: utf-8

"""
Tests for rest2 functions
"""

import os
import pytest

import openmm
from openmm import unit
import openmm.app as app

import pdb_numpy

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
        PDB_5AWL, os.path.join(tmp_path, f"{name}_fixed.pdb"), pH=7.0, overwrite=False
    )

    forcefield_files = ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    forcefield = app.ForceField(*forcefield_files)

    tools.create_water_box(
        os.path.join(tmp_path, f"{name}_fixed.pdb"),
        os.path.join(tmp_path, f"{name}_water.pdb"),
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

    pdb = app.PDBFile(os.path.join(tmp_path, f"{name}_water.pdb"))

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
    )

    tools.minimize(
        sys_rest2.simulation,
        os.path.join(tmp_path, f"{name}_em_water.pdb"),
        pdb.topology,
        maxIterations=1000,
        overwrite=False,
    )

    sys_rest2.simulation.context.setVelocitiesToTemperature(temperature)

    rest2.run_rest2(
        sys_rest2,
        os.path.join(tmp_path, f"{name}_equi_water"),
        dt=dt,
        tot_steps=nsteps,
        save_step_dcd=100000,
        save_step_log=10000,
        save_step_rest2=500,
    )

    ladder_num = tools.compute_ladder_num(
        os.path.join(tmp_path, f"{name}_equi_water"), 300, 500
    )

    assert ladder_num == 3

    temp_list = sst2.compute_temperature_list(
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
    save_check_steps = 1000

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