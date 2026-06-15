#!/usr/bin/env python3
# coding: utf-8


"""
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

__author__ = "Samuel Murail"
__version__ = "0.0.1"

import openmm.unit as unit
import random
import os
from sys import stdout
import pandas as pd
import numpy as np
import logging

# Logging
logger = logging.getLogger(__name__)

from SST2.rest2 import run_rest2


class SST2Reporter(object):
    def __init__(self, sst2):
        self.sst2 = sst2

    def describeNextReport(self, simulation):
        steps1 = (
            self.sst2.tempChangeInterval
            - simulation.currentStep % self.sst2.tempChangeInterval
        )
        steps2 = (
            self.sst2.reportInterval - simulation.currentStep % self.sst2.reportInterval
        )
        steps = min(steps1, steps2)
        isUpdateAttempt = steps1 == steps
        return (steps, False, isUpdateAttempt, False, isUpdateAttempt)
    
    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The current simulation
        state : State
            The current state of the simulation
        """
        energie_group = self.sst2.rest2.compute_all_energies()

        E_frac_dict   = energie_group[0]  # {frac: unscaled_energy}
        E_pw_unscaled = energie_group[3]  # unscaled E_pw

        temp_idx = self.sst2.currentTemperature

        # ------------------------------------------------------------------ #
        # Update running averages (online mean: avg += (x - avg) / n)        #
        # ------------------------------------------------------------------ #
        self.sst2._e_num[temp_idx] += 1
        n = self.sst2._e_num[temp_idx]

        # Per-fraction unscaled averages
        for k, frac in enumerate(self.sst2._e_fracs):
            E_t = E_frac_dict.get(frac, 0 * unit.kilojoules_per_mole)
            self.sst2._e_frac_avg[temp_idx][k] += (
                E_t - self.sst2._e_frac_avg[temp_idx][k]
            ) / n

        # Solute-solvent unscaled average
        self.sst2._e_solute_solv_avg[temp_idx] += (
            E_pw_unscaled - self.sst2._e_solute_solv_avg[temp_idx]
        ) / n

        # ------------------------------------------------------------------ #
        # Report and temperature change                                       #
        # ------------------------------------------------------------------ #
        if simulation.currentStep % self.sst2.reportInterval == 0:
            self.sst2._writeReport(energie_group)

        if simulation.currentStep % self.sst2.tempChangeInterval == 0:
            self.sst2._attemptTemperatureChange(E_frac_dict, E_pw_unscaled)


class SST2(object):
    """SimulatedTempering implements the simulated tempering algorithm for
    accelerated sampling.

    It runs a simulation while allowing the temperature to vary.  At high
    temperatures, it can more easily cross energy barriers to explore a wider
    area of conformation space.  At low temperatures, it can thoroughly
    explore each local region.  For details, see Marinari, E. and Parisi, G.,
    Europhys. Lett. 19(6). pp. 451-458 (1992).

    The set of temperatures to sample can be specified in two ways.  First,
    you can explicitly provide a list
    of temperatures by using the "temperatures" argument.  Alternatively,
    you can specify the minimum and
    maximum temperatures, and the total number of temperatures to use.
    The temperatures are chosen spaced
    exponentially between the two extremes.  For example,

    st = SimulatedTempering(simulation, numTemperatures=15,
                            minTemperature=300*kelvin,
                            maxTemperature=450*kelvin)

    After creating the SimulatedTempering object, call step() on it to
    run the simulation.

    Transitions between temperatures are performed at regular intervals,
    as specified by the "tempChangeInterval" argument.  For each transition,
    a new temperature is selected using the independence sampling method, as
    described in Chodera, J. and Shirts, M., J. Chem. Phys. 135, 194110
    (2011).

    Simulated tempering requires a "weight factor" for each temperature.
    Ideally, these should be chosen so
    the simulation spends equal time at every temperature.  You can specify
    the list of weights to use with the optional "weights" argument.  If
    this is omitted, weights are selected automatically using the Wang-Landau
    algorithm as described in Wang, F. and Landau, D. P., Phys. Rev. Lett.
    86(10), pp. 2050-2053 (2001).

    To properly analyze the results of the simulation, it is important
    to know the temperature and weight factors at every point in time.
    The SimulatedTempering object functions as a reporter, writing this
    information to a file or stdout at regular intervals (which should
    match the interval at which you save frames from the simulation).
    You can specify the output file and reporting interval with the
    "reportFile" and "reportInterval" arguments.

    Parameters
    ----------
    rest2: REST2
        The REST2 object defining the System, Context, and Integrator to use
    simulation: Simulation
        The Simulation defining the System, Context, and Integrator to use

    Methods
    -------
    step(steps)
        Run a number of time steps.
    """

    def __init__(
        self,
        rest2,
        temperatures,
        refTemperature=None,
        weights=None,
        tempChangeInterval=25,
        reportInterval=1000,
        reportFile=stdout,
        restart_files=None,
        restart_files_full=None,
    ):
        """Create a new SimulatedTempering.

        Parameters
        ----------
        simulation: Simulation
            The Simulation defining the System, Context, and Integrator to use
        temperatures: list
            The list of temperatures to use for tempering, in increasing order
        refTemperature: temperature
            The reference temperature to use for tempering. If this is not specified, the first temperature in the list is used.
        weights: list
            The weight factor for each temperature.  If none, weights are selected automatically.
        tempChangeInterval: int
            The interval (in time steps) at which to attempt transitions between temperatures
        reportInterval: int
            The interval (in time steps) at which to write information to the report file
        reportFile: string or file
            The file to write reporting information to, specified as a file name or file object
        restart_files: list of strings
            Files to read restart information to, specified as a file name
        restart_files_full: string
            Full Rest2 files to read restart information to, specified as a file name
        """
        self.rest2 = rest2
        self.simulation = rest2.simulation

        numTemperatures = len(temperatures)
        self.temperatures = [
            t.in_units_of(unit.kelvin) if unit.is_quantity(t) else t * unit.kelvin
            for t in temperatures
        ]
        minTemperature = self.temperatures[0]
        maxTemperature = self.temperatures[-1]

        if refTemperature is None:
            self.refTemperature = minTemperature
        else:
            if unit.is_quantity(refTemperature):
                self.refTemperature = refTemperature.in_units_of(unit.kelvin)
            else:
                self.refTemperature = refTemperature * unit.kelvin

        assert (
            self.refTemperature in self.temperatures
        ), f"Reference temperature {self.refTemperature} not in temperatures_list {self.temperatures}"
        self.temp_ref_index = self.temperatures.index(self.refTemperature)

        if any(
            self.temperatures[i] >= self.temperatures[i + 1]
            for i in range(numTemperatures - 1)
        ):
            raise ValueError("The temperatures must be in strictly increasing order")

        logger.info(
            f"Min={minTemperature}, Ref={refTemperature}, Max={maxTemperature}, temp_list={[temp._value for temp in self.temperatures]}"
        )
        self.tempChangeInterval = tempChangeInterval
        self.reportInterval = reportInterval
        self.inverseTemperatures = [
            1.0 / (unit.MOLAR_GAS_CONSTANT_R * t) for t in self.temperatures
        ]
        self.lambdas = [
            self.temperatures[self.temp_ref_index] / T
            for T in self.temperatures
        ]
        # If necessary, open the file we will write reports to.

        self._openedFile = isinstance(reportFile, str)
        if self._openedFile:
            self._out = open(reportFile, "w", 1)
        else:
            self._out = reportFile

        # Initialize the weights.

        if weights is None:
            first_temp_index = self.compute_starting_weight(
                restart_files, restart_files_full
            )
            self._updateWeights = True
        else:
            self._weights = weights
            self._updateWeights = False

        # Select the initial temperature.
        if restart_files is None:
            self.currentTemperature = 0
        elif weights is None:
            self.currentTemperature = first_temp_index
        else:
            # Need to treat the case where weights is not None and restart_files is not None
            # TO CHANGE ! This is BAD MOKAY !!!!! :
            self.currentTemperature = 0

        # print(self.temperatures[self.currentTemperature])
        # self.simulation.integrator.setTemperature(self.temperatures[self.currentTemperature])
        self.rest2.scale_nonbonded_torsion(
            self.temperatures[self.temp_ref_index]
            / self.temperatures[self.currentTemperature]
        )
        # Add a reporter to the simulation which will handle the updates and reports.
        self.simulation.reporters.append(SST2Reporter(self))

        # Write out the header line.

        frac_headers = [
            f"E frac {frac} (kJ/mole)" for frac in self._e_fracs
        ]
        headers = (
            ["Step", "Aim Temp (K)"]
            + frac_headers
            + [
                "E solute not scaled (kJ/mole)",
                "E solvent (kJ/mole)",
                "E solvent-solute (kJ/mole)",
            ]
        )
        print(",".join(headers), file=self._out)

        
    def compute_starting_weight(self, restart_files, restart_files_full):
        """Compute the weight factor for each temperature.

        Parameters
        ----------
        restart_files : list of strings
            Files to read restart information from, specified as file names
        restart_files_full : list of strings
            Full REST2 files to read restart information from, specified as file names

        Returns
        -------
        first_temp_index : int
            Index of the last used temperature to use
        """
        numTemperatures = len(self.temperatures)

        # ------------------------------------------------------------------ #
        # Initialize energy arrays                                            #
        # ------------------------------------------------------------------ #
        self._e_num  = [0] * numTemperatures

        # Per-fraction averages: fracs from REST2 force setup
        # e.g. [0.25, 0.5, 0.75, 1.0] for torsions, [0.5] for LJ14 boundary
        self._e_fracs = self.rest2.fractional_terms
        # shape: (n_temps, n_fracs) — unscaled energy per fraction group
        self._e_frac_avg = [
            [0.0 * unit.kilojoules_per_mole] * len(self._e_fracs)
            for _ in range(numTemperatures)
        ]

        # Solute-solvent (frac=0.5) average — unscaled E_pw
        self._e_solute_solv_avg = [
            0.0 * unit.kilojoules_per_mole
        ] * numTemperatures

        # Weights
        self._weights = [0.0] * numTemperatures

        # ------------------------------------------------------------------ #
        # Restart case                                                        #
        # ------------------------------------------------------------------ #
        if restart_files is not None and restart_files_full is not None:

            # --- Load and concatenate restart CSVs ---
            df_sim  = pd.read_csv(restart_files[0])
            df_temp = pd.read_csv(restart_files_full[0])

            for i in range(1, len(restart_files)):
                logger.info(f"Reading part {i}")
                df_sim = (
                    pd.concat(
                        [df_sim, pd.read_csv(restart_files[i])],
                        axis=0, join="outer",
                    )
                    .reset_index(drop=True)
                )
                df_temp = (
                    pd.concat(
                        [df_temp, pd.read_csv(restart_files_full[i])],
                        axis=0, join="outer",
                    )
                    .reset_index(drop=True)
                )

            # Remove NaN rows (rare cases of crashes)
            df_sim  = df_sim[df_sim.iloc[:, 0].notna()].copy()
            df_temp = df_temp[df_temp.iloc[:, 0].notna()].copy()

            df_sim["Temperature (K)"] = df_temp["Aim Temp (K)"]
            temp_array = np.sort(df_sim["Temperature (K)"].unique())
            logger.info(f"Temperatures found: {temp_array}")

            # --- Per-fraction column names ---
            # Expected CSV columns: "E frac 0.25 (kJ/mole)", "E frac 0.5 (kJ/mole)", ...
            frac_col = {
                frac: f"E frac {frac} (kJ/mole)"
                for frac in self._e_fracs
            }

            # Check all fraction columns are present
            missing = [
                col for col in frac_col.values()
                if col not in df_temp.columns
            ]
            if missing:
                # Fallback: try to read legacy single "E solute scaled" column
                # and assign it entirely to frac=1.0
                logger.warning(
                    f"Per-fraction columns not found: {missing}. "
                    f"Falling back to legacy 'E solute scaled (kJ/mole)' column."
                )
                use_legacy = True
            else:
                use_legacy = False

            # --- Fill per-temperature averages ---
            for temp_index, temp in enumerate(temp_array):
                df_local = df_temp[df_temp["Aim Temp (K)"] == temp]
                self._e_num[temp_index] = len(df_local)

                if use_legacy:
                    # Legacy: "E solute scaled (kJ/mole)" stored λ·E_pp (scaled).
                    # Divide by λ_i to recover the unscaled energy expected here.
                    legacy_avg = (
                        df_local["E solute scaled (kJ/mole)"].mean()
                        / self.lambdas[temp_index]
                        * unit.kilojoules_per_mole
                    )
                    for k, frac in enumerate(self._e_fracs):
                        self._e_frac_avg[temp_index][k] = (
                            legacy_avg if frac == 1.0
                            else 0.0 * unit.kilojoules_per_mole
                        )
                else:
                    # Per-fraction averages from CSV
                    for k, frac in enumerate(self._e_fracs):
                        self._e_frac_avg[temp_index][k] = (
                            df_local[frac_col[frac]].mean()
                            * unit.kilojoules_per_mole
                        )

                # Solute-solvent average
                self._e_solute_solv_avg[temp_index] = (
                    df_local["E solvent-solute (kJ/mole)"].mean()
                    * unit.kilojoules_per_mole
                )

            # --- Find last used temperature index ---
            first_temp_index = 0
            for index, row in df_sim.iloc[::-1].iterrows():
                if index % (50 * 10) == 0:
                    first_temp_index = np.where(
                        temp_array == row["Temperature (K)"]
                    )[0][0]
                    break

            # --- Logging ---
            logger.info(f"Sample counts per temperature: {self._e_num}")
            for k, frac in enumerate(self._e_fracs):
                logger.info(
                    f"E frac {frac} averages: "
                    f"{[avg._value for avg in [row[k] for row in self._e_frac_avg]]}"
                )
            logger.info(
                f"E solute-solvent averages: "
                f"{[e._value for e in self._e_solute_solv_avg]}"
            )
            logger.info(
                f"Last temperature: {temp_array[first_temp_index]} K"
            )

            return first_temp_index

        else:
            return 0
        
    def _writeReport(self, energie_group):
        """Write out a line to the report."""

        E_frac_dict         = energie_group[0]  # {frac: unscaled_energy}
        E_solute_not_scaled = energie_group[1]
        E_solvent           = energie_group[2]
        E_solute_solvent    = energie_group[3]

        temperature = self.temperatures[self.currentTemperature].value_in_unit(
            unit.kelvin
        )

        # Per-fraction values in sorted fraction order
        frac_values = [
            E_frac_dict.get(frac, 0 * unit.kilojoules_per_mole).value_in_unit(
                unit.kilojoule_per_mole
            )
            for frac in self._e_fracs
        ]

        values = (
            [temperature]
            + frac_values
            + [
                E_solute_not_scaled.value_in_unit(unit.kilojoule_per_mole),
                E_solvent.value_in_unit(unit.kilojoule_per_mole),
                E_solute_solvent.value_in_unit(unit.kilojoule_per_mole),
            ]
        )

        print(
            ("%d," % self.simulation.currentStep)
            + ",".join("%g" % v for v in values),
            file=self._out,
        )
    
    def __del__(self):
        if self._openedFile:
            self._out.close()

    @property
    def weights(self):
        return [x - self._weights[0] for x in self._weights]

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self.simulation.step(steps)

    def _compute_weight(self, i, j):
        r"""Compute the difference of weight $w_j - w_i$:

        $$(w_j - w_i) = \sum_t (\lambda_j^{f_t} - \lambda_i^{f_t})
        \frac{\langle \tilde{E}_t \rangle_i + \langle \tilde{E}_t \rangle_j}{2}
        + (\lambda_j^{1/2} - \lambda_i^{1/2})
        \frac{\langle \tilde{E}_{pw} \rangle_i + \langle \tilde{E}_{pw} \rangle_j}{2}$$

        If state j has no samples yet, only state i averages are used.

        Parameters
        ----------
        i : int
            Current temperature index
        j : int
            Target temperature index

        Returns
        -------
        weight : unit.Quantity
            Weight difference w_j - w_i
        """
        lambda_i = self.lambdas[i]
        lambda_j = self.lambdas[j]
        j_has_samples = self._e_num[j] != 0

        weight = 0 * unit.kilojoules_per_mole

        # Sum over all fractional energy terms
        for k, frac in enumerate(self._e_fracs):
            avg_i = self._e_frac_avg[i][k]
            avg_j = self._e_frac_avg[j][k] if j_has_samples else avg_i
            avg   = (avg_i + avg_j) / 2 if j_has_samples else avg_i

            weight += (lambda_j ** frac - lambda_i ** frac) * avg

        # Solute-solvent term (frac = 0.5)
        avg_pw_i = self._e_solute_solv_avg[i]
        avg_pw_j = self._e_solute_solv_avg[j] if j_has_samples else avg_pw_i
        avg_pw   = (avg_pw_i + avg_pw_j) / 2 if j_has_samples else avg_pw_i

        weight += (lambda_j ** 0.5 - lambda_i ** 0.5) * avg_pw

        return weight



    def _attemptTemperatureChange(self, E_frac_dict, E_pw_unscaled):
        """Attempt to move to a different temperature.

        Acceptance log probability:
        $$\\Delta_{i \\to j} = \\sum_t (\\lambda_i^{f_t} - \\lambda_j^{f_t})
        \\cdot \\tilde{E}_t(X)
        + (\\lambda_i^{1/2} - \\lambda_j^{1/2}) \\cdot \\tilde{E}_{pw}(X)
        + (w_j - w_i)$$

        Parameters
        ----------
        E_frac_dict : dict {frac: unit.Quantity}
            Unscaled energy per fraction group at current configuration
        E_pw_unscaled : unit.Quantity
            Unscaled solute-solvent energy at current configuration
        """
        temp_i   = self.currentTemperature
        lambda_i = self.lambdas[temp_i]

        # Build list of neighboring temperature indices
        temp_list = []
        if temp_i != 0:
            temp_list.append(temp_i - 1)
        if temp_i < len(self._weights) - 1:
            temp_list.append(temp_i + 1)

        logProbability = []
        for j in temp_list:
            lambda_j = self.lambdas[j]

            # Sum over fractional terms: (lambda_i^f - lambda_j^f) * E_t_unscaled
            log_prob = 0 * unit.kilojoules_per_mole
            for frac in self._e_fracs:
                E_t = E_frac_dict.get(frac, 0 * unit.kilojoules_per_mole)
                log_prob += (lambda_i ** frac - lambda_j ** frac) * E_t

            # Solute-solvent contribution: (lambda_i^0.5 - lambda_j^0.5) * E_pw
            log_prob += (lambda_i ** 0.5 - lambda_j ** 0.5) * E_pw_unscaled

            # Weight difference w_j - w_i
            log_prob += self._compute_weight(temp_i, j)

            logProbability.append(log_prob._value)

        # Metropolis criterion: p = exp(log_prob / kT)
        # Note: log_prob is already in energy units (kJ/mol),
        # kT = kB * T_ref in kJ/mol
        kT = unit.BOLTZMANN_CONSTANT_kB * self.temperatures[self.temp_ref_index] * unit.AVOGADRO_CONSTANT_NA
        kT = kT.value_in_unit(unit.kilojoule_per_mole)

        probability = [np.exp(lp / kT) for lp in logProbability]

        # Shuffle to avoid systematic bias toward lower temperature
        index_list = list(range(len(temp_list)))
        random.shuffle(index_list)
        r = random.random()

        for i in index_list:
            if r < probability[i]:
                self.currentTemperature = temp_list[i]
                self.rest2.scale_nonbonded_torsion(
                    self.lambdas[temp_list[i]]
                )
                break

def run_sst2(
    sys_rest2,
    generic_name,
    tot_steps,
    dt,
    temperatures,
    ref_temp,
    save_step_dcd=100000,
    save_step_log=500,
    save_step_rest2=500,
    tempChangeInterval=500,
    reportInterval=500,
    overwrite=False,
    save_checkpoint_steps=None,
):
    """
    Run a SST2 simulation.

    Parameters
    ----------
    sys_rest2 : Rest2 object
        The system to simulate.
    generic_name : str
        Generic name for the output files.
    tot_steps : int
        Total number of steps to run.
    dt : float
        Time step in fs.
    temperatures : list of float
        List of temperatures to simulate.
    ref_temp : float
        Reference temperature.
    save_step_dcd : int, optional
        Number of steps between each DCD save. The default is 100000.
    save_step_log : int, optional
        Number of steps between each log save. The default is 500.
    save_step_rest2 : int, optional
        Number of steps between each Rest2 save. The default is 500.
    tempChangeInterval : int, optional
        Number of steps between each temperature change. The default is 500.
    reportInterval : int, optional
        Number of steps between each report. The default is 500.
    overwrite : bool, optional
        Overwrite the previous simulation. The default is True.
    save_checkpoint_steps : int, optional
        Number of steps between each checkpoint save. The default is None.


    """

    if unit.is_quantity(ref_temp):
        ref_temp = ref_temp.in_units_of(unit.kelvin)
    else:
        ref_temp *= unit.kelvin

    assert (
        ref_temp in temperatures
    ), f"Reference temperature {ref_temp} not in temperatures_list {temperatures}"

    report_sst2 = f"{generic_name}_sst2_full.csv"
    restart_files = None
    restart_files_full = None
    tot_steps = np.ceil(tot_steps)

    if not overwrite and os.path.isfile(report_sst2):
        logger.info(
            f"File {generic_name}_sst2_full.csv exists already, restart run_sst2() step"
        )
        # Get part number
        part = 2

        report_sst2 = f"{generic_name}_sst2_full_part_{part}.csv"
        report_simple_sst2 = f"{generic_name}_sst2_part_{part}.csv"

        restart_files = [f"{generic_name}_sst2.csv"]
        restart_files_full = [f"{generic_name}_sst2_full.csv"]

        while os.path.isfile(report_sst2):
            restart_files.append(report_simple_sst2)
            restart_files_full.append(report_sst2)
            report_sst2 = f"{generic_name}_sst2_full_part_{part}.csv"
            report_simple_sst2 = f"{generic_name}_sst2_part_{part}.csv"
            part += 1

        if part != 2:
            restart_files = restart_files[:-1]
            restart_files_full = restart_files_full[:-1]

        logger.info(f"Using restart file : {restart_files}")

    sys_rest2.simulation.reporters = []
    sys_rest2.simulation.currentStep = 0

    sst2 = SST2(
        sys_rest2,
        temperatures=temperatures,
        refTemperature=ref_temp,
        tempChangeInterval=tempChangeInterval,
        reportFile=report_sst2,
        reportInterval=reportInterval,
        restart_files=restart_files,
        restart_files_full=restart_files_full,
    )

    logger.info(f"- Launch SST2")
    run_rest2(
        sst2.rest2,
        f"{generic_name}_sst2",
        tot_steps=tot_steps,
        dt=dt,
        save_step_dcd=save_step_dcd,
        save_step_log=save_step_log,
        save_step_rest2=reportInterval,
        add_REST2_reporter=False,
        remove_reporters=False,
        save_checkpoint_steps=save_checkpoint_steps,
    )
