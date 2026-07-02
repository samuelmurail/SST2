#!/usr/bin/env python3
# coding: utf-8

import os
import math
import logging
import copy

import pandas as pd
import numpy as np
import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt
import openmm.unit as unit
from scipy.ndimage import gaussian_filter1d

# Logging
logger = logging.getLogger(__name__)

def read_SST2_data(
    generic_name, dt=0.004, full_sep=",", save_step_dcd=100000, lambda_T_ref=300.0
):
    """Read the SST2 data from the csv files.

    Supports both the new per-fraction format:
        Step, Aim Temp (K), E frac 0.25 (kJ/mole), ..., E solute not scaled (kJ/mole), ...
    and the legacy single-column format:
        Step, Aim Temp (K), E solute scaled (kJ/mole), ...

    Parameters
    ----------
    generic_name : str
        Generic name of the csv files (without the `.csv` extension).
    dt : float, optional
        Time step in ps. Default is 0.004 ps.
    full_sep : str, optional
        Separator used in the full csv file. Default is ",".
    save_step_dcd : int, optional
        Step number used in the dcd file. Default is 100000.
    lambda_T_ref : float, optional
        Reference temperature for lambda. Default is 300.0.

    Returns
    -------
    df_all : pandas.DataFrame
        Dataframe with all the data.
    """
    return read_ST_data(
        generic_name=generic_name,
        dt=dt,
        full_sep=full_sep,
        save_step_dcd=save_step_dcd,
        lambda_T_ref=lambda_T_ref,
    )


def read_ST_data(
    generic_name,
    dt=0.004,
    fields=None,
    full_sep=",",
    save_step_dcd=100000,
    lambda_T_ref=None,
):
    """Read SST2/ST data from csv files, merging restart parts if present.

    Supports both the new per-fraction format and the legacy single-column
    format. If `fields` is None, all columns are read.

    Parameters
    ----------
    generic_name : str
        Generic name of the csv files (without the `.csv` extension).
    dt : float, optional
        Time step in ps. Default is 0.004 ps.
    fields : list or None, optional
        List of columns to read. If None, all columns are read.
        Default is None.
    full_sep : str, optional
        Separator used in the full csv file. Default is ",".
    save_step_dcd : int, optional
        Step number used in the dcd file. Default is 100000.
    lambda_T_ref : float or None, optional
        Reference temperature for lambda. Default is None.

    Returns
    -------
    df_all : pandas.DataFrame
        Dataframe with all the data.
    """

    def read_full(path):
        """Read full csv, selecting fields if specified."""
        if fields is not None:
            # Only load columns that actually exist in the file
            available = pd.read_csv(path, sep=full_sep, nrows=0).columns.tolist()
            cols = [f for f in fields if f in available]
            return pd.read_csv(path, sep=full_sep, usecols=cols)
        return pd.read_csv(path, sep=full_sep)

    # ------------------------------------------------------------------ #
    # Discover restart parts                                               #
    # ------------------------------------------------------------------ #
    part = 1
    while os.path.isfile(f"{generic_name}_part_{part + 1}.csv"):
        part += 1
    logger.info(f"Found {part} csv part(s) to read.")

    # ------------------------------------------------------------------ #
    # Load part 1                                                          #
    # ------------------------------------------------------------------ #
    logger.info("Reading part 1")
    df_temp_list = [read_full(f"{generic_name}_full.csv")]
    df_sim_list  = [pd.read_csv(f"{generic_name}.csv")]

    # ------------------------------------------------------------------ #
    # Load remaining parts                                                 #
    # ------------------------------------------------------------------ #
    for i in range(2, part + 1):
        logger.info(f"Reading part {i}")

        df_temp_part = read_full(f"{generic_name}_full_part_{i}.csv")
        df_sim_part  = pd.read_csv(f"{generic_name}_part_{i}.csv")

        # Detect step-reset restarts (DCD format limitation)
        last_old_step  = df_temp_list[-1].iloc[-1, 0]
        first_new_step = df_temp_part.iloc[0, 0]

        if first_new_step < last_old_step - save_step_dcd:
            step_col = df_temp_part.columns[0]
            sim_step_col = (
                '#"Step"' if '#"Step"' in df_sim_part.columns else "Step"
            )
            chk_step = (
                df_sim_list[-1][sim_step_col][
                    df_sim_list[-1][sim_step_col] % save_step_dcd == 0
                ].iloc[-1]
                - first_new_step
            )
            df_temp_part[step_col]    += chk_step
            df_sim_part[sim_step_col] += chk_step
            logger.info(
                f"Step offset of {chk_step} applied to part {i}"
            )

        df_temp_list.append(df_temp_part)
        df_sim_list.append(df_sim_part)

    # ------------------------------------------------------------------ #
    # Concatenate all parts                                                #
    # ------------------------------------------------------------------ #
    df_temp = pd.concat(df_temp_list, axis=0, join="outer", ignore_index=True)
    df_sim  = pd.concat(df_sim_list,  axis=0, join="outer", ignore_index=True)
    del df_temp_list, df_sim_list

    logger.info(f"df_sim  length: {len(df_sim)}")
    logger.info(f"df_temp length: {len(df_temp)}")

    # ------------------------------------------------------------------ #
    # Normalize column names (legacy ST format compatibility)              #
    # ------------------------------------------------------------------ #
    df_temp = df_temp.rename(columns={
        '#"Steps"': "Step",
        "Temperature (K)": "Aim Temp (K)",
    })
    df_sim = df_sim.rename(columns={'#"Step"': "Step"})

    # ------------------------------------------------------------------ #
    # Legacy format: rename single scaled column to frac 1.0              #
    # ------------------------------------------------------------------ #
    if "E solute scaled (kJ/mole)" in df_temp.columns:
        logger.info(
            "Legacy format detected: renaming "
            "'E solute scaled (kJ/mole)' → 'E frac 1.0 (kJ/mole)'"
        )
        df_temp = df_temp.rename(
            columns={"E solute scaled (kJ/mole)": "E frac 1.0 (kJ/mole)"}
        )

    # ------------------------------------------------------------------ #
    # Align lengths                                                        #
    # ------------------------------------------------------------------ #
    max_step = min(len(df_sim), len(df_temp))
    logger.info(f"Using length: {max_step}")
    df_sim  = df_sim.iloc[:max_step]
    df_temp = df_temp.iloc[:max_step]

    # Drop Step from df_temp before merge (it comes from df_sim)
    if "Step" in df_temp.columns:
        df_temp = df_temp.drop(columns=["Step"])

    # ------------------------------------------------------------------ #
    # Merge sim + temp dataframes                                          #
    # ------------------------------------------------------------------ #
    df_all = pd.concat([df_temp, df_sim], axis=1)
    del df_temp, df_sim

    # ------------------------------------------------------------------ #
    # Add derived columns                                                  #
    # ------------------------------------------------------------------ #
    df_all[r"$Time\;(\mu s)$"] = df_all["Step"] * dt / 1e6

    df_all["Temp (K)"] = pd.Categorical(df_all["Aim Temp (K)"])

    if lambda_T_ref is not None:
        df_all["lambda"] = lambda_T_ref / df_all["Aim Temp (K)"]
        df_all[r"$\lambda$"] = pd.Categorical(
            df_all["lambda"].round(2)
        )
        df_all[r"$\lambda$"] = df_all[r"$\lambda$"].cat.reorder_categories(
            df_all[r"$\lambda$"].cat.categories[::-1]
        )

    # Drop NaN rows from rare crashes
    if df_all["Step"].isna().any():
        logger.info("Removing NaN rows")
        df_all = df_all.dropna()

    return df_all

def compute_exchange_prob(
    df, temp_col="Aim Temp (K)", time_ax_name=r"$Time\;(\mu s)$", exchange_time=2
):
    """
    Compute the exchange probability and the round trip time
    for a given dataframe.
    The dataframe should be the result of a SST2 simulation.
    The dataframe should have a column with the temperature
    and a column with the time.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the SST2 simulation data.
    temp_col : str, optional
        Name of the column with the temperature. The default is "Aim Temp (K)".
    time_ax_name : str, optional
        Name of the column with the time. The default is r"$Time\;(\mu s)$".
    exchange_time : float, optional
        Time in ps between two exchange. The default is 2 ps.

    Returns
    -------
    ex_prob : float
        Exchange probability.
    trip_time : float
        Round trip time in ns.
    """

    temp_list = df[temp_col].unique()
    min_temp = temp_list[0]
    max_temp = temp_list[-1]
    logger.info(
        f"Min_temp = {min_temp:.2f}, Max temp. = {max_temp:.2f}, #Rungs = {len(temp_list)}"
    )

    # Compute exchange probability:
    last_temp = min_temp
    temp_change_num = 0

    # Compute Round trip time
    # or time to go from min to
    # max and back to min temp:

    time_list = []
    target_temp = last_temp
    step_num = 0

    step_time = df.loc[1, time_ax_name] - df.loc[0, time_ax_name]
    step_time *= 1e6  # ps
    logger.info(f"Step time = {step_time:.2f} ps")

    trip_flag = False

    all_temp_change_num = {temp: 0 for temp in temp_list}
    all_temp_num = {temp: 0 for temp in temp_list}
    temp_change_index = []
    change_index = 0
    sign = 1

    for temp in df[temp_col]:
        all_temp_num[temp] += 1
        # round trip time
        if temp == min_temp:
            if trip_flag:
                time_list.append(step_num * step_time)
            step_num = 0
            trip_flag = False
        elif temp == max_temp:
            step_num += 1
            trip_flag = True
        else:
            step_num += 1
        # Exchange prob
        if temp != last_temp:
            all_temp_change_num[temp] += 1
            temp_change_num += 1
            sign = temp - last_temp
            last_temp = temp
            change_index = 1
        else:
            change_index += 1
        temp_change_index.append(math.copysign(change_index, sign))

    # print(temp_change_index)
    df["Temp Change index"] = temp_change_index

    for temp in temp_list:
        all_temp_change_num[temp] /= all_temp_num[temp]
        all_temp_change_num[temp] *= exchange_time / step_time

    keys = list(all_temp_change_num.keys())
    # get values in the same order as keys, and parse percentage values
    vals = [all_temp_change_num[k] for k in keys]

    ax = sns.barplot(x=temp_list, y=vals, hue=temp_list)
    plt.xlabel(temp_col)
    plt.ylabel(r"$p()$")
    # plt.title(r"Transition probability at each rung")
    ax.get_legend().remove()

    # print(all_temp_change_num)
    logger.info(f'exchange step = {exchange_time / step_time}')

    ex_prob = temp_change_num / len(df) * exchange_time / step_time

    if len(time_list) != 0:
        trip_time = sum(time_list) / len(time_list) / 1000  # ns
    else:
        trip_time = None

    return ex_prob, trip_time


def plot_lineplot_avg(
    df, x, y, quant=None, color="black", max_data=50000, avg_win=1000, alpha=0.3
):
    """
    Plot a lineplot with a gaussian filter on the y axis.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the data to plot.
    x : str
        Name of the column with the x axis data.
    y : str
        Name of the column with the y axis data.
    quant : float, optional
        Quantile to use to filter the data. The default is None.
    color : str, optional
        Color of the line. The default is "black".
    max_data : int, optional
        Maximum number of data point to plot. The default is 50000.
    avg_win : int, optional
        Window size of the gaussian filter. The default is 1000.

    Returns
    -------
    g : matplotlib.axes._subplots.AxesSubplot
        Axes of the plot.
    """

    local_df = filter_df(df, max_data)

    local_df["avg"] = gaussian_filter1d(local_df[y], avg_win)

    g = sns.lineplot(data=local_df, x=x, y=y, lw=0.1, color=color, alpha=alpha)
    sns.lineplot(data=local_df, x=x, y="avg", lw=2, color=color)
    if quant is not None:
        x_min, x_max = get_quant_min_max(local_df[x], quant)
        y_min, y_max = get_quant_min_max(local_df[y], quant)
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

    return g


def filter_df(df, max_point_number):
    """
    Filter a dataframe to keep a maximum number of data point.
    The dataframe is filtered with a step size computed
    to keep the maximum number of data point.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to filter.
    max_point_number : int
        Maximum number of data point to keep.

    Returns
    -------
    local_df : pandas.DataFrame
        Filtered dataframe.
    """

    steps_num = len(df)

    md_step = 1

    while steps_num > max_point_number:
        md_step += 1
        steps_num = len(df) // md_step

    local_df = df.loc[::md_step]
    local_df = local_df.reset_index(drop=True)

    return local_df


def get_quant_min_max(pd_serie, quant=0.001):
    """
    Get the min and max value of a pandas serie
    from a quantile.

    Parameters
    ----------
    pd_serie : pandas.Series
        Pandas serie to analyze.
    quant : float, optional
        Quantile to use. The default is 0.001.

    Returns
    -------
    val_min : float
        Minimum value.
    val_max : float
        Maximum value.
    """

    val_min = pd_serie.quantile(quant)
    val_max = pd_serie.quantile(1 - quant)

    return (val_min, val_max)


def plot_distri_norm(
    df,
    x,
    hue,
    x_label=None,
    max_data=50000,
    bins=100,
    element="step",
    quant=None,
    bw_adjust=None,
):
    """
    Plot a distribution plot with a gaussian filter on the y axis.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the data to plot.
    x : str
        Name of the column with the x axis data.
    hue : str
        Name of the column with the hue data.
    x_label : str, optional
        Label of the x axis. The default is None.
    max_data : int, optional
        Maximum number of data point to plot. The default is 20000.
    bins : int, optional
        Number of bins. The default is 100.
    element : str, optional
        Element of the plot. The default is "step".
    quant : float, optional
        Quantile to use to filter the data. The default is None.
    bw_adjust : float, optional
        Bandwidth adjustment for the kernel density estimate. The default is None.

    Returns
    -------
    ax1 : matplotlib.axes._subplots.AxesSubplot
        Axes of the plot.
    """

    local_df = filter_df(df, max_data)

    if x_label is None:
        x_label = x

    if bw_adjust is None:
        kde_kws = {}
    else:
        kde_kws = {"bw_adjust": bw_adjust}

    fig, ax1 = plt.subplots()

    g = sns.histplot(
        local_df,
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=x,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        hue=hue,
        element=element,
        ax=ax1,
        kde_kws=kde_kws,
    )

    if quant is not None:
        x_min, x_max = get_quant_min_max(local_df[x], quant)
        plt.xlim(x_min, x_max)
    # g.axes.flat[0].xaxis.set_major_formatter(ticker.EngFormatter())
    # plt.legend(bbox_to_anchor=(1.01, 1.0))

    legend = ax1.get_legend()
    # handles = legend.legendHandles # Deprecated
    handles = legend.legend_handles

    # Remove alpha in legend
    for h in handles:
        h.set_alpha(1.0)

    if len(handles) > 15:
        ncol = 2
        x_gap = 1.3
    else:
        ncol = 1
        x_gap = 1.2

    sns.move_legend(
        g, "lower center", bbox_to_anchor=(x_gap, 0.0), ncol=ncol, title_fontsize=14
    )

    plt.xlabel(x_label)
    return ax1


def plot_scatter(
    df,
    x,
    y,
    hue=None,
    x_label=None,
    y_label=None,
    quant=None,
    s=10,
    color=None,
    linewidth=0,
    label=None,
    legend="auto",
    alpha=None,
    max_data=50000,
):
    """
    Plot a scatter plot.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the data to plot.
    x : str
        Name of the column with the x axis data.
    y : str
        Name of the column with the y axis data.
    hue : str, optional
        Name of the column with the hue data. The default is None.
    x_label : str, optional
        Label of the x axis. The default is None.
    y_label : str, optional
        Label of the y axis. The default is None.
    quant : float, optional
        Quantile to use to filter the data. The default is None.
    s : int, optional
        Size of the points. The default is 10.
    color : str, optional
        Color of the points. The default is None.
    linewidth : float, optional
        Width of the points. The default is 0.
    label : str, optional
        Label of the plot. The default is None.
    legend : str, optional
        Position of the legend. The default is "auto".
    alpha : float, optional
        Transparency of the points. The default is None.
    max_data : int, optional
        Maximum number of data point to plot. The default is 50000.

    Returns
    -------
    g : matplotlib.axes._subplots.AxesSubplot
        Axes of the plot.
    """

    local_df = filter_df(df, max_data)

    g = sns.scatterplot(
        data=local_df,
        x=x,
        y=y,
        s=s,
        linewidth=linewidth,
        legend=legend,
        color=color,
        alpha=alpha,
        hue=hue,
    )

    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)

    if quant is not None:
        x_min, x_max = get_quant_min_max(local_df[x], quant)
        y_min, y_max = get_quant_min_max(local_df[y], quant)
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

    plt.legend(title=hue, bbox_to_anchor=(1.01, 1.0))

    return g


def plot_weight_RMSD(
    df,
    x=r"$Time\;(\mu s)$",
    hue="Temp (K)",
    ener="new_pot",
    time_ax_name=r"$Time\;(\mu s)$",
    final_weight_dict=None,
    max_data=50000,
    plot_weights=False,
):

    local_df = filter_df(df, max_data)
    temp_list = local_df["Aim Temp (K)"].unique()
    # gaussian_filter1d(local_df[y], avg_win)
    local_df = compute_moving_average(local_df, ener=ener)

    local_df = compute_weight_RMSD(
        local_df, final_weight_dict=final_weight_dict, ener=ener
    )

    if plot_weights:

        ax1 = sns.lineplot(data=local_df, x=x, y="avg_ener", hue=hue, lw=2)
        plt.ylabel(r"E $(KJ.mol^{-1})$")

        for temp in temp_list:
            avg_temp = local_df[local_df["Aim Temp (K)"] == temp][ener].mean()
            plt.axhline(avg_temp, lw=1, c="gray", linestyle=":")
        plt.legend(title=hue, bbox_to_anchor=(1.01, 1.0))
        plt.show()

    ax2 = sns.lineplot(data=local_df, x=time_ax_name, y="Weight RMSD", lw=2)
    plt.ylabel(r"RMSD $(KJ.mol^{-1})$")
    plt.title(r"Weights $f_i$ RMSD")

    return


def compute_moving_average(df, ener="new_pot", col_name="avg_ener"):

    temp_list = list(df["Aim Temp (K)"].unique())
    temp_list.sort()
    temp_list_avg = {temp: 0 for temp in temp_list}
    temp_list_num = {temp: 0 for temp in temp_list}
    mov_avg = []

    for temp, new_pot in zip(df["Aim Temp (K)"], df[ener]):
        temp_list_num[temp] += 1
        temp_list_avg[temp] += (new_pot - temp_list_avg[temp]) / temp_list_num[temp]
        mov_avg.append(temp_list_avg[temp])

    df.loc[:, col_name] = mov_avg

    return df


def compute_weight_RMSD(df, final_weight_dict=None, ener="new_pot"):

    df.loc[:, "Weight RMSD"] = 0

    if final_weight_dict is None:
        temp_final_avg = {}
    else:
        temp_final_avg = final_weight_dict

    temp_list = list(df["Aim Temp (K)"].unique())
    temp_list.sort()

    for temp in temp_list:

        if final_weight_dict is None:
            tmp_df = df[df["Aim Temp (K)"] == temp]
            temp_final_avg[temp] = tmp_df[ener].mean()

        last_avg_ener = np.nan
        weight_list = []
        for for_temp, avg_ener in zip(df["Aim Temp (K)"], df["avg_ener"]):
            if for_temp == temp:
                last_avg_ener = avg_ener
            weight_list.append(last_avg_ener)

        df.loc[:, f"weight {temp}"] = weight_list
        df["Weight RMSD"] += (df[f"weight {temp}"] - temp_final_avg[temp]) ** 2

    df["Weight RMSD"] = (df["Weight RMSD"] / len(temp_list)) ** 0.5

    return df


def plot_weight_RMSD(
    df,
    x=r"$Time\;(\mu s)$",
    hue="Temp (K)",
    ener="new_pot",
    time_ax_name=r"$Time\;(\mu s)$",
    final_weight_dict=None,
    max_data=50000,
    plot_weights=False,
):

    local_df = filter_df(df, max_data)
    temp_list = local_df["Aim Temp (K)"].unique()
    local_df = compute_moving_average(local_df, ener=ener)

    local_df = compute_weight_RMSD(
        local_df, final_weight_dict=final_weight_dict, ener=ener
    )

    if plot_weights:

        ax1 = sns.lineplot(data=local_df, x=x, y="avg_ener", hue=hue, lw=2)
        plt.ylabel(r"E $(KJ.mol^{-1})$")

        for temp in temp_list:
            avg_temp = local_df[local_df["Aim Temp (K)"] == temp][ener].mean()
            plt.axhline(avg_temp, lw=1, c="gray", linestyle=":")
        plt.legend(title=hue, bbox_to_anchor=(1.01, 1.0))
        plt.show()

    ax2 = sns.lineplot(data=local_df, x=time_ax_name, y="Weight RMSD", lw=2)
    plt.ylabel(r"RMSD $(KJ.mol^{-1})$")
    plt.title(r"Weights $f_i$ RMSD")

    return


def plot_free_energy(
    xall,
    yall,
    weights=None,
    ax=None,
    nbins=100,
    ncontours=100,
    avoid_zero_count=False,
    minener_zero=True,
    kT=2.479,
    vmin=None,
    vmax=None,
    cmap="nipy_spectral",
    cbar=True,
    cbar_label="free energy (kJ/mol)",
    cax=None,
    levels=None,
    cbar_orientation="vertical",
    norm=None,
    range=None,
    level_gap=None,
):
    """
    Plot the free energy of a 2D histogram.

    Adapted from;
    https://github.com/markovmodel/PyEMMA/blob/devel/pyemma/plots/plots2d.py

    Parameters
    ----------
    xall : np.array
        Array with the x data.
    yall : np.array
        Array with the y data.
    weights : np.array, optional
        Array with the weights. The default is None.
    ax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes of the plot. The default is None.
    nbins : int, optional
        Number of bins. The default is 100.
    ncontours : int, optional
        Number of contours. The default is 100.
    avoid_zero_count : bool, optional
        Avoid zero count. The default is False.
    minener_zero : bool, optional
        Minimum energy to zero. The default is True.
    kT : float, optional
        kT value. The default is 2.479.
    vmin : float, optional
        Minimum value. The default is None.
    vmax : float, optional
        Maximum value. The default is None.
    cmap : str, optional
        Colormap. The default is 'nipy_spectral'.
    cbar : bool, optional
        Add colorbar. The default is True.
    cbar_label : str, optional
        Label of the colorbar. The default is 'free energy (kJ/mol)'.
    cax : matplotlib.axes._subplots.AxesSubplot, optional
        Axes of the colorbar. The default is None.
    levels : int, optional
        Number of levels. The default is None.
    cbar_orientation : str, optional
        Orientation of the colorbar. The default is 'vertical'.
    norm : matplotlib.colors.Normalize, optional
        Normalize object. The default is None.
    range : list, optional
        Range of the data. The default is None.
    level_gap : float, optional
        Gap between levels. The default is None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure of the plot.
    ax : matplotlib.axes._subplots.AxesSubplot
        Axes of the plot.
    misc : dict
        Dictionary with the colorbar.

    """
    z, xedge, yedge = np.histogram2d(
        xall, yall, bins=nbins, weights=weights, range=range
    )
    if avoid_zero_count:
        z = np.maximum(z, np.min(z[z.nonzero()]))
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])

    pi = z.T / float(z.sum())
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = pi.nonzero()
    zero = np.nonzero(pi == 0)
    free_energy[nonzero] = -np.log(pi[nonzero])
    # if minener_zero:
    free_energy[nonzero] -= np.min(free_energy[nonzero])
    free_energy *= kT

    # to show the highest free energy zones in the map,
    # replace infinity by a value slightly above the maximum free energy:

    if vmax is None:
        vmax = np.max(free_energy[nonzero]) + 0.5
    # free_energy[zero] = np.max(free_energy[nonzero])+1.0
    free_energy[zero] = vmax + 0.5
    if vmin is None:
        vmin = 0
    # and fix the levels for the colormap:
    if levels is None:
        cbar_ticks = np.linspace(vmin, vmax, ncontours)
    else:
        cbar_ticks = np.linspace(vmin, vmax, levels + 1)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    if levels is None and level_gap is not None:
        max_free = np.max(free_energy[nonzero])
        levels = int(max_free // level_gap)
        logger.info(levels)

    mappable = ax.contourf(
        x,
        y,
        free_energy,
        ncontours,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        levels=cbar_ticks,
    )

    misc = dict(mappable=mappable)
    if cbar:
        if levels is not None:
            # print("ticks", levels, cbar_ticks)
            cbar = fig.colorbar(
                mappable,
                ax=ax,
                orientation=cbar_orientation,
                ticks=cbar_ticks,
                # extend='both'
            )
        else:
            cbar = fig.colorbar(mappable, ax=ax, orientation=cbar_orientation)
        cbar.set_label(cbar_label)
        misc.update(cbar=cbar)

    # if cbar:
    #    cbar = fig.colorbar(
    #        mappable,
    #        ax=ax,
    #        orientation=cbar_orientation)
    #    cbar.set_label(cbar_label)
    #    misc.update(cbar=cbar)

    return fig, ax, misc


def compute_cluster_hdbscan(pca_df, min_cluster_size=50, min_samples=50):
    """
    Cluster the data using the HDBSCAN algorithm.

    Parameters
    ----------
    pca_df : pandas.DataFrame
        Dataframe with the data to cluster.
    min_cluster_size : int, optional
        Minimum cluster size. The default is 50.
    min_samples : int, optional
        Minimum number of samples. The default is 50.

    Returns
    -------
    clust_serie : pandas.Categorical
        Categorical serie with the cluster.
    """
    import hdbscan

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size, min_samples=min_samples
    ).fit(pca_df)
    labels = clusterer.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    logger.info(
        "Number of cluster : {}, perc of non clustered points : {:.1f}%".format(
            n_clusters_, 100 * n_noise_ / len(pca_df)
        )
    )

    # count each cluster
    clust_dict = {}

    for i in range(-1, n_clusters_):
        if sum(labels == i) > 0:
            clust_dict[i] = sum(labels == i)

    # sort cluster as function of clust pop

    sorted_dict = {
        k: v
        for k, v in sorted(clust_dict.items(), key=lambda item: item[1], reverse=True)
    }

    # Create new cluster list
    new_label = np.copy(labels)

    select_list = []
    new_value_list = []
    cat_to_remove = []
    clust_new = 1

    for i, clust in enumerate(sorted_dict):

        if clust != -1:
            logger.info(
                f"Cluster:{clust_new:3}   {sum(new_label == clust):5} | {sum(new_label == clust)/len(labels):.3f}"
            )
            new_value_list.append(clust_new)
            clust_new += 1
        else:
            logger.info(
                f"Not Clustered {sum(new_label == clust):5} | {sum(new_label == clust)/len(labels):.3f}"
            )
            new_value_list.append(0)
            cat_to_remove = [0]
        select_list.append(new_label == clust)

    new_label = np.select(select_list, new_value_list, new_label)
    clust_serie = pd.Categorical(pd.Series(data=new_label)).remove_categories(
        cat_to_remove
    )

    return clust_serie


def compute_cluster_kmean(pca_df, max_cluster=20, random_state=0):
    """
    Cluster the data using the KMeans algorithm.

    Parameters
    ----------
    pca_df : pandas.DataFrame
        Dataframe with the data to cluster.
    max_cluster : int, optional
        Maximum number of cluster to test. The default is 20.
    random_state : int, optional
        Random state for the algorithm. The default is 0.

    Returns
    -------
    clust_serie : pandas.Categorical
        Categorical serie with the cluster.
    kmeans.cluster_centers_ : numpy.ndarray
        Cluster centers.

    """

    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    kmeans_kwargs = {
        "init": "random",
        "n_init": 10,
        "max_iter": 50,
        "random_state": random_state,
    }

    # A list holds the silhouette_coefficients for each k
    silhouette_coefficients = []
    kmeans_list = []

    for k in range(2, max_cluster + 1):
        logger.info(f"{k}/{max_cluster}")
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(pca_df)
        kmeans_list.append(kmeans)
        score = silhouette_score(pca_df, kmeans.labels_)
        silhouette_coefficients.append(score)

    plt.plot(range(2, max_cluster + 1), silhouette_coefficients)
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    plt.show()

    index_min = np.argmin(silhouette_coefficients)
    logger.info(f"{index_min+2} clusters optimal")

    clust_serie = pd.Categorical(
        pd.Series(data=kmeans_list[index_min].labels_) + 1
    ).remove_categories([])

    return clust_serie, kmeans_list[index_min].cluster_centers_


def compute_Tm(temperatures, folding_fraction):
    """
    Compute the melting temperature using a sigmoidal curve fit.

    Parameters
    ----------
    temperatures : list
        List of temperatures.
    folding_fraction : list
        List of folding fraction.

    Returns
    -------
    Tm : float
        Melting temperature.
    """

    from scipy.optimize import curve_fit

    # Define the Four Parameter Logistic Regression (4PL)
    def sigmoidal_curve(x, A, B, C, D):
        return D + ((A - D) / (1.0 + ((x / C) ** B)))

    p0 = [1.0, 0.0, 340, 0.05]
    try:
        popt, pcov = curve_fit(sigmoidal_curve, temperatures, folding_fraction, p0=p0)
        logger.info(f"Melting Temperature (Tm): {popt[2]:.2f} K")
        return popt[2]
    except RuntimeError:
        logger.error("Error - curve_fit failed")
        return None


def plot_folding_fraction(
    df,
    col="RMSD (nm)",
    cutoff=0.18,
    label=None,
    start_time=0,
    time_ax_name=r"$Time\;(\mu s)$",
    recompute_temp_flag=True,
    temp_col="Aim Temp (K)",
    ref_temp=300.0,
):
    """
    Plot the fraction of folded protein as a function of the temperature.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the data to plot.
    col : str, optional
        Column to compute folding fraction. The default is "RMSD (nm)".
    cutoff : float, optional
        Cutoff value for the folding fraction. The default is 0.18 nm.
    label : str, optional
        Label of the plot. The default is None.
    start_time : float, optional
        Start time of the simulation. The default is 0 us.
    time_ax_name : str, optional
        Name of the time axis. The default is r"$Time\;(\mu s)$".
    recompute_temp_flag : bool, optional
        Recompute the temperature. The default is True.
    temp_col : str, optional
        Column with the temperature. The default is 'Aim Temp (K)'.
    ref_temp : float, optional
        Reference temperature. The default is 300.0.

    Returns
    -------
    Tm : float
        Melting temperature
    """

    df_time = df[(df[time_ax_name] > start_time)]

    fold_frac = compute_folding_fraction(df_time, col, cutoff)

    temp_list = df[temp_col].unique()
    temp_list.sort()

    if recompute_temp_flag:
        if ref_temp not in temp_list:
            ref_temp = temp_list[np.argmin(abs(temp_list - 300))]
        temp_list = recompute_temp(df, ref_temp=ref_temp)

    plt.plot(temp_list, fold_frac, label=label)
    plt.scatter(temp_list, fold_frac)
    plt.xlabel("Temperature (K)")
    plt.ylabel("fraction folded")
    plt.ylim((0, 1.0))

    return compute_Tm(temp_list, fold_frac)


def compute_folding_fraction(
    df, col="RMSD (nm)", cutoff=0.18, temp_col="Aim Temp (K)", temp_list=None
):
    """Compute the fraction of folded protein.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the data to plot.
    col : str, optional
        Column to compute folding fraction. The default is "RMSD (nm)".
    cutoff : float, optional
        Cutoff value for the folding fraction. The default is 0.18 nm.
    temp_col : str, optional
        Column with the temperature. The default is 'Aim Temp (K)'.
    temp_list : list, optional
        List of temperature to use. The default is None.

    Returns
    -------
    fold_frac : list
        List of folding fraction.
    """

    if temp_list is None:
        temp_list = df[temp_col].unique()
        temp_list.sort()
    else:
        temp_list.sort()
    fold_frac = []

    for temp in temp_list:

        local_df = df[(df[temp_col] == temp)]

        num_frame = len(local_df)
        num_cutoff = sum(local_df[col] < cutoff)

        if num_frame == 0:
            fold_frac.append(0)
        else:
            fold_frac.append(num_cutoff / num_frame)

    return fold_frac


def recompute_temp(df, ref_temp=300.0):

    # Compute real temp using Stinermann et al. JCTC 2015

    # Bi' = Bi(1 + ((Bref/Bi)**0.5 -1 ) ( Epw/(Epp+Epw)) )
    temp_list_traj = df["Aim Temp (K)"].unique()
    temp_list_traj.sort()

    new_temp_list = []
    ref_temp_index = list(temp_list_traj).index(ref_temp)

    inverseTemperatures = [
        1.0 / (unit.MOLAR_GAS_CONSTANT_R * t) for t in temp_list_traj
    ]

    for i, temp in enumerate(temp_list_traj):
        local_df = df[df["Aim Temp (K)"] == temp]
        # Epw
        Epp = abs(local_df["E solute scaled (kJ/mole)"].mean())
        Epw = abs(local_df["E solvent-solute (kJ/mole)"].mean())

        new_Bi = inverseTemperatures[i] * (
            1
            + (
                (inverseTemperatures[ref_temp_index] / inverseTemperatures[i]) ** 0.5
                - 1
            )
            * (Epw / (Epp + Epw))
        )
        new_temp = 1.0 / (unit.MOLAR_GAS_CONSTANT_R * new_Bi)
        # print(i, temp, new_temp)
        new_temp_list.append(new_temp)

    return new_temp_list


def compute_folding_fraction_RMSD(
    df,
    col="RMSD (nm)",
    temp_col="Aim Temp (K)",
    cutoff=0.18,
    start_time=0,
    time_ax_name=r"$Time\;(\mu s)$",
    ref_fold_frac=None,
    time_interval=2.0,
):

    if ref_fold_frac is None:
        df_time = df[(df[time_ax_name] > start_time)]
        ref_fold_frac = compute_folding_fraction(df_time, col, cutoff)
        ref_fold_frac = np.array(ref_fold_frac)

    temp_list = df["Aim Temp (K)"].unique()
    temp_list.sort()

    max_time = df[time_ax_name].max()

    RMSD_list = []
    time_list = []
    for i in range(int((max_time - start_time) / time_interval) + 1):
        # print(f"{i}  {start_time:.1f}  {start_time + (i + 1) * time_interval:.1f}")

        df_time = df[
            (df[time_ax_name] > start_time)
            & (df[time_ax_name] < start_time + (i + 1) * time_interval)
        ]

        fold_frac = compute_folding_fraction(
            df=df_time, col=col, cutoff=cutoff, temp_col=temp_col, temp_list=temp_list
        )

        fold_frac = np.array(fold_frac)

        RMSD_list.append(np.sum((ref_fold_frac - fold_frac) ** 2))
        time_list.append(start_time + (i + 1) * time_interval)

    return time_list, RMSD_list


def plot_folding_fraction_RMSD(
    df,
    col="RMSD (nm)",
    cutoff=0.18,
    label=None,
    start_time=0,
    time_ax_name=r"$Time\;(\mu s)$",
    ref_fold_frac=None,
    color=None,
    ls="-",
    s=20,
    alpha=1.0,
    time_interval=2.0,
):

    time_list, RMSD_list = compute_folding_fraction_RMSD(
        df=df,
        col=col,
        cutoff=cutoff,
        start_time=start_time,
        time_ax_name=time_ax_name,
        ref_fold_frac=ref_fold_frac,
        time_interval=time_interval
    )

    if color is None:
        plt.plot(time_list, RMSD_list, label=label, ls=ls, alpha=alpha)
        plt.scatter(time_list, RMSD_list, ls=ls, s=s, alpha=alpha)
    else:
        plt.plot(time_list, RMSD_list, label=label, color=color, ls=ls, alpha=alpha)
        plt.scatter(time_list, RMSD_list, color=color, ls=ls, s=s, alpha=alpha)

    plt.xlabel(time_ax_name)
    plt.ylabel("fraction folded RMSD")
    plt.legend()

    return time_list, RMSD_list


def plot_folding_fraction_convergence(
    df,
    col="RMSD (nm)",
    cutoff=0.18,
    label=None,
    start_time=0,
    time_ax_name=r"$Time\;(\mu s)$",
    recompute_temp_flag=False,
    ref_temp=300.0,
    time_interval=2.0,
):

    temp_list = df["Aim Temp (K)"].unique()
    temp_list.sort()

    if recompute_temp_flag:
        if ref_temp not in temp_list:
            ref_temp = temp_list[np.argmin(abs(temp_list - 300))]
        temp_list_plot = recompute_temp(df, ref_temp=ref_temp)
    else:
        temp_list_plot = temp_list

    max_time = df[time_ax_name].max()
    logger.info(max_time)

    for i in range(int((max_time - start_time) / time_interval) + 1):
        # print(f"{i}  {start_time:.1f}  {start_time + (i + 1) * time_interval:.1f}")

        df_time = df[
            (df[time_ax_name] > start_time)
            & (df[time_ax_name] < start_time + (i + 1) * time_interval)
        ]

        fold_frac = compute_folding_fraction(df_time, col, cutoff)

        plt.plot(
            temp_list_plot,
            fold_frac,
            label=f"{start_time}-{start_time + (i + 1) * time_interval:.1f}",
        )
        plt.scatter(temp_list_plot, fold_frac)

    plt.xlabel("Temperature (K)")
    plt.ylabel("fraction folded")
    plt.ylim((0, 1.0))
    plt.legend()


def plot_rung_occupancy(df, hue="group"):

    hue_list = df[hue].unique()

    for hue_val in hue_list:

        sim_df = df[df[hue] == hue_val]

        # for temp in enumerate(temp_list):
        temp_count = sim_df["Aim Temp (K)"].value_counts(normalize=True)
        temp_count = temp_count.sort_index()
        # print(len(temp_count))
        # print(temp_count)
        temp_count = temp_count / (1 / len(temp_count)) - 1
        # print(temp_count_rel)

        plt.plot(temp_count.index, temp_count.values, label=hue_val)
        plt.scatter(temp_count.index, temp_count.values)

    plt.xlabel("Temperature (K)")
    plt.ylabel(r"$\Delta$ Rung Occupancy")
    plt.legend(bbox_to_anchor=(1.01, 1.0))


def count_rmsd_transition(
    df,
    rmsd_fold=0.2,
    rmsd_unfold=0.4,
    dt=None,
    sim_name_col="sim",
    rmsd_col="RMSD (nm)",
    time_ax_name=r"$Time\;(\mu s)$",
):

    sim_list = df[sim_name_col].unique()
    trans_list = []

    for sim in sim_list:

        sim_df = df[df[sim_name_col] == sim]

        trans_num = 0
        fold_state = True if sim_df[rmsd_col].iloc[0] < rmsd_fold else False

        for rmsd in sim_df[rmsd_col]:
            if (rmsd < rmsd_fold) and not fold_state:
                trans_num += 1
                fold_state = True
            if (rmsd > rmsd_unfold) and fold_state:
                trans_num += 1
                fold_state = False
        if dt is None:
            dt_sim = sim_df[time_ax_name].iloc[1] - sim_df[time_ax_name].iloc[0]
            dt_steps = dt_sim * 1e6
        else:
            dt_sim = dt
            step_gap = sim_df["Step"].iloc[1] - sim_df["Step"].iloc[0]
            dt_steps = step_gap * dt_sim

        max_time = len(sim_df) * dt_steps / 1e6
        logger.info(
            f"Computed max time: {max_time:.3f}, df max time: {sim_df[time_ax_name].iloc[-1]:.3f}"
        )

        trans_freq = trans_num / max_time
        trans_list.append(trans_freq)

        logger.info(
            f"sim: {sim:20}  trans num={trans_num:6}  clust_freq = {trans_freq:.2f}/us, {max_time:.2f} {len(sim_df):.2f} "
        )

    df_trans = pd.DataFrame({"sim": sim_list, "trans": trans_list})
    return df_trans


def count_clust_transition(
    df, dt=None, sim_name_col="sim", clust_col="clust", time_ax_name=r"$Time\;(\mu s)$"
):

    sim_list = df[sim_name_col].unique()
    trans_list = []

    for sim in sim_list:

        sim_df = df[df[sim_name_col] == sim]

        last_clust = sim_df[clust_col].iloc[0]
        clust_num = 0

        for clust in sim_df[clust_col]:
            if (not np.isnan(clust)) and (clust != last_clust):
                clust_num += 1
                last_clust = clust
        if dt is None:
            dt_sim = sim_df[time_ax_name].iloc[1] - sim_df[time_ax_name].iloc[0]
            dt_steps = dt_sim * 1e6
        else:
            dt_sim = dt
            step_gap = sim_df["Step"].iloc[1] - sim_df["Step"].iloc[0]
            dt_steps = step_gap * dt_sim

        max_time = len(sim_df) * dt_steps / 1e6
        logger.info(
            f"Computed max time: {max_time:.3f}, df max time: {sim_df[time_ax_name].iloc[-1]:.3f}"
        )

        trans_freq = clust_num / max_time

        logger.info(
            f"sim: {sim:20}  clust num={clust_num:6}  clust_freq = {trans_freq:.2f}/us"
        )
        trans_list.append(trans_freq)

    df_trans = pd.DataFrame({"sim": sim_list, "trans": trans_list})
    return df_trans


def compare_weight_RMSD(
    df,
    x=r"$Time\;(\mu s)$",
    hue="sim",
    ener="new_pot",
    time_ax_name=r"$Time\;(\mu s)$",
    max_data=50000,
):

    local_df = filter_df(df, max_data)
    temp_list = local_df["Aim Temp (K)"].unique()
    temp_list.sort()

    # Compute overall avg
    temp_avg_dict = {}
    for temp in temp_list:
        temp_avg_dict[temp] = local_df[local_df["Aim Temp (K)"] == temp][ener].mean()
    logger.info(temp_avg_dict)

    group_list = local_df["group"].unique()

    rmsd_df = pd.DataFrame()

    for group in group_list:
        group_df = local_df[local_df["group"] == group]
        sim_list = group_df[hue].unique()

        logger.info(group, sim_list)

        for sim in sim_list:
            logger.info("    ", sim, group)
            sim_df = group_df[group_df[hue] == sim]
            sim_df = compute_moving_average(sim_df, ener=ener)

            compute_weight_RMSD(sim_df, final_weight_dict=temp_avg_dict, ener=ener)

            # sns.lineplot(
            #    data=sim_df,
            #    x=time_ax_name,
            #    y='Weight RMSD',
            #    label=sim,
            #    lw=2)

            local_weight_df = pd.DataFrame(
                {
                    "Weight RMSD": sim_df[sim_df["Weight RMSD"].notna()]["Weight RMSD"],
                    time_ax_name: sim_df[sim_df["Weight RMSD"].notna()][time_ax_name],
                    "group": group,
                    "sim": sim,
                }
            )

            local_weight_df[time_ax_name] = local_weight_df[time_ax_name].round(2)
            local_weight_df = local_weight_df.drop_duplicates(subset=[time_ax_name])

            sns.lineplot(
                data=local_weight_df, x=time_ax_name, y="Weight RMSD", label=sim, lw=0.5
            )

            rmsd_df = pd.concat([rmsd_df, local_weight_df])

            # print(sim_df[sim_df['Weight RMSD'].notna()])

    # Need to round time column, for a better averaging
    # rmsd_df[time_ax_name] = rmsd_df[time_ax_name].round(2)
    sns.lineplot(
        data=rmsd_df,
        x=time_ax_name,
        y="Weight RMSD",
        hue="group",
        # label=sim,
        lw=2,
    )

    plt.ylabel(r"RMSD $(KJ.mol^{-1})$")
    plt.title(r"Weights $f_i$ RMSD")

    return rmsd_df


def plot_energie_swap_convergence(
    df,
    ener_name="new_pot",
    lag_num=4,
    time_ax_name=r"$Time\;(\mu s)$",
    ylabel=r"$E_{p}$",
    split_graph=False,
    ci=95,
    avg_start=None,
):

    time_step = df.loc[1, time_ax_name] - df.loc[0, time_ax_name]
    temp_list = df["Aim Temp (K)"].unique()
    temp_list.sort()

    for temp_index in range(len(temp_list)):

        if avg_start is None:
            avg_ener = df[df["Aim Temp (K)"] == temp_list[temp_index]][ener_name].mean()
        else:
            avg_ener = df[
                (df["Aim Temp (K)"] == temp_list[temp_index])
                & (
                    (df["Temp Change index"] <= avg_start)
                    | (df["Temp Change index"] >= -avg_start)
                )
            ][ener_name].mean()

        # lag_num = 300
        df_local = df.loc[
            (df["Aim Temp (K)"] == temp_list[temp_index])
            & (df["Temp Change index"] <= lag_num)
            & (df["Temp Change index"] >= -lag_num),
            ["Temp Change index", ener_name],
        ]

        df_local.loc[:, "Time change"] = df_local["Temp Change index"] * time_step * 1e6
        df_local.loc.__setitem__(
            (slice(None), ("Time change")),
            df_local["Temp Change index"] * time_step * 1e6,
        )
        df_local_pos = df_local[df_local["Temp Change index"] > 0]
        df_local_neg = df_local[df_local["Temp Change index"] < 0]
        # df_local_neg.loc[:, "Time change"] = -1 * df_local_neg["Time change"]
        # To avoid warning, Replace by:
        df_local_neg.loc.__setitem__(
            (slice(None), ("Time change")), -1 * df_local_neg["Time change"]
        )

        if temp_index > 0:
            sns.lineplot(
                data=df_local_pos,
                x="Time change",
                y=ener_name,
                markers=True,
                dashes=False,
                ci=ci,
                label=f"E at {temp_list[temp_index]:.2f} K from {temp_list[temp_index-1]:.2f} K",
            )

        if temp_index < (len(temp_list) - 1):
            sns.lineplot(
                data=df_local_neg,
                x="Time change",
                markers=True,
                dashes=False,
                y=ener_name,
                ci=ci,
                label=f"E at {temp_list[temp_index]:.2f} K from {temp_list[temp_index+1]:.2f} K",
            )

        plt.xlabel("time (ps)")
        plt.ylabel(ylabel)
        plt.title(f"Energy after temperature change")

        plt.axhline(
            avg_ener,
            label=f"avg Epot at {temp_list[temp_index]:.2f} K",
            linestyle=":",
            c="gray",
        )
        plt.legend(bbox_to_anchor=(1.01, 1.0))

        if split_graph:
            plt.show()


def plot_energie_swap_convergence_diff(
    df,
    ener_name="new_pot",
    lag_num=4,
    time_ax_name=r"$Time\;(\mu s)$",
    ylabel=r"$E_{p}$",
    hue=None,
    color=None,
    label=r"$T_{m-1}$ update to $T_{m}$",
    errorbar=("ci", 95),
    avg_start=None,
):

    time_step = df.loc[1, time_ax_name] - df.loc[0, time_ax_name]
    temp_list = df["Aim Temp (K)"].unique()
    temp_list.sort()

    avg_ener_dict = {}

    print("Compute average energy at each temperature")

    for temp in temp_list:

        if avg_start is None:
            avg_ener = df[df["Aim Temp (K)"] == temp][ener_name].mean()
        else:
            avg_ener = df[
                (df["Aim Temp (K)"] == temp)
                & (
                    (df["Temp Change index"] <= avg_start)
                    | (df["Temp Change index"] >= -avg_start)
                )
            ][ener_name].mean()

        avg_ener_dict[temp] = avg_ener

    print(avg_ener_dict)

    print("Compute temperature change index")
    # Compute time change:
    time_change_list = []
    for temp_change_index in df["Temp Change index"]:
        time_change_list.append(temp_change_index * time_step * 1e6)

    df["Time change"] = time_change_list

    print("Compute Energy difference to average")
    energie_diff_list = []
    for temp, ener in zip(df["Aim Temp (K)"], df[ener_name]):

        ener_diff = ener - avg_ener_dict[temp]
        energie_diff_list.append(ener_diff)

    df["Energie Diff"] = energie_diff_list

    df_pos = df[(df["Temp Change index"] > 0) & (df["Temp Change index"] < lag_num)]
    df_neg = df[(df["Temp Change index"] < 0) & (df["Temp Change index"] > -lag_num)]

    df_neg.loc[:, "Time change"] *= -1

    print("Plot graph")

    sns.lineplot(
        data=df_pos,
        x="Time change",
        y="Energie Diff",
        hue=hue,
        markers=True,
        dashes=False,
        errorbar=errorbar,
        linestyle="--",
        color=color,
        label=None,
    )

    sns.lineplot(
        data=df_neg,
        x="Time change",
        y="Energie Diff",
        hue=hue,
        markers=True,
        dashes=False,
        errorbar=errorbar,
        color=color,
        label=label,
    )

    plt.xlabel(r"time $(ps)$")
    plt.ylabel(ylabel)
    plt.title(f"Energy evolution after temperature change")
    plt.legend(bbox_to_anchor=(1.01, 1.0))


def plot_energie_swap_distri_diff(
    df,
    lag_num_list,
    ener_name="new_pot",
    time_ax_name=r"$Time\;(\mu s)$",
    temp_index=1,
    ylabel=r"$E_{p}$",
    hue=None,
    bins=100,
    element="step",
    ci=95,
    avg_start=0,
):

    time_step = df.loc[1, time_ax_name] - df.loc[0, time_ax_name]
    temp_list = df["Aim Temp (K)"].unique()
    temp_list.sort()

    # Compute time change:
    time_change_list = []
    for temp_change_index in df["Temp Change index"]:
        time_change_list.append(round(temp_change_index * time_step * 1e6, 1))

    df["update (ps)"] = time_change_list

    print("Plot graph")

    fig, ax1 = plt.subplots()

    sns.histplot(
        df[
            (df["Aim Temp (K)"] == temp_list[temp_index])
            & (
                (df["Temp Change index"] > avg_start)
                | (df["Temp Change index"] < -avg_start)
            )
        ],
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=ener_name,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        element="step",
        color="black",
        label=temp_list[temp_index],
        line_kws={"linestyle": "-"},
        ax=ax1,
    )

    sns.histplot(
        df[
            (df["Aim Temp (K)"] == temp_list[temp_index - 1])
            & (
                (df["Temp Change index"] > avg_start)
                | (df["Temp Change index"] < -avg_start)
            )
        ],
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=ener_name,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        element="step",
        color="black",
        label=temp_list[temp_index - 1],
        line_kws={"linestyle": "--"},
        ax=ax1,
    )

    sns.histplot(
        df[
            (df["Aim Temp (K)"] == temp_list[temp_index + 1])
            & (
                (df["Temp Change index"] > avg_start)
                | (df["Temp Change index"] < -avg_start)
            )
        ],
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=ener_name,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        element="step",
        color="black",
        line_kws={"linestyle": ":"},
        label=temp_list[temp_index + 1],
        ax=ax1,
    )

    df_pos = df[
        (df["Aim Temp (K)"] == temp_list[temp_index])
        & (df["Temp Change index"].isin(lag_num_list))
    ]
    df_pos["update (ps)"] = pd.Categorical(df_pos["update (ps)"])

    sns.histplot(
        df_pos,
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=ener_name,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        hue="update (ps)",
        element="step",
        line_kws={"linestyle": "--"},
        ax=ax1,
    )

    df_neg = df.loc[
        (df["Aim Temp (K)"] == temp_list[temp_index])
        & (df["Temp Change index"].isin([-lag for lag in lag_num_list]))
    ]
    df_neg["update (ps)"] *= -1
    df_neg["update (ps)"] = pd.Categorical(df_neg["update (ps)"])

    sns.histplot(
        df_neg,
        stat="density",
        kde=True,
        bins=bins,
        fill=False,
        x=ener_name,
        common_norm=False,
        linewidth=1,
        alpha=0.3,
        hue="update (ps)",
        element="step",
        line_kws={"linestyle": ":"},
        ax=ax1,
    )

    sns.move_legend(ax1, "upper left", bbox_to_anchor=(1, 1), title_fontsize=9)
    # ax1.legend(bbox_to_anchor=(1.01, 1.0))

    legend = ax1.get_legend()
    # handles = legend.legendHandles # Deprecated
    handles = legend.legend_handles

    # Remove alpha in legend
    for h in handles:
        h.set_alpha(1.0)

    e_min, e_max = get_quant_min_max(
        df[df["Aim Temp (K)"] == temp_list[temp_index]][ener_name], quant=0.01
    )

    ymin, ymax = ax1.get_ylim()
    y_text = ymax - (ymax - ymin) / 20

    plt.annotate(f"{round(temp_list[temp_index-1],1)} K", xy=(e_min, y_text))
    plt.annotate(
        f"{round(temp_list[temp_index],1)} K",
        xy=(e_min + (e_max - e_min) / 2 - (e_max - e_min) / 10, y_text),
    )
    plt.annotate(
        f"{round(temp_list[temp_index+1],1)} K",
        xy=(e_max - (e_max - e_min) / 5.5, y_text),
    )

    plt.xlim(e_min, e_max)
    plt.xlabel(ylabel)
    plt.tight_layout()

    return ax1
