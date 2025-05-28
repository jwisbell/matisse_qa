"""
Author: Jacob Isbell
Date: October 3, 2024

Script which contains functions to plot closure phase data.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import triplet_idx_from_statriplet, bcd_flip, bcd_marker_dict, mask_wls

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True


normpa = mpl.colors.Normalize(vmin=-180, vmax=180)
normvis = mpl.colors.Normalize(vmin=0, vmax=1.1)


def plot_cphase(
    data_dict_all, output_dir: str = "/.", save_fig: bool = False, verbose: int = 0
):
    """
    data_dict: dictionary containing the following keys
    'cphase':{"t3phi":[],"t3phi_err":[], "u1":[],"u2":[], "v1":[],"v2":[],"wl_t3":[], "t3_sta":[], 'bcd':[]}

    Used to plot the closure phase
    """
    data_dict = data_dict_all["cphase"]
    targname = data_dict_all["inst"]["targname"]
    band = data_dict_all["inst"]["band"]

    fig1, axarr1 = plt.subplots(2, 5, figsize=(11, 8.5), width_ratios=[1, 1, 1, 1, 0.1])
    fig2, axarr2 = plt.subplots(1, 3, figsize=(8.5, 4.25), width_ratios=[1, 1, 0.1])
    axarr1[0, -1].axis("off")
    axarr1[1, -1].axis("off")
    axarr2[-1].axis("off")

    loop_names = ["2-3-4", "1-2-3", "1-2-4", "1-3-4"]

    collected_vals = [[], [], [], []]

    for i in range(0, len(data_dict["t3phi"])):
        sta = data_dict["t3_sta"][i]
        u1 = data_dict["u1"][i]
        v1 = data_dict["v1"][i]
        u2 = data_dict["u2"][i]
        v2 = data_dict["v2"][i]
        u3 = u2 - u1
        v3 = v2 - v1
        perimeter = (
            np.sqrt((u1 - u2) ** 2 + (v1 - v2) ** 2)
            + np.sqrt((u2 - u3) ** 2 + (v2 - v3) ** 2)
            + np.sqrt((u3 - u1) ** 2 + (v3 - v1) ** 2)
        )
        perimeter /= 3
        pa = np.degrees(np.arctan2(v1, u1))

        bcd = data_dict["bcd"][i]

        ydata = data_dict["t3phi"][i]
        yerr = data_dict["t3phi_err"][i]
        ydata = bcd_flip(ydata, sta, bcd)
        xdata = data_dict["wl_t3"][i] * 1e6

        # s = np.where(np.logical_and(xdata > 4.0, xdata < 4.5))[0]
        # if band == "N":
        #     # s = np.where(np.logical_and(xdata > 13.0, xdata < 7.5))[0]
        #     s = []  # no wavelengths to mask
        # ydata[s] = np.nan

        wl_mask = mask_wls(xdata, band)
        s = wl_mask
        ydata[~s] = np.nan

        zdata = [pa] * len(xdata)

        j = triplet_idx_from_statriplet(sta)
        marker = bcd_marker_dict[bcd]

        im = axarr1.flatten()[j % 4].scatter(
            xdata,
            ydata,
            s=2,
            c=zdata,
            norm=normpa,
            cmap="twilight",
            zorder=1,
            alpha=0.95,
        )

        axarr1.flatten()[j % 4].errorbar(
            xdata,
            ydata,
            yerr=yerr,
            ls="none",
            marker=".",
            color="k",
            alpha=0.1,
            zorder=0,
            errorevery=10,
        )
        axarr1.flatten()[j % 4].set_ylim([-180, 180])

        axarr1[1, j % 4].scatter(
            xdata,
            ydata,
            s=5,
            c=zdata,
            norm=normpa,
            cmap="twilight",
            zorder=1,
            marker=marker,
            alpha=0.5,
        )
        axarr1[1, j % 4].errorbar(
            xdata, ydata, yerr=yerr, ls="-", marker=".", color="k", alpha=0.25, zorder=0
        )

        axarr1[1, j % 4].set_ylim([-30, 30])
        axarr1[1, j % 4].set_title(loop_names[j])

        collected_vals[j].append(ydata)

        axarr2[0].scatter(
            perimeter / xdata,
            ydata,
            s=5,
            c=zdata,
            norm=normpa,
            cmap="twilight",
            zorder=1,
            marker=marker,
            alpha=0.5,
        )
        axarr2[0].errorbar(
            perimeter / xdata,
            ydata,
            yerr=yerr,
            ls="-",
            marker=".",
            color="k",
            alpha=0.25,
            zorder=0,
            errorevery=10,
        )

        axarr2[1].scatter(
            perimeter / xdata,
            ydata,
            s=5,
            c=zdata,
            norm=normpa,
            cmap="twilight",
            zorder=1,
            marker=marker,
            alpha=0.5,
        )
        axarr2[1].errorbar(
            perimeter / xdata,
            ydata,
            yerr=yerr,
            ls="-",
            marker=".",
            color="k",
            alpha=0.25,
            zorder=0,
            errorevery=10,
        )
        axarr2.flatten()[0].set_ylim([-180, 180])
        axarr2.flatten()[1].set_ylim([-30, 30])

    for i in range(0, 4):
        try:
            axarr1[0, i].errorbar(
                xdata,
                np.median(collected_vals[i], 0),
                yerr=np.std(collected_vals[i], 0) * 0,
                ls="--",
                marker=".",
                color="r",
                alpha=0.9,
                zorder=3,
                errorevery=10,
            )
            axarr1[1, i].errorbar(
                xdata,
                np.median(collected_vals[i], 0),
                yerr=np.std(collected_vals[i], 0) * 0,
                ls="--",
                marker=".",
                color="r",
                alpha=0.9,
                zorder=3,
                errorevery=10,
            )
        except:
            continue
    plt.colorbar(im, ax=axarr1[-1, -1], label="PA [deg]")
    plt.colorbar(im, ax=axarr2[-1], label="PA [deg]")

    axarr1[1, 0].set_ylabel("T3 [deg]")
    axarr1[1, 0].set_xlabel("Wavelength [micron]")

    axarr2[0].set_ylabel("T3 [deg]")
    axarr2[0].set_xlabel("Mean Leg Length Spatial Freq. [Mlambda]")

    fig1.tight_layout()
    fig2.tight_layout()

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/{targname}_{band}_t3_phi.png")
        fig2.savefig(f"{output_dir}/{targname}_{band}_t3_phi_perimeter_spatialfreq.png")

    if verbose > 1:
        plt.show()
    plt.close("all")

    _compute_t3_stats(data_dict_all)

    return None


def _compute_t3_stats(data_dict_all):
    # For now, just compute the exposure to exposure dispersion
    t3phi = [[] for _ in range(4)]
    stations = [0 for _ in range(4)]

    for i in range(len(data_dict_all["cphase"]["t3phi"])):
        sta = data_dict_all["cphase"]["t3_sta"][i]
        idx = triplet_idx_from_statriplet(sta)
        stations[idx] = sta
        t3phi[idx].append(data_dict_all["cphase"]["t3phi"][i])

    band = data_dict_all["inst"]["band"]
    wl_c = 3.5e-6
    delta_wl = 0.1e-6
    if band == "N":
        wl_c = 11e-6

    s = np.where(
        np.logical_and(
            np.array(data_dict_all["cphase"]["wl_t3"][0]) >= wl_c - delta_wl / 2,
            np.array(data_dict_all["cphase"]["wl_t3"][0]) <= wl_c + delta_wl / 2,
        )
    )[0]

    data_dict_all["qcparams"]["custom"]["cphase"] = {
        "stations": stations,
        "wavelength": [wl_c] * len(stations),
        "mean": [np.mean(np.mean(np.array(x), 0)[s]) for x in t3phi],
        "std": [np.mean(np.std(np.array(x), 0)[s]) for x in t3phi],
    }
