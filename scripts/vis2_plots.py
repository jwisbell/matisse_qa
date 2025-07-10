"""
Author: Jacob Isbell
Date: October 1, 2024

Script which contains functions to plot visibility and correlated flux data.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import (
    apply_magic_numbers,
    baseline_idx_from_stapair,
    bcd_color_dict,
    bcd_magic_numbers,
    mask_wls,
    make_legend,
)

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True


normpa = mpl.colors.Normalize(vmin=-180, vmax=180)
normphi = mpl.colors.Normalize(vmin=-180, vmax=180)
normvis = mpl.colors.Normalize(vmin=0, vmax=1.1)


def plot_vis(
    data_dict_all,
    output_dir: str = "/.",
    save_fig: bool = False,
    verbose: int = 0,
    do_magic_numbers: bool = True,
):
    """
    data_dict: dictionary containing the following keys
    'vis':{"cflux":[],"cflux_err":[],"u":[], "v":[],"wl_vis":[],"vis2":[], "vis2_err":[],'bcd':[],'vis2_sta':[]

    Used to plot the correlated flux and visibility
    """
    data_dict = data_dict_all["vis"]
    targname = data_dict_all["inst"]["targname"]
    band = data_dict_all["inst"]["band"]
    fig1, axarr1 = plt.subplots(3, 3, figsize=(8.5, 8.5))
    fig2, axarr2 = plt.subplots(3, 3, figsize=(8.5, 8.5))
    fig4, axarr4 = plt.subplots(3, 3, figsize=(8.5, 8.5))
    fig3, axarr3 = plt.subplots(4, 3, figsize=(8.5, 8.5))

    bl_lengths = [[], [], [], [], [], []]
    bl_names = ["UT3-UT4", "UT1-UT2", "UT2-UT3", "UT2-UT4", "UT1-UT3", "UT1-UT4"]
    exposure_counter = {k: 0 for k in bcd_color_dict.keys()}
    marker_dict = {
        1: "s",
        2: "o",
        3: "^",
        4: "v",
        5: "d",
        6: "x",
        7: "+",
        8: ".",
        0: "none",
    }
    for i in range(len(data_dict["cflux"])):
        u = data_dict["u"][i]
        v = data_dict["v"][i]
        bl = np.sqrt(u**2 + v**2)
        sta = data_dict["vis2_sta"][i]
        idx = baseline_idx_from_stapair(sta)
        xdata = data_dict["wl_vis"][i] * 1e6
        ydata = data_dict["cflux"][i]

        bcd = data_dict["bcd"][i]
        yupper = 1.1
        exposure_counter[bcd] += 1
        marker = marker_dict[exposure_counter[bcd] // 6]

        wl_mask = mask_wls(xdata, band)
        s = wl_mask
        ydata[~s] = np.nan

        if band == "N":
            yupper = np.nanmax(data_dict["cflux"])
        normvis = mpl.colors.Normalize(vmin=0, vmax=yupper)

        bl_lengths[idx].append(bl)
        pa = np.degrees(np.arctan2(v, u))
        im = axarr1.flatten()[idx].scatter(
            xdata,
            ydata,
            s=2,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            fc=bcd_color_dict[bcd],
            # ec="k",
            marker=marker,
            zorder=1,
        )
        axarr1.flatten()[idx].errorbar(
            xdata,
            ydata,
            yerr=data_dict["cflux_err"][i],
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.0,
            errorevery=10,
        )
        axarr1.flatten()[idx].set_ylim([0, yupper])

        axarr1.flatten()[6].scatter(
            bl / xdata,
            ydata,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            fc=bcd_color_dict[bcd],
            # ec="k",
            marker=marker,
            zorder=1,
            s=2,
        )
        axarr1.flatten()[6].errorbar(
            bl / xdata,
            ydata,
            yerr=data_dict["cflux_err"][i],
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.1,
            errorevery=10,
        )
        axarr1.flatten()[6].set_ylim([0, yupper])

        im_uv1 = axarr1.flatten()[7].scatter(
            u / xdata[::4],
            v / xdata[::4],
            c=ydata[::4],
            norm=normvis,
            cmap="rainbow",
            zorder=1,
            alpha=0.5,
            s=2,
        )
        axarr1.flatten()[7].scatter(
            -u / xdata[::4],
            -v / xdata[::4],
            c=ydata[::4],
            norm=normvis,
            cmap="rainbow",
            zorder=1,
            alpha=0.5,
            s=2,
        )

        uvscale = 45
        if band == "N":
            uvscale = 20
        axarr1.flatten()[7].set_xlim([uvscale, -uvscale])
        axarr1.flatten()[7].set_ylim([-uvscale, uvscale])
        axarr1.flatten()[7].set_xlabel("u [Mlambda]")
        axarr1.flatten()[7].set_ylabel("v [Mlambda]")
        axarr1.flatten()[7].set_aspect("equal")

    for num, marker in marker_dict.items():
        if num == 0:
            continue
        axarr1[0, 0].scatter(
            xdata[0], 200, color="k", marker=marker, label=f"Exp. {num}"
        )
    axarr1[0, 0].legend(fontsize="x-small", ncol=4)

    bcd_sorted = [
        [[], [], [], [], [], []],
        [[], [], [], [], [], []],
        [[], [], [], [], [], []],
        [[], [], [], [], [], []],
        [[], [], [], [], [], []],
        [[], [], [], [], [], []],
    ]

    exposure_counter = {k: 0 for k in bcd_color_dict.keys()}
    marker_dict = {
        1: "s",
        2: "o",
        3: "^",
        4: "v",
        5: "d",
        6: "x",
        7: "+",
        8: ".",
        0: "none",
    }
    # vis 2 plots!
    for i in range(len(data_dict["vis2"])):
        u = data_dict["u"][i]
        v = data_dict["v"][i]
        bl = np.sqrt(u**2 + v**2)
        sta = data_dict["vis2_sta"][i]
        idx = baseline_idx_from_stapair(sta)

        xdata = data_dict["wl_vis"][i] * 1e6
        ydata = data_dict["vis2"][i]
        bcd = data_dict["bcd"][i]
        exposure_counter[bcd] += 1
        marker = marker_dict[exposure_counter[bcd] // 6]

        if do_magic_numbers:
            ydata = apply_magic_numbers(sta, bcd, band, ydata, xdata)

        yupper = 1.1

        wl_mask = mask_wls(xdata, band)
        s = wl_mask
        ydata[~s] = np.nan
        # ydata2[~s] = np.nan

        if band == "N":
            yupper = np.nanmax(data_dict["vis2"])
        normvis = mpl.colors.Normalize(vmin=0, vmax=yupper)

        bl_lengths[idx].append(bl)
        pa = np.degrees(np.arctan2(v, u))

        im2 = axarr2.flatten()[idx].scatter(
            xdata,
            ydata,
            s=2,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            fc=bcd_color_dict[bcd],
            # ec="k",
            marker=marker,
            zorder=1,
        )

        axarr2.flatten()[idx].errorbar(
            xdata,
            ydata,
            yerr=data_dict["vis2_err"][i],
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.0,
            errorevery=10,
        )
        axarr2.flatten()[idx].set_ylim([0, yupper])

        axarr2.flatten()[6].scatter(
            bl / xdata,
            ydata,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            fc=bcd_color_dict[bcd],
            marker=marker,
            # ec="k",
            zorder=1,
            s=2,
        )
        axarr2.flatten()[6].errorbar(
            bl / xdata,
            ydata,
            yerr=data_dict["vis2_err"][i],
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.1,
            errorevery=10,
        )
        axarr2.flatten()[6].set_ylim([0, yupper])

        im2_uv1 = axarr2.flatten()[7].scatter(
            u / xdata[::4],
            v / xdata[::4],
            c=ydata[::4],
            norm=normvis,
            cmap="rainbow",
            zorder=1,
            alpha=0.5,
            s=2,
        )
        axarr2.flatten()[7].scatter(
            -u / xdata[::4],
            -v / xdata[::4],
            c=ydata[::4],
            norm=normvis,
            cmap="rainbow",
            zorder=1,
            alpha=0.5,
            s=2,
        )

        uvscale = 45
        if band == "N":
            uvscale = 20
        axarr2.flatten()[7].set_xlim([uvscale, -uvscale])
        axarr2.flatten()[7].set_ylim([-uvscale, uvscale])
        axarr2.flatten()[7].set_xlabel("u [Mlambda]")
        axarr2.flatten()[7].set_ylabel("v [Mlambda]")
        axarr2.flatten()[7].set_aspect("equal")

        bcd = data_dict["bcd"][i]
        if bcd == "oo":
            bcd_sorted[0][idx].append(ydata)
        elif bcd == "ii":
            bcd_sorted[1][idx].append(ydata)
        elif bcd == "oi":
            bcd_sorted[2][idx].append(ydata)
        elif bcd == "io":
            bcd_sorted[3][idx].append(ydata)
        elif bcd == "oo_phot":
            bcd_sorted[4][idx].append(ydata)
        elif bcd == "ii_phot":
            bcd_sorted[5][idx].append(ydata)
    for num, marker in marker_dict.items():
        if num == 0:
            continue
        axarr2[0, 0].scatter(
            xdata[0], 200, color="k", marker=marker, label=f"Exp. {num}"
        )
    axarr2[0, 0].legend(fontsize="x-small", ncol=4)

    exposure_counter = {k: 0 for k in bcd_color_dict.keys()}
    marker_dict = {
        1: "s",
        2: "o",
        3: "^",
        4: "v",
        5: "d",
        6: "x",
        7: "+",
        8: ".",
        0: "none",
    }
    for i in range(len(data_dict["diff_phase"])):
        u = data_dict["u"][i]
        v = data_dict["v"][i]
        bl = np.sqrt(u**2 + v**2)
        sta = data_dict["vis2_sta"][i]
        idx = baseline_idx_from_stapair(sta)

        xdata = data_dict["wl_vis"][i] * 1e6
        ydata = data_dict["diff_phase"][i]
        yerr = data_dict["diff_phase_err"][i]
        bcd = data_dict["bcd"][i]
        exposure_counter[bcd] += 1
        marker = marker_dict[exposure_counter[bcd] // 6]

        wl_mask = mask_wls(xdata, band)
        s = wl_mask
        ydata[~s] = np.nan

        bl_lengths[idx].append(bl)

        im2 = axarr4.flatten()[idx].scatter(
            xdata,
            ydata,
            s=2,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            color=bcd_color_dict[bcd],
            marker=marker,
            zorder=1,
        )
        axarr4.flatten()[idx].errorbar(
            xdata,
            ydata,
            yerr=yerr,
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.0,
            errorevery=10,
        )
        axarr4.flatten()[idx].set_ylim([-180, 180])

        axarr4.flatten()[6].scatter(
            bl / xdata,
            ydata,
            # c=[pa] * len(data_dict["wl_vis"][i]),
            # norm=normpa,
            # cmap="twilight",
            color=bcd_color_dict[bcd],
            marker=marker,
            zorder=1,
            s=2,
        )
        axarr4.flatten()[6].errorbar(
            bl / xdata,
            ydata,
            yerr=yerr,
            zorder=0,
            color="k",
            ls="none",
            marker=".",
            alpha=0.0,
            errorevery=10,
        )
        axarr4.flatten()[6].set_ylim([-180, 180])

        im3_uv1 = axarr4.flatten()[7].scatter(
            u / xdata[::4],
            v / xdata[::4],
            c=ydata[::4],
            norm=normphi,
            cmap="coolwarm",
            zorder=1,
            alpha=0.5,
            s=2,
        )
        axarr4.flatten()[7].scatter(
            -u / xdata[::4],
            -v / xdata[::4],
            c=ydata[::4],
            norm=normphi,
            cmap="coolwarm",
            zorder=1,
            alpha=0.5,
            s=2,
        )

        uvscale = 45
        if band == "N":
            uvscale = 20
        axarr4.flatten()[7].set_xlim([uvscale, -uvscale])
        axarr4.flatten()[7].set_ylim([-uvscale, uvscale])
        axarr4.flatten()[7].set_xlabel("u [Mlambda]")
        axarr4.flatten()[7].set_ylabel("v [Mlambda]")
        axarr4.flatten()[7].set_aspect("equal")

    for num, marker in marker_dict.items():
        if num == 0:
            continue
        axarr4[0, 0].scatter(
            xdata[0], 200, color="k", marker=marker, label=f"Exp. {num}"
        )
    axarr4[0, 0].legend(fontsize="x-small", ncol=4)

    for i in range(6):
        label = "in-in"
        yupper = 1.1
        if band == "N":
            yupper = np.nanmax(bcd_sorted[0])
        if i != 0:
            label = None
        axarr3.flatten()[i].errorbar(
            np.median(bcd_sorted[0][i], 0),
            np.median(bcd_sorted[1][i], 0),
            yerr=np.std(bcd_sorted[1][i], 0),
            c=bcd_color_dict["ii"],
            zorder=1,
            marker="o",
            ls="none",
            alpha=0.75,
            label=label,
            errorevery=10,
        )
        for jdx in range(len(bcd_sorted[1][i])):
            axarr3.flatten()[i].errorbar(
                bcd_sorted[0][i][jdx],
                bcd_sorted[1][i][jdx],
                yerr=0,
                color="gray",
                zorder=0.5,
                marker="o",
                ls="none",
                alpha=0.25,
                errorevery=10,
            )
        label = "out-in"
        if i != 0:
            label = None
        axarr3.flatten()[i].errorbar(
            np.median(bcd_sorted[0][i], 0),
            np.median(bcd_sorted[2][i], 0),
            yerr=np.std(bcd_sorted[2][i], 0),
            c=bcd_color_dict["oi"],
            zorder=1,
            marker="s",
            ls="none",
            alpha=0.75,
            label=label,
            errorevery=10,
        )

        for jdx in range(len(bcd_sorted[2][i])):
            axarr3.flatten()[i].errorbar(
                bcd_sorted[0][i][jdx],
                bcd_sorted[2][i][jdx],
                yerr=0,
                color="gray",
                zorder=0.5,
                marker="o",
                ls="none",
                alpha=0.25,
                errorevery=10,
            )
        label = "in-out"
        if i != 0:
            label = None
        axarr3.flatten()[i].errorbar(
            np.median(bcd_sorted[0][i], 0),
            np.median(bcd_sorted[3][i], 0),
            yerr=np.std(bcd_sorted[3][i], 0),
            c=bcd_color_dict["io"],
            zorder=1,
            marker="^",
            ls="none",
            alpha=0.75,
            label=label,
        )
        for jdx in range(len(bcd_sorted[3][i])):
            axarr3.flatten()[i].errorbar(
                bcd_sorted[0][i][jdx],
                bcd_sorted[3][i][jdx],
                yerr=0,
                color="gray",
                zorder=0.5,
                marker="o",
                ls="none",
                alpha=0.25,
                errorevery=10,
            )

        try:
            axarr3.flatten()[i + 6].errorbar(
                np.median(bcd_sorted[4][i], 0),
                np.median(bcd_sorted[5][i], 0),
                yerr=np.std(bcd_sorted[4][i], 0),
                c="r",
                zorder=1,
                marker="o",
                ls="none",
                alpha=0.75,
            )
            for jdx in range(len(bcd_sorted[4][i])):
                axarr3.flatten()[i + 6].errorbar(
                    bcd_sorted[4][i + 6][jdx],
                    bcd_sorted[5][i + 6][jdx],
                    yerr=0,
                    color="gray",
                    zorder=0.5,
                    marker="o",
                    ls="none",
                    alpha=0.25,
                    errorevery=10,
                )
        except (IndexError, ValueError) as _:
            print("WARNING: No chopping data for this source.")

        axarr3.flatten()[i].plot([0, 1.1], [0, 1.1], c="k", ls="--", zorder=0)
        axarr3.flatten()[i + 6].plot([0, 1.1], [0, 1.1], c="k", ls="--", zorder=0)
        axarr3.flatten()[i].set_ylim([0, yupper])
        axarr3.flatten()[i].set_xlim([0, yupper])
        axarr3.flatten()[i].set_title("UNCHOPPED - " + bl_names[i])

        axarr3.flatten()[i + 6].set_ylim([0, yupper])
        axarr3.flatten()[i + 6].set_xlim([0, yupper])
        axarr3.flatten()[i + 6].set_title("CHOPPED - " + bl_names[i])
    axarr3.flatten()[0].legend()
    axarr3[3, 0].set_xlabel("Out-Out Vis2")
    axarr3[3, 0].set_ylabel("Other BCD Vis2")

    axarr1.flatten()[-1].axis("off")
    axarr2.flatten()[-1].axis("off")
    axarr4.flatten()[-1].axis("off")

    for i in range(6):
        axarr1.flatten()[i].set_title(
            f"{bl_names[i]}\nMean Proj. BL: {np.mean(bl_lengths[i]):.1f} m"
        )
        axarr2.flatten()[i].set_title(
            f"{bl_names[i]}\nMean Proj. BL: {np.mean(bl_lengths[i]):.1f} m"
        )
        axarr4.flatten()[i].set_title(
            f"{bl_names[i]}\nMean Proj. BL: {np.mean(bl_lengths[i]):.1f} m"
        )

    # axarr1[1,0].set_ylabel("Flux [cts]")
    axarr1[2, 0].set_ylabel("Coherent Vis [cts or None]")
    axarr1[1, 0].set_xlabel("Wavelength [micron]")
    axarr1[2, 0].set_xlabel("Spatial Frequency [MLambda]")
    axarr2[2, 0].set_ylabel("Incoherent Vis2")
    axarr2[1, 0].set_xlabel("Wavelength [micron]")
    axarr2[2, 0].set_xlabel("Spatial Frequency [MLambda]")

    axarr4[2, 0].set_ylabel("Differential Phase")
    axarr4[1, 0].set_xlabel("Wavelength [micron]")
    axarr4[2, 0].set_xlabel("Spatial Frequency [MLambda]")

    make_legend(axarr2[-1, -1])
    make_legend(axarr1[-1, -1])
    make_legend(axarr4[-1, -1])
    # axarr2.suptitle.set_text(f'')

    # plt.colorbar(
    #    im,
    #    ax=axarr1.flatten()[-1],
    #    orientation="horizontal",
    #    label="Position Angle [deg]",
    # )
    plt.colorbar(im_uv1, ax=axarr1.flatten()[-2], orientation="vertical", label="Vis")
    # plt.colorbar(
    #    im2,
    #    ax=axarr2.flatten()[-1],
    #    orientation="horizontal",
    #    label="Position Angle [deg]",
    # )
    plt.colorbar(im2_uv1, ax=axarr2.flatten()[-2], orientation="vertical", label="Vis")
    plt.colorbar(
        im3_uv1, ax=axarr4.flatten()[-2], orientation="vertical", label="Diff. Phase"
    )

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/{targname}_{band}_coher_vis.pdf")
        fig2.savefig(f"{output_dir}/{targname}_{band}_incoher_vis.pdf")
        fig3.savefig(f"{output_dir}/{targname}_{band}_bcd_compare_vis.pdf")
        fig4.savefig(f"{output_dir}/{targname}_{band}_diff_phase.pdf")

    if verbose > 1:
        plt.show()
    plt.close("all")

    _compute_vis2_stats(data_dict_all)

    return None


def _compute_vis2_stats(data_dict_all):
    # For now, just compute the exposure to exposure dispersion
    vis2 = [[] for _ in range(6)]
    dphase = [[] for _ in range(6)]
    cflux = [[] for _ in range(6)]
    stations = [0 for _ in range(6)]

    for i in range(len(data_dict_all["vis"]["vis2"])):
        sta = data_dict_all["vis"]["vis2_sta"][i]
        idx = baseline_idx_from_stapair(sta)
        stations[idx] = sta
        vis2[idx].append(data_dict_all["vis"]["vis2"][i])
        dphase[idx].append(data_dict_all["vis"]["diff_phase"][i])
        cflux[idx].append(data_dict_all["vis"]["cflux"][i])

    band = data_dict_all["inst"]["band"]
    wl_c = 3.5e-6
    delta_wl = 0.1e-6
    if band == "N":
        wl_c = 11e-6

    s = np.where(
        np.logical_and(
            np.array(data_dict_all["vis"]["wl_vis"][0]) >= wl_c - delta_wl / 2,
            np.array(data_dict_all["vis"]["wl_vis"][0]) <= wl_c + delta_wl / 2,
        )
    )[0]

    data_dict_all["qcparams"]["custom"]["vis2"] = {
        "stations": stations,
        "wavelength": [wl_c] * len(stations),
        "mean": [np.mean(np.mean(np.array(x), 0)[s]) for x in vis2],
        "std": [np.mean(np.std(np.array(x), 0)[s]) for x in vis2],
    }
    data_dict_all["qcparams"]["custom"]["dphase"] = {
        "stations": stations,
        "wavelength": [wl_c] * len(stations),
        "mean": [np.mean(np.mean(np.array(x), 0)[s]) for x in dphase],
        "std": [np.mean(np.std(np.array(x), 0)[s]) for x in dphase],
    }
    data_dict_all["qcparams"]["custom"]["visamp"] = {
        "stations": stations,
        "wavelength": [wl_c] * len(stations),
        "mean": [np.mean(np.mean(np.array(x), 0)[s]) for x in cflux],
        "std": [np.mean(np.std(np.array(x), 0)[s]) for x in cflux],
    }
