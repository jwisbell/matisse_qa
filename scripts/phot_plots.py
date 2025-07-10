"""
Author: Jacob Isbell
Date: November 22, 2024

Script which contains functions to plot the spectra
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import bcd_color_dict

mpl.rcParams["font.family"] = "serif"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True


normpa = mpl.colors.Normalize(vmin=-180, vmax=180)
normvis = mpl.colors.Normalize(vmin=0, vmax=1.1)
normtime = mpl.colors.Normalize(vmin=0, vmax=16)


def plot_spectra(
    data_dict, band, output_dir: str = "/.", save_fig: bool = False, verbose: int = 0
):
    fig, axarr = plt.subplots(3, 2, sharex=False, figsize=(11, 8.5))
    mapping = {32: 0, 33: 1, 34: 2, 35: 3}

    tel_stats = {0: [], 1: [], 2: [], 3: []}

    normtime = mpl.colors.Normalize(
        vmin=0,
        vmax=12,  # np.max([np.max(x) for x in data_dict["oo"]["time"]])
    )
    targname = data_dict["targname"]
    # band = data_dict["inst"]["band"]
    mywls = 0
    maxflux = 0

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

    for key, value in data_dict.items():
        if key == "targname":
            continue
        fluxes = value["spectra"]
        fluxerr = value["spectra_err"]
        wls = np.array(value["wls"])
        mywls = wls
        tels = value["tels"]
        time = value["time"]
        bcd = key

        color = bcd_color_dict[key]
        exposure_counter[bcd] += 1

        marker = marker_dict[exposure_counter[bcd]]

        for idx in range(len(fluxes)):
            for jdx, t in enumerate(tels[idx]):
                axarr.flatten()[mapping[t]].plot(
                    wls[idx] * 1e6,
                    fluxes[idx][jdx],
                    color=color,
                    label=key,
                    marker=marker,
                )
                axarr.flatten()[mapping[t]].set_xticklabels("")
                if np.max(fluxes[idx][jdx]) > maxflux:
                    maxflux = np.max(fluxes[idx][jdx])
                tel_stats[mapping[t]].append(fluxes[idx][jdx])
                im = axarr.flatten()[-1].scatter(
                    wls[idx] * 1e6,
                    fluxes[idx][jdx],
                    c=[time[idx]] * len(wls[idx]),
                    norm=normtime,
                    s=1,
                )

    plt.colorbar(im, ax=axarr.flatten()[-1], label="Time")
    t1 = np.array(tel_stats[0])[:, ::-1]
    # Create a column of zeros
    zeros_col = np.zeros((t1.shape[0], 5)) * np.nan
    zeros_row = np.zeros((5, t1.shape[1])) * np.nan
    # Append the column of zeros
    # t1 = np.hstack((t1, zeros_col))
    t1 = np.append(t1, zeros_row, axis=0)

    t2 = np.array(tel_stats[1])[:, ::-1]
    # t2 = np.hstack((t2, zeros_co))
    t2 = np.append(t2, zeros_row, axis=0)

    t3 = np.array(tel_stats[2])[:, ::-1]
    # t3 = np.hstack((t3, zeros_col))
    t3 = np.append(t3, zeros_row, axis=0)

    t4 = np.array(tel_stats[3])[:, ::-1]
    t4 = np.append(t4, zeros_row, axis=0)

    test = np.append(t1, t2, axis=0)
    test = np.append(test, t3, axis=0)
    test = np.append(test, t4, axis=0)

    axarr.flatten()[-2].imshow(
        test,
        origin="lower",
    )
    # print(data_dict)
    # print(mywls, "here")
    axarr.flatten()[-2].set(
        xticks=np.linspace(0, t1.shape[1], 7),
        xticklabels=[
            f"{x:.1f}" for x in np.linspace(np.min(mywls), np.max(mywls), 7) * 1e6
        ],
    )
    axarr.flatten()[-2].set_ylabel("Time")
    axarr.flatten()[-2].set_xlabel("Wavelength")
    axarr.flatten()[-1].set_xlabel("Wavelength")

    axarr.flatten()[0].set_title("UT1")
    axarr.flatten()[1].set_title("UT2")
    axarr.flatten()[2].set_title("UT3")
    axarr.flatten()[3].set_title("UT4")
    axarr.flatten()[0].legend(fontsize="small", ncol=2)

    axarr.flatten()[0].set_ylim(0, maxflux)
    axarr.flatten()[1].set_ylim(0, maxflux)
    axarr.flatten()[2].set_ylim(0, maxflux)
    axarr.flatten()[3].set_ylim(0, maxflux)
    axarr.flatten()[-1].set_ylim(0, maxflux)

    h = t1.shape[0]
    spc = 5
    axarr.flatten()[-2].text(5, h * 1 - spc * 1 + 0, "UT1")
    axarr.flatten()[-2].text(5, h * 2 - spc * 1 + 0, "UT2")
    axarr.flatten()[-2].text(5, h * 3 - spc * 1 + 0, "UT3")
    axarr.flatten()[-2].text(5, h * 4 - spc * 1 + 0, "UT4")

    plt.tight_layout()

    if output_dir is not None and save_fig:
        plt.savefig(f"{output_dir}/{targname}_{band}_singledish_spectra.pdf")

    if verbose > 1:
        plt.show()

    plt.close("all")

    return None
