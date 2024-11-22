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
    data_dict, output_dir: str = "/.", save_fig: bool = False, verbose: int = 0
):
    fig, axarr = plt.subplots(3, 2, sharex=False, figsize=(11, 8.5))
    mapping = {32: 0, 33: 1, 34: 2, 35: 3}

    tel_stats = {0: [], 1: [], 2: [], 3: []}

    normtime = mpl.colors.Normalize(
        vmin=0,
        vmax=12,  # np.max([np.max(x) for x in data_dict["oo"]["time"]])
    )
    for key, value in data_dict.items():
        fluxes = value["spectra"]
        fluxerr = value["spectra_err"]
        wls = np.array(value["wls"])
        tels = value["tels"]
        time = value["time"]

        color = bcd_color_dict[key]
        for idx in range(len(fluxes)):
            for jdx, t in enumerate(tels[idx]):
                axarr.flatten()[mapping[t]].plot(
                    wls[idx] * 1e6, fluxes[idx][jdx], color=color
                )
                axarr.flatten()[mapping[t]].set_xticklabels("")

                tel_stats[mapping[t]].append(fluxes[idx][jdx])
                im = axarr.flatten()[-1].scatter(
                    wls[idx] * 1e6,
                    fluxes[idx][jdx],
                    c=[time[idx]] * len(wls[idx]),
                    norm=normtime,
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
    print(t1.shape)
    axarr.flatten()[-2].set(
        xticks=np.linspace(0, t1.shape[1], 7),
        xticklabels=[
            f"{x:.1f}" for x in np.linspace(np.min(wls), np.max(wls), 7) * 1e6
        ],
    )
    axarr.flatten()[-2].set_ylabel("Time")
    axarr.flatten()[-2].set_xlabel("Wavelength")
    axarr.flatten()[-1].set_xlabel("Wavelength")

    h = t1.shape[0]
    spc = 5
    axarr.flatten()[-2].text(5, h * 1 - spc * 1 + 0, "UT1")
    axarr.flatten()[-2].text(5, h * 2 - spc * 1 + 0, "UT2")
    axarr.flatten()[-2].text(5, h * 3 - spc * 1 + 0, "UT3")
    axarr.flatten()[-2].text(5, h * 4 - spc * 1 + 0, "UT4")

    plt.show()

    return None
