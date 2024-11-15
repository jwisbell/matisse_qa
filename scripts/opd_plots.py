import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import baseline_idx_from_stapair, bcd_color_dict


mpl.rcParams["font.family"] = "serif"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True


def plot_opd(data_dict, opd_dict, save_fig=True, verbose=0, output_dir=None):
    fig1, axarr = plt.subplots(6, 1, figsize=(11, 8.5))
    fig2, axarr2 = plt.subplots(6, 1, figsize=(11, 8.5))
    targname = data_dict["inst"]["targname"]

    for bcd in opd_dict.keys():
        times = opd_dict[bcd]["time"]
        opds = opd_dict[bcd]["opd"]
        stas = opd_dict[bcd]["tel"]

        for idx in range(6):
            sta_idx = baseline_idx_from_stapair(stas[idx])
            std = np.std(opds[:, idx])
            mn = np.mean(opds[:, idx])
            axarr.flatten()[sta_idx].plot(
                times, opds[:, idx], color=bcd_color_dict[bcd]
            )
            axarr.flatten()[sta_idx].fill_between(
                times, mn - std, mn + std, color="lightgray", alpha=0.5
            )
            s = np.where(np.abs(opds[:, idx] - mn) > (2 * std))
            axarr.flatten()[sta_idx].scatter(times[s], opds[:, idx][s], color="r")
            axarr2.flatten()[sta_idx].plot(
                times[1:], np.abs(np.diff(opds[:, idx])), color=bcd_color_dict[bcd]
            )

    for ax in axarr.flatten():
        ax.set_ylim([-35, 35])

    for ax in axarr2.flatten():
        ax.set_ylim([0, 15])

    axarr.flatten()[-1].set_xlabel("Time [MJD]")
    axarr.flatten()[-1].set_ylabel("OPD [micron]")

    axarr2.flatten()[-1].set_xlabel("Time [MJD]")
    axarr2.flatten()[-1].set_ylabel("Delta OPD [micron]")

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/{targname}_opd_values.png")
        fig2.savefig(f"{output_dir}/{targname}_diffopd_values.png")

    if verbose > 1:
        plt.show()

    plt.close("all")
