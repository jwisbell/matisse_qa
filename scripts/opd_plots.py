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
    fig1, axarr = plt.subplots(
        6,
        2,
        figsize=(11, 8.5),
        sharex=False,
        gridspec_kw={"width_ratios": [4, 1]},
        sharey=True,
    )
    fig2, axarr2 = plt.subplots(6, 1, figsize=(11, 8.5), sharex=True)
    targname = data_dict["inst"]["targname"]

    baselines = {k: [] for k in range(6)}

    for bcd in opd_dict.keys():
        for times, opds, stas in zip(
            opd_dict[bcd]["time"], opd_dict[bcd]["opd"], opd_dict[bcd]["tels"]
        ):
            """
            times = opd_dict[bcd]["time"]
            opds = opd_dict[bcd]["opd"]
            stas = opd_dict[bcd]["tel"]
            """
            for idx in range(6):
                baselines[idx].append(opds)
                sta_idx = baseline_idx_from_stapair(stas[idx])
                std = np.nanstd(opds[:, idx])
                mn = np.nanmean(opds[:, idx])
                print(times[0], mn, std, "here")
                axarr[sta_idx, 0].plot(times, opds[:, idx], color=bcd_color_dict[bcd])
                axarr[sta_idx, 0].fill_between(
                    times, mn - std, mn + std, color="lightgray", alpha=0.5
                )
                s = np.where(np.abs(opds[:, idx] - mn) > (2 * std))
                axarr[sta_idx, 0].scatter(times[s], opds[:, idx][s], color="r")
                axarr2.flatten()[sta_idx].plot(
                    times[1:], np.abs(np.diff(opds[:, idx])), color=bcd_color_dict[bcd]
                )

    for k, v in baselines.items():
        vals = np.array([])
        for entry in v:
            vals = np.append(vals, entry)
        axarr[k, 1].hist(vals, bins=25, orientation="horizontal", density=True)

    for ax in axarr[:, 0]:
        ax.set_ylim([-35, 35])

    for ax in axarr2.flatten():
        ax.set_ylim([0, 15])

    axarr.flatten()[-2].set_xlabel("Time [MJD]")
    axarr.flatten()[-2].set_ylabel("OPD [micron]")

    axarr2.flatten()[-1].set_xlabel("Time [MJD]")
    axarr2.flatten()[-1].set_ylabel("Delta OPD [micron]")
    plt.tight_layout()
    fig1.subplots_adjust(hspace=0.05, wspace=0.05)
    fig2.subplots_adjust(hspace=0.05, wspace=0.05)

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/{targname}_opd_values.png")
        fig2.savefig(f"{output_dir}/{targname}_diffopd_values.png")

    if verbose > 1:
        plt.show()

    plt.close("all")
