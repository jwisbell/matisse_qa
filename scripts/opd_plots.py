import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from utils import baseline_idx_from_stapair, bcd_color_dict


mpl.rcParams["font.family"] = "serif"
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["ytick.right"] = True


def _compute_mask_stats(masks, cutoff=3):
    flattened_masks = {k: [] for k, _ in masks.items()}
    max_time = 0
    counts = {}
    for k, v in masks.items():
        for val in v:
            for loc in val:
                flattened_masks[k].append(loc)
                if loc > max_time:
                    max_time = loc
                try:
                    counts[loc] += 1
                except KeyError:
                    counts[loc] = 1

    final_mask = []
    for idx, val in counts.items():
        if val >= cutoff:
            # print(idx, val)
            final_mask.append(idx)

    return final_mask


def plot_opd(
    data_dict, opd_dict, save_fig=True, verbose=0, output_dir=None, opd_cutoff=3
):
    fig1, axarr = plt.subplots(
        6,
        2,
        figsize=(8.5, 8.5),
        sharex=False,
        gridspec_kw={"width_ratios": [4, 1]},
        sharey=True,
    )
    fig2, axarr2 = plt.subplots(6, 1, figsize=(8.5, 8.5), sharex=True)
    targname = data_dict["inst"]["targname"]
    band = data_dict["inst"]["band"]
    yscale = 35
    y2scale = 15
    if band == "N":
        yscale = 50
        y2scale = 30

    baselines = {k: [] for k in range(6)}
    masks = {k: [] for k in range(6)}  # mask with bad values

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
                default_std = 4
                if band == "N":
                    default_std = 12
                mn = np.nanmean(opds[:, idx])
                # print(times[0], mn, std, "here")

                axarr[sta_idx, 0].scatter(
                    times,
                    opds[:, idx],
                    color=bcd_color_dict[bcd],
                    # label=bcd,
                    marker=".",
                )
                axarr[sta_idx, 0].fill_between(
                    times, mn - std, mn + std, color="lightgray", alpha=0.5
                )
                s = np.where(np.abs(opds[:, idx] - mn) > (2 * default_std))
                masks[sta_idx].append(times[s[0]])
                axarr[sta_idx, 0].scatter(
                    times[s], opds[:, idx][s], color="r", marker="x", alpha=0.25
                )
                axarr2.flatten()[sta_idx].scatter(
                    times[1:],
                    np.abs(np.diff(opds[:, idx])),
                    color=bcd_color_dict[bcd],
                    # label=bcd,
                )
    final_mask = _compute_mask_stats(masks, cutoff=opd_cutoff)
    for loc in final_mask:
        for idx in range(6):
            axarr[idx, 0].plot(
                [loc, loc], [-yscale, yscale], color="r", lw=3, alpha=0.35
            )

    for bcd, color in bcd_color_dict.items():
        axarr.flatten()[-2].scatter(-100, 2 * yscale, c=color, label=bcd)
        axarr2.flatten()[-1].scatter(-100, 2 * yscale, c=color, label=bcd)
    axarr2.flatten()[-1].legend(fontsize="small", ncol=6)
    axarr.flatten()[-2].legend(fontsize="small", ncol=6)

    bins = np.linspace(-yscale, yscale, 25)
    for k, v in baselines.items():
        vals = np.array([])
        for entry in v:
            vals = np.append(vals, entry)
        axarr[k, 1].hist(vals, bins=bins, orientation="horizontal", density=False)

    for ax in axarr[:, 0]:
        ax.set_ylim([-yscale, yscale])
        ax.set_xlim([0, None])

    for ax in axarr2.flatten():
        ax.set_ylim([0, y2scale])
        ax.set_xlim([0, None])

    axarr.flatten()[-2].set_xlabel("Time [MJD]")
    axarr.flatten()[-2].set_ylabel("OPD [micron]")

    fig1.suptitle(f"OPDs -- flagging with {opd_cutoff} simultaneous bad OPDs")

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
