import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    tplstart = data_dict["inst"]["tpl"]
    yscale = 35
    y2scale = 15
    if band == "N":
        yscale = 50
        y2scale = 30

    baselines = {k: [] for k in range(6)}
    masks = {k: [] for k in range(6)}  # mask with bad values
    start = 1e8
    position_tracker = {}
    testing_times = []

    for bcd in opd_dict.keys():
        for times, opds, stas in zip(
            opd_dict[bcd]["time"], opd_dict[bcd]["opd"], opd_dict[bcd]["tels"]
        ):
            """
            times = opd_dict[bcd]["time"]
            opds = opd_dict[bcd]["opd"]
            stas = opd_dict[bcd]["tel"]
            """
            if times[0] < start:
                start = times[0]

            for t in times:
                position_tracker[f"{t:.3f}"] = bcd
                testing_times.append(t)

            for idx in range(6):
                baselines[idx].append(opds)
                sta_idx = baseline_idx_from_stapair(stas[idx])
                std = np.nanstd(opds[:, idx])
                default_std = 4
                if band == "N":
                    default_std = 12
                mn = np.nanmean(opds[:, idx])

                axarr[sta_idx, 0].scatter(
                    times,
                    opds[:, idx],
                    color=bcd_color_dict[bcd],
                    # label=bcd,
                    marker=".",
                    zorder=1,
                )
                axarr[sta_idx, 0].fill_between(
                    times, mn - std, mn + std, color="lightgray", alpha=0.5
                )
                s = np.where(np.abs(opds[:, idx] - mn) > (2 * default_std))
                masks[sta_idx].append(times[s[0]])
                axarr[sta_idx, 0].scatter(
                    times[s],
                    opds[:, idx][s],
                    color="r",
                    marker="x",
                    alpha=0.25,
                    zorder=1,
                )
                axarr2.flatten()[sta_idx].scatter(
                    times[1:],
                    np.abs(np.diff(opds[:, idx])),
                    color=bcd_color_dict[bcd],
                    zorder=1,
                    # label=bcd,
                )
    counts_per_bcd = {k: 0 for k in bcd_color_dict.keys()}
    total_per_bcd = {k: 0 for k in bcd_color_dict.keys()}
    for val in position_tracker.values():
        total_per_bcd[val] += 1

    testing_times = np.sort(testing_times)
    test_diff = np.diff(testing_times)
    normal_difference = np.mean([np.mean(test_diff[:10]), np.mean(test_diff[-10:])])
    breaks = np.where(np.abs(test_diff) - normal_difference > 3 * normal_difference)[0]
    breaks = np.append([0], breaks)

    exposure_bins = [
        (testing_times[breaks[idx]], testing_times[breaks[idx + 1]])
        for idx in range(len(breaks) - 1)
    ] + [(testing_times[breaks[-1]], testing_times[-1])]

    exposure_mask_counts = [0 for _ in range(len(exposure_bins))]
    exposure_total_counts = [0 for _ in range(len(exposure_bins))]

    for time in testing_times:
        bin_idx = 0
        for idx, bin in enumerate(exposure_bins):
            if bin[0] <= time and time < bin[1]:
                bin_idx = idx
                break
        exposure_total_counts[bin_idx] += 1

    # debugging the exposure bin selection
    # for b in breaks:
    #     axarr[0, 0].scatter(
    #         testing_times[b + 1], 0, marker="s", s=50, color="k", zorder=100
    #     )

    final_mask = _compute_mask_stats(masks, cutoff=opd_cutoff)
    for loc in final_mask:
        counts_per_bcd[position_tracker[f"{loc:.3f}"]] += 1

        bin_idx = 0
        for bidx, bin in enumerate(exposure_bins):
            if bin[0] <= loc < bin[1]:
                bin_idx = bidx
                break
        exposure_mask_counts[bin_idx] += 1

        for idx in range(6):
            axarr[idx, 0].plot(
                [loc, loc], [-yscale, yscale], color="r", lw=3, alpha=0.15, zorder=0
            )
            axarr2[idx].plot(
                [loc, loc], [0, y2scale], color="r", lw=3, alpha=0.15, zorder=0
            )

    # percentages = {k: 0.0 for k in counts_per_bcd.keys()}
    percentages = [
        exposure_mask_counts[k] / 250 * 100 if exposure_total_counts[k] > 0 else 0.0
        for k in range(len(exposure_total_counts))
    ]

    # for k in counts_per_bcd.keys():
    #     if total_per_bcd[k] > 0:
    #         percentages[k] = counts_per_bcd[k] / total_per_bcd[k] * 100

    for bcd, color in bcd_color_dict.items():
        axarr.flatten()[-2].scatter(start, 0 * yscale, c=color, label=bcd)
        axarr2.flatten()[-1].scatter(start, 0 * yscale, c=color, label=bcd)
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
        ax.set_xlim([start, None])

    for ax in axarr2.flatten():
        ax.set_ylim([0, y2scale])
        ax.set_xlim([start, None])

    axarr.flatten()[-2].set_xlabel("Time [MJD]")
    axarr.flatten()[-2].set_ylabel("OPD [micron]")

    # fig1.suptitle(
    #     f"OPDs -- flagging with {opd_cutoff} simultaneous bad OPDs\n"
    #     + f"Unchopped flagged %:  I-I:{percentages['ii']:.1f},   O-I:{percentages['oi']:.1f},   I-O:{percentages['io']:.1f},   O-O:{percentages['oo']:.1f}\n"
    #     + f"Chopped flagged %: I-I_chop:{percentages['ii_phot']:.1f},    O-O_chop:{percentages['oo_phot']:.1f}"
    # )

    fig1.suptitle(
        f"OPDs -- flagging with {opd_cutoff} simultaneous bad OPDs\n"
        + f"flagged %: {[f'{x:.1f}%' for x in percentages]} ",
        fontsize="small",
    )

    axarr2.flatten()[-1].set_xlabel("Time [MJD]")
    axarr2.flatten()[-1].set_ylabel("Delta OPD [micron]")
    plt.tight_layout()
    fig1.subplots_adjust(hspace=0.05, wspace=0.05)
    fig2.subplots_adjust(hspace=0.05, wspace=0.05)

    if output_dir is not None and save_fig:
        fig1.savefig(f"{output_dir}/{targname}_{band}_opd_values.pdf")
        fig2.savefig(f"{output_dir}/{targname}_{band}_diffopd_values.pdf")

    if verbose > 1:
        plt.show()

    plt.close("all")

    df = pd.DataFrame.from_dict(opd_dict)
    df.to_pickle(
        f"{output_dir}/../{targname}_{tplstart.replace(":",'-').replace("_","-")}_{band}band_opd_df.pkl"
    )
    # print(f"Fraction flagged: {len(final_mask) / total_opd_measurements}")
    data_dict["qcparams"]["custom"]["frac_flagged_opd"] = percentages
