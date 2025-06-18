import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from glob import glob
from scipy.ndimage import median_filter
import pandas as pd

from utils import bcd_color_dict_long as bcd_color_dict

lm_raw_slice_locs = {
    "7": 313,
    "1": 334,
    "2": 355,
    "3": 379,
    "4": 402,
    "5": 421,
    "6": 442,
}

n_raw_slice_locs = {
    "7": 235,
    "1": 251,
    "2": 270,
    "3": 289,
    "4": 305,
    "5": 322,
    "6": 342,
}

markers_dict = {
    "t1t2": "s",
    "t1t3": "o",
    "t1t4": "^",
    "t2t3": "d",
    "t2t4": "*",
    "t3t4": "v",
}


def _compute_obj_corr_vals(fname, band, wl=3.6):
    x = fits.open(fname)

    targname = x[0].header["eso obs targ name"]
    tpl = x[0].header["eso tpl start"]
    mjd0 = x[0].header["mjd-obs"]
    bcd_config = (
        f"{x[0].header['ESO INS BCD1 ID']}-{x[0].header['ESO INS BCD2 ID']}".lower()
    )
    is_chopping = x[0].header["eso iss chop st"] == "T"
    test_gd = {int(k) - 1: [] for k in lm_raw_slice_locs}
    test_fp = {int(k) - 1: [] for k in lm_raw_slice_locs}
    if band == "N":
        # TODO: handle differently
        test_gd = {int(k) - 1: [] for k in n_raw_slice_locs}
        test_fp = {int(k) - 1: [] for k in n_raw_slice_locs}

    times = []

    wl_min = 8.0
    wl_max = 13.0  # nband is flipped
    if band == "LM":
        wl_min = 3.0
        wl_max = 5.0

    nchannels = x[1].data["corrfluxreal1"][0].shape[0]
    all_wls = np.linspace(wl_min, wl_max, nchannels)
    wl_idx = np.argmin(np.abs(all_wls - wl))
    wl_w = 10

    if is_chopping:
        bcd_config += "_phot"

    group_delay_dict = {
        "out-out": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "out-out_phot": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "out-in": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-out": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-in": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-in_phot": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
    }

    fringe_peak_dict = {
        "out-out": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "out-out_phot": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "out-in": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-out": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-in": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-in_phot": {
            "total": [],
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
    }

    for real, imag, mjd in zip(
        x[1].data["corrfluxreal1"], x[1].data["corrfluximag1"], x[1].data["time"]
    ):
        ft = real + 1j * imag
        # fig = plt.figure()
        # plt.imshow(np.log10(np.abs(ft)), origin="lower")
        gd = np.angle(ft, deg=False)
        times.append(mjd)
        for k, v in test_gd.items():
            if band == "LM":
                if k == 0:
                    v.append(
                        np.mean(np.diff(gd)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            lm_raw_slice_locs[f"{k+1}"]
                        ]
                        * wl
                    )
                    test_fp[k].append(
                        np.max(np.abs(ft)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            lm_raw_slice_locs[f"{k+1}"]
                        ]
                    )
                else:
                    v.append(
                        np.mean(np.diff(gd)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            lm_raw_slice_locs[f"{k+1}"]
                        ]
                        * wl
                    )
                    test_fp[k].append(
                        np.max(np.abs(ft)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            lm_raw_slice_locs[f"{k+1}"]
                        ]
                    )
                    # plt.scatter(lm_raw_slice_locs[f"{k+1}"], 50)
            else:
                if k == 0:
                    v.append(
                        np.mean(np.diff(gd)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            n_raw_slice_locs[f"{k+1}"]
                        ]
                        * wl
                    )
                    test_fp[k].append(
                        np.max(np.abs(ft)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            n_raw_slice_locs[f"{k+1}"]
                        ]
                    )
                    # plt.scatter(n_raw_slice_locs[f"{k+1}"], wl_idx)
                else:
                    v.append(
                        np.mean(np.diff(gd)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            n_raw_slice_locs[f"{k+1}"]
                        ]
                        * wl
                    )
                    test_fp[k].append(
                        np.max(np.abs(ft)[wl_idx - wl_w : wl_idx + wl_w], 0)[
                            n_raw_slice_locs[f"{k+1}"]
                        ]
                    )
                    # plt.scatter(n_raw_slice_locs[f"{k+1}"], wl_idx)
        # plt.show()
        # plt.close()
        # exit()

    if bcd_config == "out-out":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t2t3"] = test_gd[2]
        group_delay_dict["out-out"]["t2t4"] = test_gd[3]
        group_delay_dict["out-out"]["t1t3"] = test_gd[4]
        group_delay_dict["out-out"]["t1t4"] = test_gd[5]
    elif bcd_config == "out-out_phot":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t2t3"] = test_gd[2]
        group_delay_dict["out-out"]["t2t4"] = test_gd[3]
        group_delay_dict["out-out"]["t1t3"] = test_gd[4]
        group_delay_dict["out-out"]["t1t4"] = test_gd[5]
    elif bcd_config == "out-in":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t1t3"] = test_gd[2]
        group_delay_dict["out-out"]["t1t4"] = test_gd[3]
        group_delay_dict["out-out"]["t2t3"] = test_gd[4]
        group_delay_dict["out-out"]["t2t4"] = test_gd[5]
    elif bcd_config == "in-out":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t2t4"] = test_gd[2]
        group_delay_dict["out-out"]["t2t3"] = test_gd[3]
        group_delay_dict["out-out"]["t1t4"] = test_gd[4]
        group_delay_dict["out-out"]["t1t3"] = test_gd[5]
    elif bcd_config == "in-in":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t1t4"] = test_gd[2]
        group_delay_dict["out-out"]["t1t3"] = test_gd[3]
        group_delay_dict["out-out"]["t2t4"] = test_gd[4]
        group_delay_dict["out-out"]["t2t3"] = test_gd[5]
    elif bcd_config == "in-in_phot":
        group_delay_dict["out-out"]["total"] = test_gd[6]
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t1t4"] = test_gd[2]
        group_delay_dict["out-out"]["t1t3"] = test_gd[3]
        group_delay_dict["out-out"]["t2t4"] = test_gd[4]
        group_delay_dict["out-out"]["t2t3"] = test_gd[5]

    if bcd_config == "out-out":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[2]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[3]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[4]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[5]
    elif bcd_config == "out-out_phot":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[2]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[3]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[4]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[5]
    elif bcd_config == "out-in":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[2]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[3]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[4]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[5]
    elif bcd_config == "in-out":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[2]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[3]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[4]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[5]
    elif bcd_config == "in-in":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[2]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[3]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[4]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[5]
    elif bcd_config == "in-in_phot":
        fringe_peak_dict["out-out"]["total"] = test_fp[6]
        fringe_peak_dict["out-out"]["t3t4"] = test_fp[0]
        fringe_peak_dict["out-out"]["t1t2"] = test_fp[1]
        fringe_peak_dict["out-out"]["t1t4"] = test_fp[2]
        fringe_peak_dict["out-out"]["t1t3"] = test_fp[3]
        fringe_peak_dict["out-out"]["t2t4"] = test_fp[4]
        fringe_peak_dict["out-out"]["t2t3"] = test_fp[5]
    """_ = plt.subplots()
    for k, v in test_gd.items():
        plt.plot(times, np.array(v) + 15 * int(k), label=k)
    plt.legend()
    plt.show()

    plt.close("all")
    """
    gd_dict = {
        "times": np.array(times),
        "group_delay": group_delay_dict,
        "fringe_peak": fringe_peak_dict,
        "target": targname,
        "bcd": bcd_config,
        "tpl": tpl,
    }
    return gd_dict


def argmax2d(arr):
    # wrapper for numpy argmax to give x,y coordinates of the max in a 2d array
    m = np.nanmax(arr)
    s = np.where(arr == m)
    return s[1], s[0]


def extract_fluxes(tel_arr, wl=3.7e-6):
    # extract the spectrum from the photometric channel and return the mean value
    # eventually return flux at a certain wl?
    return_vals = []
    w = 1
    for idx, t in enumerate(tel_arr):
        im = median_filter(t, 3)
        bias = np.nanmedian(im[:20, :20])
        slice = im[80, :]
        x = np.argmax(slice)
        spectrum = np.nanmedian(t[:, x - w : x + w], 1)
        val = spectrum - bias
        return_vals.append(val)

    return return_vals


def _basic_waterfall(
    targ,
    sky_files,
    output_dir: str = "/.",
    save_fig: bool = False,
    verbose: int = 0,
    band="LM",
):
    # using rudimentary bad pixel correction and sky subtraction, use the raw interferogram to
    # make a waterfall plot for the set of exposures in the targ file
    # edit to also show the photometric values over time
    xf = fits.open(targ)
    ext = "imaging_data"

    """
    for k in list(xf[0].header.keys()):
        if "chop" in k.lower():
            print(k, xf[0].header[k])
    """

    # load some basic info for displaying
    targname = xf[0].header["eso obs targ name"]
    tpl = xf[0].header["eso tpl start"]
    bcd = f"{xf[0].header['ESO INS BCD1 ID']}-{xf[0].header['ESO INS BCD2 ID']}"
    mjds = []
    is_chopping = xf[0].header["eso iss chop st"] == "T"

    # calculate the mean sky
    sky = np.mean([_calc_mean_sky(sf, band) for sf in sky_files], 0)

    slices = {k: [] for k in lm_raw_slice_locs if k != "7"}
    w = 10

    all_fluxes = [[], [], [], []]
    # load each interferogram and process it
    for i in range(len(xf[ext].data)):
        # print(f"{i}/{len(xf[ext].data)}")
        if band == "LM":
            interferogram = xf[ext].data[i][11]
            t1 = xf[ext].data[i][9]
            t2 = xf[ext].data[i][10]
            t3 = xf[ext].data[i][12]
            t4 = xf[ext].data[i][13]
            fluxes = extract_fluxes([t1, t2, t3, t4])
        else:
            interferogram = xf[ext].data["data5"][i]  # [11]
            t1 = xf[ext].data["data4"]
            t2 = xf[ext].data["data4"]
            t3 = xf[ext].data["data6"]
            t4 = xf[ext].data["data6"]
            fluxes = [1, 1, 1, 1]

        all_fluxes[0].append(fluxes[0])
        all_fluxes[1].append(fluxes[1])
        all_fluxes[2].append(fluxes[2])
        all_fluxes[3].append(fluxes[3])

        mjds.append(xf[ext].data[i][0])

        obs = np.array(interferogram - sky, dtype="float")

        # do rudimentary bad pixel correction
        obs[obs >= 65300] = np.nan
        s = np.where(np.isnan(obs))
        for x, y in zip(s[0], s[1]):
            obs[x, y] = np.nanmedian(obs[x - 2 : x + 2, y - 2 : y + 2])
        obs[np.isnan(obs)] = 0.0
        obs -= np.min(obs)

        ft = np.fft.fftshift(np.fft.fft2(obs))

        if band == "LM":
            for key, val in lm_raw_slice_locs.items():
                if key == "7" or key == "8":
                    continue
                xc, yc = argmax2d(np.abs(ft)[60 - 5 : 60 + 5, val - 10 : val + 10])
                slc = np.abs(ft)[yc + 60 - 5, val - w : val + w]
                slc -= np.mean(slc)
                slc /= np.max(slc)
                slices[key].append(slc.flatten())
        else:
            for key, val in n_raw_slice_locs.items():
                if key == "7" or key == "8":
                    continue
                xc, yc = argmax2d(np.abs(ft)[60 - 5 : 60 + 5, val - 10 : val + 10])
                w = 50
                slc = np.abs(ft)[yc + 60 - 5, val - w : val + w]
                slc -= np.mean(slc)
                slc /= np.max(slc)
                slices[key].append(slc.flatten())

    fig1, axarr = plt.subplots(1, 6, sharey=True, figsize=(10, 10))

    for idx, (key, value) in enumerate(slices.items()):
        axarr.flatten()[idx].imshow(
            np.array(value),
            origin="upper",
            interpolation="bilinear",
            cmap="viridis",
            vmin=-0.1,
            vmax=1,
        )
    fig1.suptitle(
        f"{targname} @ {tpl}\n(BCD:{bcd}  MJD: {mjds[0]:.4f}, chopping={is_chopping}) "
    )
    axarr.flatten()[0].set_ylabel("Time [increasing downward]")
    axarr.flatten()[0].set_xlabel("OPD")
    fig1.tight_layout()
    fig1.subplots_adjust(hspace=0.01, wspace=0.01)

    fig2, axarr2 = plt.subplots(2, 2, sharey=True, figsize=(6.5, 6.5), sharex=True)
    for idx, tel in enumerate(all_fluxes):
        # disp_arr = np.zeros((len(tel), len(tel[0])))
        for jdx, spectrum in enumerate(tel):
            axarr2.flatten()[idx].scatter(
                mjds[jdx], np.median(spectrum), color="k", marker="s"
            )
        axarr2.flatten()[idx].set_title(f"UT{idx+1}")
        # disp_arr[jdx] = spectrum
        # axarr2.flatten()[idx].imshow(disp_arr, origin="lower")

    # fig2.suptitle("Photometric flux")
    fig2.suptitle(
        f"{targname}_{band}_photometry_bcd{bcd}_ch{is_chopping}_mjd{f'{mjds[0]:.4f}'.replace('.','p')}"
    )

    if output_dir is not None and save_fig:
        fig1.savefig(
            f"{output_dir}/{targname}_{band}_waterfall_bcd{bcd}_ch{is_chopping}_mjd{f'{mjds[0]:.4f}'.replace('.','p')}.png"
        )
        fig2.savefig(
            f"{output_dir}/{targname}_{band}_photometry_bcd{bcd}_ch{is_chopping}_mjd{f'{mjds[0]:.4f}'.replace('.','p')}.png"
        )

    if verbose > 1:
        plt.show()

    plt.close("all")
    return None


def _get_files_from_sof(sofname):
    # extract all the sky and targ/calib files from the sof that the pipeline generates
    sky_files = []
    targ_files = []
    with open(sofname, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "CALIB" in line or "TARG" in line:
                targ_files.append(line.split()[0])
            elif "SKY_RAW" in line:
                sky_files.append(line.split()[0])

    return targ_files, sky_files


def _calc_mean_sky(sky_file, band="LM"):
    # calculate the mean sky from the sky exposures
    sf = fits.open(sky_file)
    ext = "imaging_data"

    sky_vals = []
    for i in range(len(sf[ext].data)):
        if band == "LM":
            sky_vals.append(sf[ext].data[i][11])
        else:
            sky_vals.append(sf[ext].data["data5"][i])
    return np.mean(sky_vals, 0)


def do_obj_corr_plots(
    files, band, wl=3.5e-6, tplstart="", output_dir="./", verbose=0, save_fig=False
):
    # plot the group delay from all the OBJ_CORR_FLUX files
    all_vals = []
    files = np.sort(files)
    for fname in files:
        all_vals.append(_compute_obj_corr_vals(fname, band, wl))

    all_fps = {k: {"sum": 0, "count": 0} for k in all_vals[0]["fringe_peak"]["out-out"]}
    for entry in all_vals:
        fp = entry["fringe_peak"]
        for key, value in fp["out-out"].items():
            all_fps[key]["count"] += len(value)
            all_fps[key]["sum"] += np.sum(value)

    means = {k: v["sum"] / v["count"] for k, v in all_fps.items()}

    offset_dict = {
        "t1t2": 0,
        "t1t3": 15,
        "t1t4": 30,
        "t2t3": 45,
        "t2t4": 60,
        "t3t4": 75,
    }

    fig1, ax = plt.subplots(figsize=(8.5, 8.5))
    fig2, bx = plt.subplots(figsize=(8.5, 8.5))
    fig3, cx = plt.subplots(figsize=(8.5, 8.5))
    max_time = 0
    fake_time = 0
    for entry in all_vals:
        times = entry["times"]
        if np.max(times) > max_time:
            max_time = np.max(times)
        gd_vals = entry["group_delay"]
        fp_vals = entry["fringe_peak"]
        target = entry["target"]
        bcd = entry["bcd"]
        ax.text(times[0], 80, f"{bcd}")
        bx.text(times[0], 11.5, f"{bcd}")

        for key, val in gd_vals["out-out"].items():
            xvals = range(fake_time, fake_time + len(fp_vals["out-out"][key]))
            if key == "total":
                cx.scatter(
                    xvals,
                    fp_vals["out-out"][key],
                    color=bcd_color_dict[bcd],
                )
                # fake_time += len(fp_vals["out-out"][key])
                continue
            offset = offset_dict[key]

            fp_yval = np.array(fp_vals["out-out"][key])
            # fp_yval -= np.min(fp_yval)
            fp_yval /= means[key]  # np.mean(fp_yval)
            fp_offset = offset_dict[key] / 15 * 2

            ax.scatter(
                # times,
                xvals,
                np.array(val) + offset,
                # color=colors_dict[key],
                color=bcd_color_dict[bcd],
                alpha=0.5,
                marker=markers_dict[key],
            )
            ax.plot(xvals, [offset] * len(xvals), "k--")
            ax.fill_between(
                xvals,
                [offset + 5] * len(xvals),
                [offset - 5] * len(xvals),
                color="gray",
                alpha=0.5,
            )

            bx.scatter(
                xvals,
                fp_yval + fp_offset,
                # color=colors_dict[key],
                color=bcd_color_dict[bcd],
                alpha=0.5,
                marker=markers_dict[key],
            )
            # bx.plot(times, [np.mean(fp_yval) + fp_offset] * len(times), "k--")
            bx.fill_between(
                xvals,
                [np.mean(fp_yval) + fp_offset - 2 * np.std(fp_yval)] * len(xvals),
                [np.mean(fp_yval) + fp_offset + 2 * np.std(fp_yval)] * len(xvals),
                color="gray",
                alpha=0.25,
            )
        fake_time += len(fp_vals["out-out"]["t1t2"])
    for key, value in offset_dict.items():
        ax.text(0, value, f"{key}")
        bx.text(0, value / 15 * 2 + 1, f"{key}")
    fig1.suptitle(f"{target}")
    fig2.suptitle(f"{target}")
    fig3.suptitle(f"{target}")
    ax.set_xlabel("MJD")
    ax.set_ylabel("Group Delay [um] (+offset)")
    bx.set_xlabel("MJD")
    bx.set_ylabel("Fringe Peak / Mean(FP) ")
    cx.set_xlabel("MJD")
    cx.set_ylabel("Zero Order Fringe Peak [counts]")
    fig1.tight_layout()
    fig2.tight_layout()
    plt.tight_layout()

    if output_dir is not None and save_fig:
        # fig1.savefig(f"{output_dir}/{target}_group_delay_lambda{wl}.png")
        fig2.savefig(f"{output_dir}/{target}_{band}_fringe_peak_lambda{wl}.png")
        fig3.savefig(f"{output_dir}/{target}_{band}_zero-order-fringe_lambda{wl}.png")

    if verbose > 1:
        plt.show()
    plt.close("all")

    # Save the group delay data
    df = pd.DataFrame.from_dict(all_vals)
    df.to_pickle(
        f"{output_dir}/../{target}_{tplstart.replace(":",'-').replace("_","-")}_{band}band_groupdelay_df.pkl"
    )

    return None


def do_waterfall(sof, output_dir, verbose, save_fig):
    # find the right files and then call the plotting function
    # the main script is responsible for finding the sof

    targ_files, sky_files = _get_files_from_sof(sof)
    band = "N"
    if "HAWAII" in sof:
        band = "LM"

    print(targ_files, sky_files)

    for tf in targ_files:
        _basic_waterfall(
            tf,
            sky_files,
            output_dir=output_dir,
            verbose=verbose,
            save_fig=save_fig,
            band=band,
        )


if __name__ == "__main__":
    from sys import argv

    script, fdir = argv

    fnames = np.sort(glob(f"{fdir}/OBJ_CORR_FLUX*.fits"))
    do_obj_corr_plots(fnames, "LM", 3.6, output_dir=None, verbose=2, save_fig=False)
    # do_obj_corr_plots(fnames, "N", 8.5, output_dir=None, verbose=2, save_fig=False)

    exit()

    # below here is working okay for now

    targ_files, sky_files = _get_files_from_sof(
        "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T01:06:38.HAWAII-2RG.sof"
    )

    for tf in targ_files:
        _basic_waterfall(tf, sky_files, verbose=2)
