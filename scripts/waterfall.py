import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from glob import glob

raw_slice_locs = {"1": 334, "2": 355, "3": 379, "4": 402, "5": 421, "6": 442}

colors_dict = {
    "t1t2": "blue",
    "t1t3": "red",
    "t1t4": "green",
    "t2t3": "mediumorchid",
    "t2t4": "orange",
    "t3t4": "brown",
}


def _compute_group_delay(fname, wl=3.6):
    x = fits.open(fname)

    targname = x[0].header["eso obs targ name"]
    tpl = x[0].header["eso tpl start"]
    bcd_config = (
        f"{x[0].header['ESO INS BCD1 ID']}-{x[0].header['ESO INS BCD2 ID']}".lower()
    )
    test_gd = {int(k) - 1: [] for k in raw_slice_locs}
    times = []

    group_delay_dict = {
        "out-out": {
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "out-in": {
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-out": {
            "t3t4": [],
            "t1t2": [],
            "t1t3": [],
            "t1t4": [],
            "t2t3": [],
            "t2t4": [],
        },
        "in-in": {
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
        gd = np.angle(ft, deg=False)
        times.append(mjd)
        for k, v in test_gd.items():
            v.append(np.mean(np.diff(gd)[20:40], 0)[raw_slice_locs[f"{k+1}"]] * wl)

    if bcd_config == "out-out":
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t2t3"] = test_gd[2]
        group_delay_dict["out-out"]["t2t4"] = test_gd[3]
        group_delay_dict["out-out"]["t1t3"] = test_gd[4]
        group_delay_dict["out-out"]["t1t4"] = test_gd[5]
    elif bcd_config == "out-in":
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t1t3"] = test_gd[2]
        group_delay_dict["out-out"]["t1t4"] = test_gd[3]
        group_delay_dict["out-out"]["t2t3"] = test_gd[4]
        group_delay_dict["out-out"]["t2t4"] = test_gd[5]
    elif bcd_config == "in-out":
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t2t4"] = test_gd[2]
        group_delay_dict["out-out"]["t2t3"] = test_gd[3]
        group_delay_dict["out-out"]["t1t4"] = test_gd[4]
        group_delay_dict["out-out"]["t1t3"] = test_gd[5]
    elif bcd_config == "in-in":
        group_delay_dict["out-out"]["t3t4"] = test_gd[0]
        group_delay_dict["out-out"]["t1t2"] = test_gd[1]
        group_delay_dict["out-out"]["t1t4"] = test_gd[2]
        group_delay_dict["out-out"]["t1t3"] = test_gd[3]
        group_delay_dict["out-out"]["t2t4"] = test_gd[4]
        group_delay_dict["out-out"]["t2t3"] = test_gd[5]
    """_ = plt.subplots()
    for k, v in test_gd.items():
        plt.plot(times, np.array(v) + 15 * int(k), label=k)
    plt.legend()
    plt.show()

    plt.close("all")
    """
    gd_dict = {
        "times": np.array(times),
        "vals": group_delay_dict,
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


def _basic_waterfall(
    targ, sky_files, output_dir: str = "/.", save_fig: bool = False, verbose: int = 0
):
    # using rudimentary bad pixel correction and sky subtraction, use the raw interferogram to
    # make a waterfall plot for the set of exposures in the targ file
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
    sky = np.mean([_calc_mean_sky(sf) for sf in sky_files], 0)

    slices = {k: [] for k in raw_slice_locs}
    w = 15

    # load each interferogram and process it
    for i in range(len(xf[ext].data)):
        interferogram = xf[ext].data[i][11]
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

        for key, val in raw_slice_locs.items():
            xc, yc = argmax2d(np.abs(ft)[60 - 5 : 60 + 5, val - 10 : val + 10])
            slc = np.abs(ft)[yc + 60 - 5, val - w : val + w]
            slc -= np.mean(slc)
            slc /= np.max(slc)
            slices[key].append(slc.flatten())

    _, axarr = plt.subplots(1, 6)
    for idx, (key, value) in enumerate(slices.items()):
        axarr.flatten()[idx].imshow(
            np.array(value),
            origin="upper",
            interpolation="bilinear",
            cmap="rainbow",
            vmin=-0.1,
            vmax=1,
        )
    plt.suptitle(
        f"{targname} @ {tpl}\n(BCD:{bcd} starting MJD {mjds[0]:.4f}, chopping={is_chopping}) "
    )
    axarr.flatten()[0].set_ylabel("Time [increasing downward]")
    axarr.flatten()[0].set_xlabel("OPD")
    plt.tight_layout()

    if output_dir is not None and save_fig:
        plt.savefig(
            f"{output_dir}/{targname}_waterfall_bcd{bcd}_ch{is_chopping}_mjd{f'{mjds[0]:.4f}'.replace('.','p')}.pdf"
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


def _calc_mean_sky(sky_file):
    # calculate the mean sky from the sky exposures
    sf = fits.open(sky_file)
    ext = "imaging_data"

    sky_vals = []
    for i in range(len(sf[ext].data)):
        sky_vals.append(sf[ext].data[i][11])

    return np.mean(sky_vals, 0)


def do_group_delay(files, wl, output_dir, verbose, save_fig):
    # plot the group delay from all the OBJ_CORR_FLUX files
    all_vals = []
    for fname in files:
        all_vals.append(_compute_group_delay(fname, wl))

    offset_dict = {
        "t1t2": 0,
        "t1t3": 15,
        "t1t4": 30,
        "t2t3": 45,
        "t2t4": 60,
        "t3t4": 75,
    }

    _ = plt.figure(figsize=(11, 8.5))
    max_time = 0
    for entry in all_vals:
        times = entry["times"]
        if np.max(times) > max_time:
            max_time = np.max(times)
        gd_vals = entry["vals"]
        target = entry["target"]
        bcd = entry["bcd"]
        plt.text(times[0], 90, f"{bcd}")
        for key, val in gd_vals["out-out"].items():
            offset = offset_dict[key]
            plt.plot(
                times,
                np.array(val) + offset,
                color=colors_dict[key],
                alpha=0.5,
            )
            plt.plot(times, [offset] * len(times), "k--")
            plt.fill_between(
                times,
                [offset + 5] * len(times),
                [offset - 5] * len(times),
                color="gray",
                alpha=0.5,
            )
    for key, value in offset_dict.items():
        plt.text(max_time, value, f"{key}")
    plt.suptitle(f"{target}")
    plt.xlabel("MJD")
    plt.ylabel("Group Delay [um] (+offset)")
    plt.tight_layout()

    if output_dir is not None and save_fig:
        plt.savefig(f"{output_dir}/{target}_group_delay_lambda{wl}.pdf")

    if verbose > 1:
        plt.show()
    plt.close("all")

    return None

    plt.close()


def do_waterfall(sof, output_dir, verbose, save_fig):
    # find the right files and then call the plotting function
    # the main script is responsible for finding the sof

    targ_files, sky_files = _get_files_from_sof(sof)

    for tf in targ_files:
        _basic_waterfall(
            tf, sky_files, output_dir=output_dir, verbose=verbose, save_fig=save_fig
        )


if __name__ == "__main__":
    from sys import argv

    script, fdir = argv

    fnames = np.sort(glob(f"{fdir}/OBJ_CORR_FLUX*.fits"))
    do_group_delay(fnames, 3.6, output_dir=None, verbose=2, save_fig=False)

    exit()

    # below here is working okay for now

    targ_files, sky_files = _get_files_from_sof(
        "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T01:06:38.HAWAII-2RG.sof"
    )

    for tf in targ_files:
        _basic_waterfall(tf, sky_files, verbose=2)
