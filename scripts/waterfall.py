import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

raw_slice_locs = {"1": 334, "2": 355, "3": 378, "4": 402, "5": 425, "6": 445}


def testing(fname):
    x = fits.open(fname)

    print(x[1].header)

    ft = x[1].data["corrfluxreal1"][0] + 1j * x[1].data["corrfluximag1"][0]
    ift = np.fft.fftshift(np.fft.ifft2(ft))

    ift_shift = ift * np.exp(-1j * np.angle(ift))
    fft_shift = np.fft.fftshift(np.fft.fft2(ift_shift))

    _ = plt.figure()
    plt.imshow(
        np.abs(fft_shift),
        origin="lower",
    )
    plt.scatter(350, 63)
    plt.show()

    test = []
    test_gd = []
    for real, imag in zip(x[1].data["corrfluxreal1"], x[1].data["corrfluximag1"]):
        ft = real + 1j * imag
        gd = np.angle(ft, deg=False)

        # print(ampft, phaseft)
        ift = np.fft.fftshift(np.fft.ifft2(ft))

        ift_shift = ift * np.exp(-1j * np.angle(ift))
        fft_shift = np.fft.fftshift(np.fft.fft2(ift_shift))

        # for k,v in slice_locations.items():
        # slices[k].append(ft[63, v-width:v+width])
        test.append(
            fft_shift[
                fft_shift.shape[0] // 2,
                320 : 320 + 120,
            ],
        )
        test_gd.append(np.mean(np.diff(gd)[20:40], 0)[350] * 3.6)

    _ = plt.subplots()
    plt.imshow(np.sqrt(np.abs(test)), origin="lower")
    plt.show()

    _ = plt.subplots()
    plt.plot(test_gd)
    plt.show()

    plt.close("all")


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

    script, fname = argv

    # testing(fname)
    # basic_waterfall(
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T01:15:01.244.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T01:07:28.669.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/M.MATISSE.2022-02-15T13:08:58.613.fits",
    # )
    # basic_waterfall(
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T00:48:51.358.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T00:32:13.530.fits",
    #    "",
    # )

    # targ_files, sky_files = get_files_from_sof(
    #    "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T00:30:39.HAWAII-2RG.sof"
    # )

    targ_files, sky_files = _get_files_from_sof(
        "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T01:06:38.HAWAII-2RG.sof"
    )

    for tf in targ_files:
        _basic_waterfall(tf, sky_files)
