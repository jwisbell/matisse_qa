import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

slice_locations = {"1": 355, "2": 395}
raw_slice_locs = {"1": 334, "2": 355, "3": 378, "4": 402, "5": 425, "6": 445}


def testing(fname):
    x = fits.open(fname)

    print(x[1].header)

    ft = x[1].data["corrfluxreal1"][0] + 1j * x[1].data["corrfluximag1"][0]
    ift = np.fft.fftshift(np.fft.ifft2(ft))

    ift_shift = ift * np.exp(-1j * np.angle(ift))
    fft_shift = np.fft.fftshift(np.fft.fft2(ift_shift))

    ft_test = ft / np.exp(1j * np.angle(ft))
    fig = plt.figure()
    plt.imshow(
        np.abs(fft_shift),
        origin="lower",
    )
    plt.scatter(350, 63)
    plt.show()

    slices = {k: [] for k in slice_locations}
    width = 20
    test = []
    test_gd = []
    for real, imag in zip(x[1].data["corrfluxreal1"], x[1].data["corrfluximag1"]):
        ft = real + 1j * imag
        gd = np.angle(ft, deg=False)

        ft_test = ft * np.exp(-1j * np.angle(ft))

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

    fig = plt.subplots()
    plt.imshow(np.sqrt(np.abs(test)), origin="lower")
    plt.show()

    fig = plt.subplots()
    plt.plot(test_gd)
    plt.show()

    plt.close("all")


def argmax2d(arr):
    # wrapper for numpy argmax to give x,y coordinates of the max in a 2d array
    m = np.nanmax(arr)
    s = np.where(arr == m)
    return s[1], s[0]


def test_raw(targ, sky):
    xf = fits.open(targ)
    ext = "imaging_data"

    sky = calc_mean_sky(sky)

    slices = {k: [] for k in raw_slice_locs}
    w = 15

    fig = plt.figure()

    for i in range(len(xf[ext].data)):
        interferogram = xf[ext].data[i][11]

        obs = np.array(interferogram - sky, dtype="float")
        obs[obs >= 65300] = np.nan

        s = np.where(np.isnan(obs))
        for x, y in zip(s[0], s[1]):
            obs[x, y] = np.nanmedian(obs[x - 2 : x + 2, y - 2 : y + 2])
        obs[np.isnan(obs)] = 0.0
        obs -= np.min(obs)

        plt.scatter(i, np.nanmean(obs))

        # if np.nanmean(obs) < 50:
        #    continue

        # plt.plot(i, np.nanmean(obs))
        ft = np.fft.fftshift(np.fft.fft2(obs))
        # ft -= np.mean(ft)

        # plt.scatter(i, np.nanmean(np.abs(ft)))
        for key, val in raw_slice_locs.items():
            # plt.scatter(val, 60)
            xc, yc = argmax2d(np.abs(ft)[60 - 5 : 60 + 5, val - 10 : val + 10])
            # plt.scatter(xc + val - 10, yc + 60 - 5, marker="x")
            slc = np.abs(ft)[yc + 60 - 5, val - w : val + w]
            slc -= np.mean(slc)
            slc /= np.max(slc)
            slices[key].append(slc.flatten())

    plt.close()
    plt.show()

    fig, axarr = plt.subplots(1, 6)
    for idx, (key, value) in enumerate(slices.items()):
        print(value)
        axarr.flatten()[idx].imshow(
            np.array(value),
            origin="upper",
            interpolation="bilinear",
            cmap="rainbow",
            vmin=-0.1,
            vmax=1,
        )
    axarr.flatten()[0].set_ylabel("Time [increasing downward]")
    axarr.flatten()[0].set_xlabel("OPD")
    plt.show()

    plt.close("all")

    # print(x[ext].data["tdim1"])


def get_files_from_sof(sofname):
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


def calc_mean_sky(sky_file):
    sf = fits.open(sky_file)
    ext = "imaging_data"

    sky_vals = []
    for i in range(len(sf[ext].data)):
        sky_vals.append(sf[ext].data[i][11])

    return np.mean(sky_vals, 0)


if __name__ == "__main__":
    from sys import argv

    script, fname = argv

    # testing(fname)
    # test_raw(
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T01:15:01.244.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T01:07:28.669.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/M.MATISSE.2022-02-15T13:08:58.613.fits",
    # )
    # test_raw(
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T00:48:51.358.fits",
    #    "/Users/jwisbell/Documents/matisse/cena/rawdata/MATIS.2022-04-24T00:32:13.530.fits",
    #    "",
    # )

    # targ_files, sky_files = get_files_from_sof(
    #    "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T00:30:39.HAWAII-2RG.sof"
    # )

    targ_files, sky_files = get_files_from_sof(
        "/Users/jwisbell/Documents/matisse/cena/processeddata/lm_vis2_b5_v2/mat_raw_estimates.2022-04-24T01:06:38.HAWAII-2RG.sof"
    )

    for tf in targ_files:
        test_raw(tf, sky_files[0])
