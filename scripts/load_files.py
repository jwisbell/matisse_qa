"""
Author: Jacob Isbell
Date: October 1, 2024

Script which contains functions to load MATISSE fits files and turn them into dictionaries for later usage.
This will unify data formatting.

Specifically
1. Loads TARGET_RAW_INT/CALIB_RAW_INT
2. Loads OPD files
3. Loads photometric files
"""

import numpy as np
from astropy.io import fits
from glob import glob
import matplotlib.pyplot as plt


def load_raw_int(fnames, verbose: int = 0):
    """
    If fname is a single file load it, if fname is an array of files, load each.
    In both cases, return a data_dictionary (or pandas dataframe?)
    """
    if fnames is str:
        fnames = [fnames]

    # TODO: properly handle differential phase
    data = {
        "vis": {
            "cflux": [],
            "cflux_err": [],
            "u": [],
            "v": [],
            "wl_vis": [],
            "vis2": [],
            "vis2_err": [],
            "bcd": [],
            "vis2_sta": [],
            "diff_phase": [],
            "diff_phase_err": [],
        },
        "cphase": {
            "t3phi": [],
            "t3phi_err": [],
            "u1": [],
            "u2": [],
            "v1": [],
            "v2": [],
            "wl_t3": [],
            "t3_sta": [],
            "bcd": [],
        },
        "phot": {
            "phot_flux": [],
            "phot_flux_err": [],
            "phot_tel": [],
            "phot_wl": [],
            "bcd": [],
        },
        "inst": {"bcd": [], "tpl": "", "targname": "", "band": "N"},
        "qcparams": {"eso": {}, "custom": {}},
    }

    for f in fnames:
        band = "N"
        if "HAWAII" in f:
            band = "LM"

        data["inst"]["band"] = band
        x = fits.open(f)
        print(f"\t\t {f}")
        data["inst"]["tpl"] = x[0].header["eso tpl start"]
        data["inst"]["targname"] = (
            x[0].header["eso obs targ name"].strip().replace(" ", "").replace("_", "")
        )

        data["inst"]["tau0"] = np.mean(
            [
                float(x[0].header["ESO ISS AMBI TAU0 START"]),
                float(x[0].header["ESO ISS AMBI TAU0 END"]),
            ]
        )

        try:
            catalog_flux = float(x[0].header["HIERARCH PRO MDFC FLUX L"])
        except:
            # not a calibrator
            catalog_flux = -1
            pass
        try:
            catalog_diam = float(x[0].header["ESO PRO JSDC DIAMETER"])
        except:
            # not a calibrator
            catalog_diam = -1
            pass
        data["inst"]["cat_flux"] = catalog_flux
        data["inst"]["cat_diam"] = catalog_diam

        bcd = "oo"
        if "0002.fits" in f:
            bcd = "ii"
        elif "0003.fits" in f:
            bcd = "io"
        elif "0004.fits" in f:
            bcd = "oi"
        elif "0005.fits" in f:
            bcd = "oo_phot"
        elif "0006.fits" in f:
            bcd = "ii_phot"

        oi_key = "OI_VIS"
        for i, val in enumerate(x[oi_key].data["visamp"]):
            data["vis"]["cflux"].append(val)
            data["vis"]["diff_phase"].append(x[oi_key].data["visphi"][i])
            data["vis"]["diff_phase_err"].append(x[oi_key].data["visphierr"][i])
            data["vis"]["u"].append(x[oi_key].data["ucoord"][i])
            data["vis"]["v"].append(x[oi_key].data["vcoord"][i])
            data["vis"]["cflux_err"].append(x[oi_key].data["visamperr"][i])
            data["vis"]["wl_vis"].append(x["oi_wavelength"].data["eff_wave"])
            data["vis"]["vis2_sta"].append(x[oi_key].data["sta_index"][i])
            data["vis"]["bcd"].append(bcd)

        oi_key = "OI_VIS2"
        for i, val in enumerate(x[oi_key].data["vis2data"]):
            data["vis"]["vis2"].append(val)
            data["vis"]["vis2_err"].append(x[oi_key].data["vis2err"][i])

        oi_key = "OI_T3"
        for i, val in enumerate(x[oi_key].data["t3phi"]):
            data["cphase"]["t3phi"].append(val)
            data["cphase"]["t3phi_err"].append(x[oi_key].data["t3phierr"][i])
            data["cphase"]["u1"].append(x[oi_key].data["u1coord"][i])
            data["cphase"]["v1"].append(x[oi_key].data["v1coord"][i])
            data["cphase"]["u2"].append(x[oi_key].data["u2coord"][i])
            data["cphase"]["v2"].append(x[oi_key].data["v2coord"][i])
            data["cphase"]["t3_sta"].append(x[oi_key].data["sta_index"][i])
            data["cphase"]["wl_t3"].append(x["oi_wavelength"].data["eff_wave"])
            data["cphase"]["bcd"].append(bcd)

        oi_key = "OI_FLUX"
        try:
            for i, val in enumerate(x[oi_key].data["fluxdata"]):
                data["phot"]["phot_flux"].append(val)
                data["phot"]["phot_flux_err"].append(x[oi_key].data["fluxerr"][i])
                data["phot"]["phot_tel"].append(x[oi_key].data["sta_index"][i])
                data["phot"]["phot_wl"].append(x["oi_wavelength"].data["eff_wave"])
                data["phot"]["bcd"].append(bcd)
        except KeyError:
            print("No flux data found for", f)

        eso_qc = _parse_qc(x[0].header)
        data["qcparams"]["eso"] = eso_qc

    return data


def _parse_qc(header: dict) -> dict:
    parsed_qc = {}

    for key, value in header.items():
        if "qc" in key.lower():
            mod_key = key.replace("ESO QC", "").replace("DET1", "").strip()
            parsed_qc[mod_key] = value

    return parsed_qc


def load_opd(fnames, verbose=0):
    """
    If fname is a single file load it, if fname is an array of files, load each.
    In both cases, return a data_dictionary (or pandas dataframe?)
    """
    if fnames is str:
        fnames = [fnames]

    data = {
        k: {"time": [], "opd": [], "tels": []}
        for k in ["oo", "ii", "io", "oi", "oo_phot", "ii_phot"]
    }

    custom_time = 0
    for f in fnames:
        x = fits.open(f)
        opd = x[3].data["opd"]
        time = x[3].data["time"]
        # print(custom_time, time[-1])
        custom_time = time[-1] + 1
        # get mjd start
        # get exposure time
        # make custom time
        # test = x[3].data["time"]
        # custom_time = test

        sta_idx = x[3].data["sta_index"]
        bcd = "oo"
        bcd_config = (
            f"{x[0].header['ESO INS BCD1 ID']}-{x[0].header['ESO INS BCD2 ID']}".lower()
        )

        is_chopping = x[0].header["eso iss chop st"] == "T"

        if bcd_config == "in-out":
            bcd = "io"
        elif bcd_config == "out-in":
            bcd = "oi"
        elif bcd_config == "out-out":
            if is_chopping:
                bcd = "oo_phot"
            else:
                bcd = "oo"
        elif bcd_config == "in-in":
            if is_chopping:
                bcd = "ii_phot"
            else:
                bcd = "ii"

        sta_pairs = [(sta_idx[0][i], sta_idx[0][i + 1]) for i in range(0, 12, 2)]

        data[bcd]["time"].append(time)
        data[bcd]["opd"].append(opd)
        data[bcd]["tels"].append(sta_pairs)

    return data


def find_sof(fdir, tpl):
    # search in this fdir for a file with extension .sof that matches the tpl
    time = tpl.replace("_", ":")
    files = glob(f"{fdir}/*mat_raw_estimates*.sof")
    for file in files:
        if time in file:
            return file

    raise FileNotFoundError


def load_phot_beams(fnames, verbose=0):
    """
    If fname is a single file load it, if fname is an array of files, load each.
    In both cases, return a data_dictionary (or pandas dataframe?)
    """
    if fnames is str:
        fnames = [fnames]

    data = {
        k: {"time": [], "fluxes": [], "tels": []}
        for k in ["oo", "ii", "io", "oi", "oo_phot", "ii_phot"]
    }

    custom_time = 0
    for f in fnames:
        x = fits.open(f)

        test = x[1].data["datanosky1"][0]
        fig = plt.figure()
        plt.imshow(test, origin="lower")
        plt.show()

        exit()

        time = np.arange(custom_time, custom_time + len(opd), 1)  # x[3].data["time"]
        print(custom_time, time[-1])
        custom_time = time[-1] + 1

        sta_idx = x[3].data["sta_index"]

        bcd_config = (
            f"{x[0].header['ESO INS BCD1 ID']}-{x[0].header['ESO INS BCD2 ID']}".lower()
        )

        bcd = "oo"
        if "0002.fits" in f or bcd_config == "in-in":
            bcd = "ii"
        elif "0003.fits" in f or bcd_config == "in-out":
            bcd = "io"
        elif "0004.fits" in f or bcd_config == "out-in":
            bcd = "oi"
        elif "0005.fits" in f or (bcd_config == "out-out" and "0011.fits" in f):
            bcd = "oo_phot"
        elif "0006.fits" in f or (bcd_config == "in-in" and "0012.fits" in f):
            bcd = "ii_phot"

        sta_pairs = [(sta_idx[0][i], sta_idx[0][i + 1]) for i in range(0, 12, 2)]

        data[bcd]["time"].append(time)
        data[bcd]["opd"].append(opd)
        data[bcd]["tels"].append(sta_pairs)

    return data


def load_spectrum(fnames, verbose=0):
    """
    If fname is a single file load it, if fname is an array of files, load each.
    In both cases, return a data_dictionary (or pandas dataframe?)
    """
    if fnames is str:
        fnames = [fnames]

    data = {
        k: {"time": [], "spectra": [], "spectra_err": [], "tels": [], "wls": []}
        for k in ["oo", "ii", "io", "oi", "oo_phot", "ii_phot"]
    }

    custom_time = 0
    for f in fnames:
        x = fits.open(f)

        #'oi_wavelength'
        #'oi_flux'
        wls = x["oi_wavelength"].data["eff_wave"]
        fluxes = x["oi_flux"].data["fluxdata"]
        fluxerr = x["oi_flux"].data["fluxerr"]
        tels = x["oi_flux"].data["sta_index"]
        data["targname"] = x[0].header["eso obs targ name"].strip().replace(" ", "")
        # print(data["targname"], "targname")
        # print(tels)
        # print(fluxes.shape)

        is_chopping = x[0].header["eso iss chop st"] == "T"
        bcd_config = f"{x[0].header['ESO INS BCD1 ID'][0]}{x[0].header['ESO INS BCD2 ID'][0]}".lower()
        # TODO: figure out chopping status from header

        bcd = bcd_config
        if is_chopping:
            bcd += "_phot"

        # print(bcd)

        data[bcd]["time"].append(custom_time)
        data[bcd]["spectra"].append(fluxes)
        data[bcd]["spectra_err"].append(fluxerr)
        data[bcd]["tels"].append(tels)
        data[bcd]["wls"].append(wls)
        custom_time += 1
    return data
