import numpy as np
import pandas as pd

_ut_mapping = {
    "[32 33]": 1,
    "[32 34]": 4,
    "[32 35]": 5,
    "[33 34]": 2,
    "[33 35]": 3,
    "[34 35]": 0,
}

_ut_triplet_mapping = {
    "[32 33 34]": 1,
    "[32 33 35]": 2,
    "[33 34 35]": 0,
    "[32 34 35]": 3,
}


def apply_magic_numbers(sta_pair, bcd: str, band: str, vis2, wls):
    print(sta_pair)
    idx = baseline_idx_from_stapair(sta_pair)
    try:
        slope = bcd_magic_numbers[band][bcd][idx][0]
        intercept = bcd_magic_numbers[band][bcd][idx][1]
        divmult = bcd_magic_numbers[band][bcd][idx][2]
        magic_interp = slope * wls / 1e6 + intercept
    except KeyError:
        magic_interp = 1
        divmult = "m"

    if divmult == "m":
        return vis2 * magic_interp
    else:
        return vis2 / magic_interp


bcd_magic_numbers = {
    "LM": {
        "oo": {5: [0, 1, "m"], 2: [0, 1, "m"]},
        "oo_phot": {5: [0, 1, "m"], 2: [0, 1, "m"]},
        "oi": {
            4: [-112293.6863, 1.547707906, "d"],
            3: [-169956.9351, 1.818169518, "m"],
        },
        "io": {
            1: [-101396.8417, 1.634890426, "m"],
            0: [-158400.8854, 1.892069706, "d"],
        },
        "ii": {
            5: [-316934.6898, 2.653648283, "d"],
            2: [-20401.00783, 0.996200466, "d"],
        },
        "ii_phot": {
            5: [-316934.6898, 2.653648283, "d"],
            2: [-20401.00783, 0.996200466, "d"],
        },
    },
    "N": {
        "oo": {5: [0, 1], 2: [0, 1]},
        "oo_phot": {5: [0, 1], 2: [0, 1]},
        "oi": {
            4: [-52769.1174440545, 1.71871289158438],
            3: [-68604.0869050424, 1.90429960230199],
        },
        "io": {
            0: [-15776.0528733816, 1.2607215283921],
            1: [-31184.7101419972, 1.43912074965831],
        },
        "ii": {
            5: [-97803.7261954151, 2.30693829426419],
            2: [-31974.458991784, 1.39183795726365],
        },
    },
}
bcd_marker_dict = {
    "oo": "o",
    "ii": "s",
    "oi": "^",
    "io": "x",
    "oo_phot": "d",
    "ii_phot": "v",
}
bcd_color_dict = {
    "oo": "dodgerblue",
    "ii": "firebrick",
    "oi": "green",
    "io": "orange",
    "oo_phot": "black",
    "ii_phot": "mediumorchid",
}

bcd_color_dict_long = {
    "out-out": "dodgerblue",
    "in-in": "firebrick",
    "out-in": "green",
    "in-out": "orange",
    "out-out_phot": "black",
    "in-in_phot": "mediumorchid",
}


def baseline_idx_from_stapair(sta_pair):
    return _ut_mapping[str(np.sort(sta_pair))]


def baseline_name_from_stapair(sta_pair):
    idx = _ut_mapping[str(np.sort(sta_pair))]
    names = ["U1-U2", "U1-U3", "U1-U4", "U2-U3", "U2-U4", "U3-U4"]

    return names[idx]


def triplet_idx_from_statriplet(sta_triplet):
    return _ut_triplet_mapping[str(np.sort(sta_triplet))]


def bcd_flip(ydata, sta_triplet, bcd):
    if bcd == "oi":
        if str(np.sort(sta_triplet)) == "[32 33 34]":
            return -1 * ydata
        elif str(np.sort(sta_triplet)) == "[32 33 35]":
            return -1 * ydata
    elif bcd == "io":
        if str(np.sort(sta_triplet)) == "[33 34 35]":
            return -1 * ydata
        elif str(np.sort(sta_triplet)) == "[32 34 35]":
            return -1 * ydata
    elif bcd == "ii" or bcd == "ii_phot":
        return -1 * ydata

    return ydata


def mask_wls(wls, band):
    if band == "LM":
        s1 = np.where(np.logical_and(wls > 3.1, wls < 3.95))[0]
        s2 = np.where(np.logical_and(wls > 4.5, wls < 5))[0]
        kept = np.append(s1, s2)
        return np.array([i in kept for i in range(len(wls))])
    else:
        kept = np.where(np.logical_and(wls < 13.0, wls > 8.0))[0]
        return np.array([i in kept for i in range(len(wls))])


def export_dict_to_df(data_dict, outdir):
    dfvis = pd.DataFrame.from_dict(data_dict["vis"])
    dft3 = pd.DataFrame.from_dict(data_dict["cphase"])
    dfphot = pd.DataFrame.from_dict(data_dict["phot"])
    dfqc = pd.DataFrame.from_dict(data_dict["qcparams"])

    inst = data_dict["inst"]
    tpl = inst["tpl"]
    targ = inst["targname"]
    band = inst["band"]
    visfname = (
        f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_vis_df.pkl"
    )
    t3fname = f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_cphase_df.pkl"
    photfname = f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_phot_df.pkl"
    qcfname = (
        f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_qc_df.pkl"
    )

    gdfname = f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_groupdelay_df.pkl"
    opdfname = (
        f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_opd_df.pkl"
    )
    fnames = {
        "vis": visfname,
        "cphase": t3fname,
        "phot": photfname,
        "qc": qcfname,
        "gd": gdfname,
        "opd": opdfname,
    }
    parent_df = pd.DataFrame.from_dict({"files": fnames})
    parent_fname = f"{outdir}/{targ}_{tpl.replace(":","-").replace('_','-')}_{band}band_parent_df.pkl"

    # TODO: save these

    for df, fname in zip(
        [dfvis, dft3, dfphot, dfqc, parent_df],
        [visfname, t3fname, photfname, qcfname, parent_fname],
    ):
        df.to_pickle(fname)

    print("Wrote data to dataframes")
