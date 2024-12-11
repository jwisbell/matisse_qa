import numpy as np

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
