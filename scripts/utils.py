import numpy as np

_ut_mapping = {
        '[32 33]':1,
        '[32 34]':4,
        '[32 35]':5,
        '[33 34]':2,
        '[33 35]':3,
        '[34 35]':0
    }


def baseline_idx_from_stapair(sta_pair):
    return _ut_mapping[str(np.sort(sta_pair) )]

def baseline_name_from_stapair(sta_pair):

    idx = _ut_mapping[str(np.sort(sta_pair) )]
    names = ['U1-U2', 'U1-U3', 'U1-U4', 'U2-U3', 'U2-U4', 'U3-U4']


    return names[idx]