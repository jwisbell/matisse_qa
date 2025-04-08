data_dict = {
    "info": {
        "stats": {"qc_data": [], "weather": []},
        "headers": [list_of_headers],
        "raw_files": [list_of_raw_files],
        "detector": "hawaii",
    },
    "data": {
        "cycles": [
            {
                "cycle_number": 1,
                "out-out": {
                    "chopped": {
                        "opd": {
                            "t1-t2": [],
                            "t1-t3": [],
                        },
                        "uv": {
                            "u": [],
                            "v": [],
                            "u1": [],
                            "u2": [],
                            "v1": [],
                            "v2": [],
                        },
                        "wl": wl_arr,
                        "vis2": {"vis2data": [], "vis2err": []},
                        "vis": {"visdata": [], "viserr": []},
                        "t3phi": {"t3phi": [], "t3phi_err": []},
                        "diff_phase": {"diff_phase": [], "diff_phase_err": []},
                        "spectra": {
                            "t1": [],
                            "t2": [],
                            "t3": [],
                            "t4": [],
                        },
                        "fringes": {"fringe_peak": []},
                    },
                    "unchopped": {...},
                },
                "in-in": {...},
                "in-out": {...},
                "out-in": {...},
            },
            {"cycle_number":2,
             ...}
        ]
    },
}
