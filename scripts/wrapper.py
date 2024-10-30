"""
Author: Jacob Isbell
Date: October 1, 2024

Script which contains primary logic for MATISSE QA.
"""

import json
from glob import glob

from load_files import load_raw_int
from vis2_plots import plot_vis
from t3_plots import plot_cphase
from sys import argv
from os import mkdir

if __name__ == "__main__":
    script, configfile = argv
    with open(configfile, "r") as f:
        cf = json.load(f)
    data_dir = cf["data_dir"]
    output_dir = cf["output_dir"]
    verbose = int(cf["verbose"])  # > 0

    directories = glob(data_dir + "/*raw_estimates*.rb")
    print(directories, data_dir)
    for d in directories:
        files = glob(f"{d}/*RAW_INT_00*.fits")
        print("\n\n" + "#" * 128)
        print(f"Working on {d}")
        try:
            print(f"Loading files ... ")
            data_dict = load_raw_int(files, verbose=verbose)
            print(f"Files loaded successfully!")
        except:
            print(f"No files found in {d}")
            continue
        # TODO: using the tpl information, make subdirs in the output dir

        if len(data_dict["inst"]["tpl"]) == 0:
            continue

        tpl_start = data_dict["inst"]["tpl"]
        formatted_outdir = f"{output_dir}/{tpl_start.replace(':','_')}/"
        print(data_dict["inst"]["targname"])
        try:
            mkdir(formatted_outdir)
        except FileExistsError:
            print(f"subdir {formatted_outdir} already exists")
        except FileNotFoundError:
            print(f"subdir {formatted_outdir} is not a valid path. Exiting...")
            sys.exit()

        try:
            print(f"Plotting visibilities ... ")
            plot_vis(
                data_dict,
                save_fig=True,
                output_dir=formatted_outdir,
                verbose=verbose,
            )
            print(f"Visibilities plotted successfully!")
        except UnboundLocalError:
            print(f"Something wrong with {d} when plotting visibilities")
            continue

        try:
            print("Plotting closure phases ...")
            plot_cphase(
                data_dict,
                output_dir=formatted_outdir,
                save_fig=True,
                verbose=verbose,
            )
            print(f"Closure phases plotted successfully!")
        except:
            print(f"Something wrong with {d} when plotting closure phases")
            continue
