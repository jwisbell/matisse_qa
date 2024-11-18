"""
Author: Jacob Isbell
Date: October 1, 2024

Script which contains primary logic for MATISSE QA.
"""

import json
from glob import glob

from load_files import load_opd, load_raw_int, find_sof
from waterfall import do_obj_corr_plots, do_waterfall
from vis2_plots import plot_vis
from t3_plots import plot_cphase
from opd_plots import plot_opd
from db import create_database, get_obs, insert_observation
from sys import argv
from os import mkdir

if __name__ == "__main__":
    script, configfile = argv
    with open(configfile, "r") as f:
        cf = json.load(f)
    data_dir = cf["data_dir"]
    output_dir = cf["output_dir"]
    verbose = int(cf["verbose"])  # > 0
    save_fig = True

    db_name = f"{output_dir}/all_obs.db"
    # Initialize the database
    create_database(db_name)

    band = "LM"
    directories = glob(data_dir + "/*raw_estimates*.rb")
    print(directories, data_dir)
    for d in directories:
        if "HAWAII" in d:
            band = "LM"
        else:
            band = "N"

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
            exit()

        insert_observation(
            data_dict["inst"]["targname"],
            tpl_start,
            data_dict["inst"]["tau0"] * 1000,
            db_name,
        )

        try:
            print(f"Plotting visibilities ... ")
            plot_vis(
                data_dict,
                save_fig=save_fig,
                output_dir=formatted_outdir,
                verbose=verbose,
            )
            print("Visibilities plotted successfully!")
        except UnboundLocalError:
            print(f"Something wrong with {d} when plotting visibilities")
            continue

        try:
            print("Plotting closure phases ...")
            plot_cphase(
                data_dict,
                output_dir=formatted_outdir,
                save_fig=save_fig,
                verbose=verbose,
            )
            print("Closure phases plotted successfully!")
        except UnboundLocalError:
            print(f"Something wrong with {d} when plotting closure phases")
            continue

        opd_files = glob(f"{d}/*OPD*00*.fits")
        objcorr_files = glob(f"{d}/*OBJ_CORR*.fits")

        try:
            print("Plotting group delay ... ")
            mywl = 3.6
            if band == "N":
                mywl = 8.5
            do_obj_corr_plots(
                objcorr_files,
                band,
                output_dir=formatted_outdir,
                verbose=verbose,
                wl=mywl,
                save_fig=save_fig,
            )
            print("Closure phases plotted successfully!")
        except:
            print("Something went wrong while plotting group delay...")

        try:
            print("Processing OPDs...")
            opd_dict = load_opd(opd_files, verbose=verbose)
            plot_opd(data_dict, opd_dict, verbose=verbose, output_dir=formatted_outdir)
        except:
            print("Error processing OPDs")

        try:
            sof = find_sof(data_dir, tpl_start)
            do_waterfall(
                sof, output_dir=formatted_outdir, verbose=verbose, save_fig=save_fig
            )
        except (KeyError, FileNotFoundError) as e:
            print("SOF file missing, skipping waterfall plots, ", e)

    print("Now showing all targets that have been processed!")
    get_obs(db_name)
