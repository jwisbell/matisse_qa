"""
Author: Jacob Isbell
Date: October 1, 2024

Script which contains primary logic for MATISSE QA.
"""

import json
from glob import glob
import numpy as np

from load_files import load_opd, load_phot_beams, load_raw_int, find_sof, load_spectrum
from phot_plots import plot_spectra
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
    opd_cutoff = int(cf["opd_cutoff"])
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
            print("Loading files ... ")
            data_dict = load_raw_int(files, verbose=verbose)
            opd_files = np.sort(glob(f"{d}/*OPD*00*.fits"))
            objcorr_files = glob(f"{d}/*OBJ_CORR*.fits")
            phot_files = np.sort(glob(f"{d}/*PHOT*00*.fits"))
            spectra_files = np.sort(glob(f"{d}/*SPECTRUM*00*.fits"))
            print("Files loaded successfully!")
        except:
            print(f"No files found in {d}")
            continue
        # TODO: using the tpl information, make subdirs in the output dir

        if len(data_dict["inst"]["tpl"]) == 0:
            continue

        tpl_start = data_dict["inst"]["tpl"]
        formatted_outdir = f"{output_dir}/{tpl_start.replace(':','_')}/"
        print(data_dict["inst"]["targname"])

        # make the directories
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

        # start with the opds so that a mask can be generated
        if True:
            print("Processing OPDs...")
            opd_dict = load_opd(opd_files, verbose=verbose)
            mask = plot_opd(
                data_dict,
                opd_dict,
                verbose=verbose,
                output_dir=formatted_outdir,
                opd_cutoff=opd_cutoff,
            )
        else:
            print("Error processing OPDs")

        # plot the visibilities
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

        try:
            print("Plotting fringe peaks ... ")
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
        except:
            print("Something went wrong while plotting group delay...")

        sof = ""
        try:
            sof = find_sof(data_dir, tpl_start)
            print(sof)
        except:
            print("SOF not found... ")

        if len(sof) > 0:  # True:
            do_waterfall(
                sof, output_dir=formatted_outdir, verbose=verbose, save_fig=save_fig
            )
        else:  # except (KeyError, FileNotFoundError) as e:
            print("SOF file missing, skipping waterfall plots, ")

        spectral_dict = load_spectrum(spectra_files, verbose=verbose)
        # print(spectral_dict, "test")
        try:
            plot_spectra(
                spectral_dict,
                verbose=verbose,
                save_fig=save_fig,
                output_dir=formatted_outdir,
            )
        except (KeyError, ValueError) as e:
            print(f"Something went wrong while creating the spectral plots... {e}")

        # phot_dict = load_phot_beams(phot_files, verbose=verbose)
        # then plot the photometries
        #
        #
        #
    print("Now showing all targets that have been processed!")
    get_obs(db_name)
