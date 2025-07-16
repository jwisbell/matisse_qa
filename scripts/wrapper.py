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
from make_pdf import merge_pdfs, create_title_page
from utils import export_dict_to_df, cleanup_plots
from waterfall import do_obj_corr_plots, do_waterfall
from vis2_plots import plot_vis
from t3_plots import plot_cphase
from opd_plots import plot_opd
from sys import argv
from os import mkdir, _exit
from concurrent.futures import ProcessPoolExecutor


def process_directory(x):
    # TODO: up until N_cores, do each directory in parallel
    # script, configfile = argv
    # with open(configfile, "r") as f:
    #    cf = json.load(f)
    # data_dir = cf["data_dir"]
    # output_dir = cf["output_dir"]
    # verbose = int(cf["verbose"])  # > 0
    # opd_cutoff = int(cf["opd_cutoff"])
    # skip_plotting = cf["tables_only"]
    save_fig = True
    d, data_dir, output_dir, verbose, opd_cutoff, skip_plotting = x
    logfile = f"{output_dir}/qc_log_{d.split('/')[-1]}.txt"
    with open(logfile, "w") as f:
        f.write(f"Log for {d}")
        print(f"Making log file at {logfile}")

        try:
            if "HAWAII" in d:
                band = "LM"
            else:
                band = "N"

            files = glob(f"{d}/*RAW_INT_00*.fits")
            f.write("\n\n" + "#" * 128)
            f.write(f"Working on {d}\n")
            try:
                f.write("Loading files ... \n")
                data_dict = load_raw_int(files, verbose=verbose)
                opd_files = np.sort(glob(f"{d}/*OPD*00*.fits"))
                objcorr_files = glob(f"{d}/*OBJ_CORR*.fits")
                phot_files = np.sort(glob(f"{d}/*PHOT*00*.fits"))
                spectra_files = np.sort(glob(f"{d}/*SPECTRUM*00*.fits"))
                f.write("Files loaded successfully!\n")
            except FileNotFoundError:
                f.write(f"No files found in {d}\n")
                return
            # TODO: using the tpl information, make subdirs in the output dir

            f.write(objcorr_files[0])
            f.write("\n")

            try:
                # Do this here and then again later when the qcparams['custom'] has been updated
                export_dict_to_df(data_dict, output_dir)
                f.write("Wrote data to dataframes \n")
            except FileNotFoundError:
                f.write("Unable to write dataframes\n")

            if len(data_dict["inst"]["tpl"]) == 0:
                return

            tpl_start = data_dict["inst"]["tpl"]
            formatted_outdir = f"{output_dir}/{tpl_start.replace(':','_')}/"
            f.write(data_dict["inst"]["targname"])

            # make the directories
            try:
                mkdir(formatted_outdir)
            except FileExistsError:
                f.write(
                    f"subdir {formatted_outdir} already exists. Clearing existing plots... \n"
                )
                cleanup_plots(formatted_outdir)
            except FileNotFoundError:
                f.write(f"subdir {formatted_outdir} is not a valid path. Exiting... \n")
                exit()

            # insert_observation(
            #     data_dict["inst"]["targname"],
            #     tpl_start,
            #     data_dict["inst"]["tau0"] * 1000,
            #     db_name,
            # )

            # start with the opds so that a mask can be generated
            if True:
                f.write("Loading OPDs...\n")
                opd_dict = load_opd(opd_files, verbose=verbose)

                f.write("Plotting OPDs...\n")
                mask = plot_opd(
                    data_dict,
                    opd_dict,
                    verbose=verbose,
                    output_dir=formatted_outdir,
                    opd_cutoff=opd_cutoff,
                )
            else:
                f.write("Error processing OPDs")

            # plot the visibilities
            try:
                f.write("Plotting visibilities ... ")
                plot_vis(
                    data_dict,
                    save_fig=save_fig,
                    output_dir=formatted_outdir,
                    verbose=verbose,
                )
                f.write("Visibilities plotted successfully!\n")
            except UnboundLocalError:
                f.write(f"Something wrong with {d} when plotting visibilities \n")

            try:
                f.write("Plotting closure phases ...")
                plot_cphase(
                    data_dict,
                    output_dir=formatted_outdir,
                    save_fig=save_fig,
                    verbose=verbose,
                )
                f.write("Closure phases plotted successfully!\n")
            except UnboundLocalError:
                f.write(f"Something wrong with {d} when plotting closure phases\n")

            if not skip_plotting:
                try:
                    f.write("Plotting fringe peaks ... ")
                    mywl = 3.5
                    if band == "N":
                        mywl = 8.5
                    do_obj_corr_plots(
                        objcorr_files,
                        band,
                        tplstart=tpl_start,
                        output_dir=formatted_outdir,
                        verbose=verbose,
                        wl=mywl,
                        save_fig=save_fig,
                    )
                    f.write("Plotting fringe peaks successful \n ")
                except Exception as e:
                    f.write(
                        "Something went wrong while plotting group delay or fringe peaks...\n"
                    )
                    print(e)

                sof = ""
                try:
                    sof = find_sof(data_dir, tpl_start, band=band)
                    f.write(sof)
                except:
                    f.write("SOF not found, not making waterfalls \n")

                if len(sof) > 0:  # True:
                    do_waterfall(
                        sof,
                        output_dir=formatted_outdir,
                        verbose=verbose,
                        save_fig=save_fig,
                    )
                else:  # except (KeyError, FileNotFoundError) as e:
                    f.write("SOF file missing, skipping waterfall plots\n")

                spectral_dict = load_spectrum(spectra_files, verbose=verbose)
                # f.write(spectral_dict, "test")
                try:
                    plot_spectra(
                        spectral_dict,
                        band,
                        verbose=verbose,
                        save_fig=save_fig,
                        output_dir=formatted_outdir,
                    )
                except (KeyError, ValueError) as e:
                    f.write(
                        f"Something went wrong while creating the spectral plots... {e}\n"
                    )

                # phot_dict = load_phot_beams(phot_files, verbose=verbose)
                # then plot the photometries
                #

            pdfs = np.sort(glob(formatted_outdir + "/*.pdf"))
            qc1 = np.mean(
                np.array(data_dict["qcparams"]["custom"]["vis2"]["mean"])
                / np.array(data_dict["qcparams"]["custom"]["vis2"]["std"])
            )
            qc2 = np.mean(
                np.array(data_dict["qcparams"]["custom"]["cphase"]["mean"])
                / np.array(data_dict["qcparams"]["custom"]["cphase"]["std"])
            )
            qc3 = 0.0
            for j in range(1, 7):
                qc3 += data_dict["qcparams"]["eso"][f"DPHASE{j} STDEV"]
            qc3 /= 6
            create_title_page(
                output_path=f"{formatted_outdir}/_title_page.pdf",
                title=f"Target: {data_dict['inst']['targname']}",
                subtitle=f"TPL Start: {tpl_start}",
                tau0=f"Average tau0: {data_dict['inst']['tau0']*1000 :0.2f} ms",
                qc1=f"Vis2 SNR @{data_dict['qcparams']['custom']['vis2']['wavelength'][0]*1e6:.2f} um: {qc1:.2f}",
                qc2=f"t3phi SNR @{data_dict['qcparams']['custom']['cphase']['wavelength'][0]*1e6:.2f} um: {qc2:.2f}",
                qc3=f"Mean DPHASE STD: {qc3:.2f}",
                flux=f"Catalog Flux={data_dict['inst']['cat_flux']:.1f}",
                diameter=f"Catalog Diam.={data_dict['inst']['cat_diam']:.1f}",
            )
            pdfs = [x for x in pdfs if "_title_page" not in x]
            pdfs = np.append([f"{formatted_outdir}/_title_page.pdf"], pdfs)
            merge_pdfs(
                pdfs, formatted_outdir + f"/{tpl_start.replace(':','_')}_output.pdf"
            )
            try:
                # Do this here and then again later when the qcparams['custom'] has been updated
                export_dict_to_df(data_dict, output_dir)
            except FileNotFoundError:
                f.write("Unable to write dataframes")

        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                _exit(1)
            else:
                print(e)
                return


if __name__ == "__main__":
    script, configfile = argv
    with open(configfile, "r") as f:
        cf = json.load(f)
    data_dir = cf["data_dir"]
    # output_dir = cf["output_dir"]
    verbose = int(cf["verbose"])  # > 0
    opd_cutoff = int(cf["opd_cutoff"])
    skip_plotting = cf["tables_only"]
    nb_cores = cf["nb_cores"]
    save_fig = True
    try:
        output_dir = cf["output_dir"]
        inout = False
    except KeyError:
        output_dir = data_dir
        inout = True

    if type(data_dir) is list:
        print("list input")
        data_dirs = cf["data_dir"]

        for data_dir in data_dirs:
            if inout:
                output_dir = data_dir

            # db_name = f"{output_dir}/all_obs.db"
            # Initialize the database
            # create_database(db_name)

            band = "LM"
            directories = glob(data_dir + "/*raw_estimates*.rb")
            if ".rb" in data_dir:
                # it's already a directory!
                directories = [data_dir]
            x = []
            for direc in directories:
                x.append(
                    [direc, data_dir, output_dir, verbose, opd_cutoff, skip_plotting]
                )
            print(directories, data_dir)
            with ProcessPoolExecutor(max_workers=nb_cores) as executor:
                executor.map(process_directory, x)

    # print("Now showing all targets that have been processed!")
    # get_obs(db_name)
    # db_name = f"{output_dir}/all_obs.db"
    # Initialize the database
    # create_database(db_name)
    else:
        band = "LM"
        print("here?")
        directories = glob(data_dir + "/*raw_estimates*.rb")
        if ".rb" in data_dir:
            # it's already a directory!
            directories = [data_dir]

        x = []
        for direc in directories:
            x.append([direc, data_dir, output_dir, verbose, opd_cutoff, skip_plotting])
        print(directories, data_dir)
        with ProcessPoolExecutor(max_workers=nb_cores) as executor:
            executor.map(process_directory, x)

    # print("Now showing all targets that have been processed!")
    # get_obs(db_name)
