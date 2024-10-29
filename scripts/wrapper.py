"""
Author: Jacob Isbell 
Date: October 1, 2024

Script which contains primary logic for MATISSE QA.
"""

import json 
from glob import glob
#from icecream import ic

from load_files import load_raw_int
from vis2_plots import plot_vis 
from t3_plots import plot_cphase
from sys import argv

if __name__ =="__main__":
    script, configfile = argv
    with open(configfile,'r') as f:
        cf = json.load(f)
    data_dir = cf['data_dir']
    output_dir = cf['output_dir']
    verbose = int(cf['verbose']) > 0

    directories = glob(data_dir + '/*raw_estimates*.rb')
    print(directories,data_dir)
    for d in directories:
        files = glob(f'{d}/*RAW_INT_00*.fits')
        print('\n\n' + '#'*128)
        print(f'Working on {d}')
        try:
            print(f'Loading files ... ')
            data_dict = load_raw_int(files, verbose=verbose)
            print(f'Files loaded successfully!')
        except:
            print(f'No files found in {d}')
            continue
        #TODO: using the tpl information, make subdirs in the output dir 
        try:
            print(f'Plotting visibilities ... ')
            plot_vis(data_dict['vis'], save_fig=True, output_dir=output_dir, verbose=verbose)
            print(f'Visibilities plotted successfully!')
        except UnboundLocalError:
            print(f'Something wrong with {d} when plotting visibilities')
            continue

        try:
            print('Plotting closure phases ...')
            plot_cphase(data_dict['cphase'], output_dir=output_dir, save_fig=True, verbose=verbose)
            print(f'Closure phases plotted successfully!')
        except:
            print(f'Something wrong with {d} when plotting closure phases')
            continue


    
