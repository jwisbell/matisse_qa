# matisse_qa
MATISSE Quality Assessment Scripts

These scripts analyze the final and intermediate products of the standard MATISSE data reduction pipeline. By plotting and measuring the squared visibilities, closure phases, and other observables, a summary pdf is made for each observed object giving an overview of the quality of the observation. 

# Installation
After cloning the repository, run 
```
pip install -r requirements.txt
```
to ensure that all necessary packages are installed. 

# Usage
First, modify the config.json file
```
{
  "data_dir": "/path/to/raw/data/",
  "output_dir": "/path/to/destination/for/plots/",
  "verbose": "1",
  "opd_cutoff": 3,
  "tables_only": false,
  "nb_cores": 4
}
```
The `data_dir` should contain one or more *.rb files produced by the MATISSE reduction pipeline. 

`nb_cores` specifies how many cores to use for naive multiprocessing (i.e., running each observation TPL in a thread).

*NOTE:* For faster execution at the expense of not computing the waterfall plots, set `tables_only` to `true`. 

Then execution is straightforward with
```
python wrapper.py config.json
```

All produced plots, logs, and pandas dataframes are saved in the specified `output_dir`


# Outputs 
## Logs 
Each TPL start (individual observation) has an associated log for inspection in the case of issues or missing plots. 

## Plots 
A large number of plots are generated showing
1. Squared visibilities 
2. Visibilities
3. Closure phases
4. Differential Phases
5. Optical Path Delays (OPDs)
6. Photometry (for each telescope)
7. Interferometric "waterfall" plots showing the fringe-tracking information throughout the observation

## Pandas Dataframes
ESO pipeline-created as well as custom quality control metrics for various observables are saved as (pickled) pandas dataframes for further inspection.

