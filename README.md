A wrapper to launch any number of &Uacute;a jobs on a Linux machine or the ARCHER2 hpc cluster. The code was originally developed for a circum-Antarctic &Uacute;a configuration as part of [CONNECT](https://gtr.ukri.org/projects?ref=MR%2FW011816%2F1) and [OCEAN:ICE](https://ocean-ice.eu/), but can be used for any &Uacute;a domain.

The user can currently chose to run Inverse or Diagnostic simulations; Transient simulations are in development. In each case, a list of experiments needs to be provided in a csv file (RunTable.csv) in the corresponding ANT_Inverse or ANT_Diagnostic folder. The csv file contains information about model parameters and geometry for each experiment, which is passed on to &Uacute;a. An example csv file can be found in the ANT_Inverse and ANT_Diagnostic folders.

## Linux
Launch new jobs from the top directory (/ANT) using `python ua_jobs_master.py`. 

## ARCHER2
Launch new jobs from within the ANT_Inverse, ... directory using `sbatch --export=ALL -A budget-code ./run_ua_archer2.sh`
