A wrapper to launch any number of &Uacute;a jobs on a Linux machine or ARCHER2 hpc cluster. The code was originally developed for a circum-Antarctic &Uacute;a configuration as part of [CONNECT](https://gtr.ukri.org/projects?ref=MR%2FW011816%2F1) and [OCEAN:ICE](https://ocean-ice.eu/), but can be used for any &Uacute;a domain.

The user can currently choose to run Inverse or Diagnostic simulations; Transient simulations are in development. In each case, a list of experiments needs to be provided in a csv file (RunTable.csv) in the corresponding ANT_Inverse or ANT_Diagnostic folder. The csv file contains information about model parameters and geometry for each experiment, which is passed on to &Uacute;a. Example csv files can be found in the ANT_Inverse and ANT_Diagnostic folders.

## Linux desktop
Launch new jobs from the top directory (/ANT) using `python ua_jobs_master.py`. Jobs will run in the background so users can close the terminal and logout.
The user will be asked whether they want to launch inverse or diagnostic simulations, and how many. The code will read the relevant RunTable.csv file, and submit the chosen number of jobs from that table. Jobs that are already running or have finished will be ignored.

## ARCHER2
(Currently only available for inverse simulations and diagnostic simulations)
Launch new jobs from the ANT_Inverse or ANT_Diagnostic folder using `sbatch --export=ALL,UA_CONFIG=<path to UA config file>,ACC=budget-code -A budget-code -N number_of_nodes ./run_ua_archer2.sh`.
The user needs to provide a configuration file (ua_config.txt) with basic settings such as the type of experiment and requested walltime.