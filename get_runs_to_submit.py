import os
import sys
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import read_runinfo, save_runinfo

# read input from config file
config_file = sys.argv[1]
flag = sys.argv[2]

with open(config_file, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtype="):
      runtype = ln.strip("runtype=").strip("\n")

with open(config_file, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtable="):
      runtable_global = ln.strip("runtable=").strip("\n")

# Initialize the Run table, which is used to store details about each simulation
runtable = read_runinfo(runtable_global,runtype)

runs_to_submit = 0

expid = runtable['ExpID'].values
error = runtable['Error'].values
submitted = runtable['Submitted'].values
running = runtable['Running'].values
finished = runtable['Finished'].values

if flag == "new":
    # find runs without experiment id and shift table index by 1 to comply with Matlab syntax
    Inew = np.where(expid == 0)[0] + 1
    #nInew = len(Inew)
    Inew_str = " ".join(str(x) for x in Inew)
    print(Inew_str)

if flag == "existing":
    # find existing runs that do not have any errors, are not submitted, running or finished
    Iexisting = np.where((expid != 0) & (error == 0) & \
            (submitted == 0) & (running == 0) & (finished == 0))[0] + 1
    #nIexisting = len(Iexisting)
    Iexisting_str = " ".join(str(x) for x in Iexisting)
    print(Iexisting_str)
