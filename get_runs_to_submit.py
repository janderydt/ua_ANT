import os
import sys
import numpy as np
import random

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

with open(config_file, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("idrange="):
      idrange = ln.strip("idrange=[").strip("]\n")
      idrange = idrange.split(':')
      idrange = range(int(idrange[0]),int(idrange[1]))

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
    # THE INDEX HAS BEEN ADVANCED BY 1 TO COMPLY WITH MATLAB SYNTAX
    Inew = np.where(expid == 0)[0] + 1 
    #nInew = len(Inew)
    Inew_str = " ".join(str(x) for x in Inew)
    print(Inew_str)

if flag == "existing":
    # find existing runs that do not have any errors, are not submitted, running or finished
    # THE INDEX HAS BEEN ADVANCED BY 1 TO COMPLY WITH MATLAB SYNTAX
    Iexisting = np.where((expid != 0) & (error == 0) & \
            (submitted == 0) & (running == 0) & (finished == 0))[0] + 1
    #nIexisting = len(Iexisting)
    Iexisting_str = " ".join(str(x) for x in Iexisting)
    print(Iexisting_str)

if flag == "all":
    # find all runs to submit
    # THE INDEX HAS BEEN ADVANCED BY 1 TO COMPLY WITH MATLAB SYNTAX
    Iall = np.where((error == 0) & (submitted == 0) & (running == 0) & (finished == 0))[0] + 1
    #nIexisting = len(Iexisting)
    Iall_str = " ".join(str(x) for x in Iall)
    print(Iall_str)

if flag == "count":  
    # find all runs to submit
    Iall = np.where((error == 0) & (submitted == 0) & (running == 0) & (finished == 0))[0]
    #nIexisting = len(Iexisting)
    count = len(Iall)
    print(count)

if flag == "expid":
    # find all runs to submit
    Iall = np.where((error == 0) & (submitted == 0) & (running == 0) & (finished == 0))[0]
    # find corresponding experiment ids
    ExpID_tmp = expid[Iall.astype(int)]
    # find expid's that are zero and assign unique id using a random number 
    # generator and excluding any existing ids
    I_zero = np.where(ExpID_tmp == 0)[0].astype(int)
    ExpID_to_exclude = set(expid)

    for i in I_zero:
        newid = random.choice(list(set([x for x in idrange]) - ExpID_to_exclude))
        ExpID_tmp[i] = newid
        expid=np.append(expid,newid)
        ExpID_to_exclude = set(expid)

    ExpID_str = " ".join(str(x) for x in ExpID_tmp)
    print(ExpID_str)

