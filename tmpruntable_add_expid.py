import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import read_runinfo, save_runinfo

# read input from config file
config_file = sys.argv[1]
row_number = int(sys.argv[2])
expid = int(sys.argv[3])

with open(config_file, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtype="):
      runtype = ln.strip("runtype=").strip("\n")

with open(config_file, "r") as fi: # Open the file in read mode
  for ln in fi:
    if ln.startswith("runtable="):
      runtable_global = ln.strip("runtable=").strip("\n")
    
# read temporary global runtable
data_global = read_runinfo(runtable_global+".tmp",runtype)

# update experiment info within required range
data_global.at[row_number,'ExpID']=expid

# save updated global runtable to temporary file
save_runinfo(data_global, runtable_global+".tmp")


