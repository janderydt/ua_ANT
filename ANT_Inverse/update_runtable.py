import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import read_runinfo, save_runinfo

# read input from config file
config_file = sys.argv[1]

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
      idrange = ln.strip("idrange=").strip("\n").strip("[").strip("]").split(":")

# read global runtable
data_global = read_runinfo(runtable_global,runtype)

# update experiment info within required range
for i in range(data_global.shape[0]):
   ExpID = data_global['ExpID'].values[i]
   if float(idrange[0]) <= ExpID <=float(idrange[1]):
       expfolder = os.getcwd()+'/cases/'+table['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'/'
       exptable = 'RunTable_'+data_global['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'.csv'
       if os.path.isfile(expfolder+exptable):
          print('Reading data from '+exptable)
          data_exp = read_runinfo(expfolder+exptable,runtype)
          data_global[i]=data_exp[0]
          print(runtable)
       else:
          print('Cannot find '+exptable+', no changes to global run table.')

# save updated global runtable
#save_runinfo(data_global, runtable_global)
