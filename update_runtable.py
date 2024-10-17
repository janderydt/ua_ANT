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
       expfolder = os.getcwd()+'/cases/'+data_global['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'/'
       exptable = 'RunTable_'+data_global['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'.csv'
       if os.path.isfile(expfolder+exptable):
          print('Reading data from '+exptable)
          data_exp = read_runinfo(expfolder+exptable,runtype)
          data_global.loc[i]=data_exp.loc[0]
       else:
          print('Cannot find '+exptable+', no changes to global run table.')

# save intermediate version of updated global runtable
# save_runinfo(data_global, runtable_global+".tmp")

# now check for rows with pgid~=0. These are experiments that did not finish cleanly before the end of the walltime.
# we set pgid=0, error=1, restart=0
pd_copy = data_global.copy()
pgid = pd_copy['pgid'].values
for i in range(pgid.shape[0]):
   if pgid[i] != 0:
      print('non-zero pgid found:')
      print(data_global.loc[[i]])
      #correct = input("Would you like to set pgid=0, error=1, restart=0 in the Runtable? yes=1, no=0: ")
      correct=1 #int(correct)
      if correct == 1:
         data_global.at[i,'pgid']=0
         data_global.at[i,'Error']=1
         data_global.at[i,'Restart']=0
         data_global.at[i,'Running']=0
         data_global.at[i,'Submitted']=0
         data_global.at[i,'Comments']='walltime exceeded'
# save modified version of the global runtable
save_runinfo(data_global, runtable_global)
