import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from utils import read_runinfo, save_runinfo

# read input from config file
config_file = sys.argv[1]
read_table = str(sys.argv[2]) # either 'global' or 'local'
write_table = str(sys.argv[3]) # either 'global' or "local'

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

# update write_table within required range
for i in range(data_global.shape[0]):
    ExpID = data_global['ExpID'].values[i]
    if float(idrange[0]) <= ExpID <=float(idrange[1]):
        expfolder = os.getcwd()+'/cases/'+data_global['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'/'
        exptable = 'RunTable_'+data_global['Domain'].values[i]+'_'+runtype.strip().strip("\"")+'_'+str(ExpID)+'.csv'
        runtable_local = expfolder+exptable
        if os.path.isfile(runtable_local):
            data_exp = read_runinfo(runtable_local,runtype)
            if read_table == 'global' and write_table == 'local':
                print('Writing data to '+exptable)               
                data_exp.loc[0]=data_global.loc[i]
                # save modified version of the local runtable
                save_runinfo(data_exp, runtable_local)
            elif read_table == 'local' and write_table == 'global':
                print('Writing data to '+runtable_global)               
                data_global.loc[i]=data_exp.loc[0]
                # save modified version of the global runtable
                save_runinfo(data_global, runtable_global)
            else:
                print('Wrong combination of read_table and write_table. They should not be the same.')
        else:
            print('Cannot find '+exptable+', no changes to '+write_table+' run table.')


# save intermediate version of updated global runtable
# save_runinfo(data_global, runtable_global+".tmp")

# if write_table=global: check for rows with pgid~=0. These are experiments that did not finish cleanly before the end of the walltime.
# we set pgid=0, error=1, restart=0
if write_table == 'global':
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