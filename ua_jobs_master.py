from utils import submit_batch, active_jobs, read_runinfo, save_runinfo
import datetime
import logging
import pandas as pd
import sys
import os

# Top-level function that controls Ua job submissions. This is written for Linux
# and some elements will need to be modified for use with Windows/MacOS
if __name__ == "__main__":
   
    # First check what type of run needs to be submitted.
    inputselected = 0

    while inputselected == 0:
        
        runtype = input("Choose runtype (1. Inverse, 2. Diagnostic): ")
        runtype = int(runtype)

        if runtype == 1:
            runtype = 'Inverse'
            inputselected = 1
            runinfo = '> CHECKING INVERSE RUNS'
            os.chdir('./ANT_Inverse')
        elif runtype == 2:
            runtype = 'Diagnostic'
            inputselected = 1
            runinfo = '> CHECKING DIAGNOSTIC RUNS'
            os.chdir('./ANT_Diagnostic')
        else:
            print('Runtype does not exist. Choose 1 or 2.')

    # Initialize the log file
    logging.basicConfig(filename=os.getcwd()+'/jobs_master.log', level=logging.INFO, encoding='utf-8', format='%(message)s')

    logging.info(runinfo)

    logging.info('============================')
    logging.info(str(datetime.datetime.now()))
    logging.info('============================')
    logging.info('> Reading table with run info')

    # Initialize the Run table, which is used to store details about each simulation
    runtable = read_runinfo('RunTable.csv',runtype)

    # Check that processes that are marked as running are indeed running
    logging.info('> Check that processes that are marked as running are indeed running')

    run_counter = 0

    for i in range(runtable.shape[0]):
        if runtable['Running'].values[i] == 1:

            pgid = runtable['pgid'].values[i]
            command = "ps -g " + str(pgid) + " -o pgid,comm"
            out = active_jobs(command)

            if str(pgid) in out:
                logging.info('   ...Process '+str(pgid)+' '+'(ExpID '+str(runtable['ExpID'][i])+') is running')
                run_counter +=1

            else:
                logging.info('   ...Table indicated process is running, but could not find the pgid')
                runtable.loc[i,'Running'] = 0
                runtable.loc[i,'Submitted'] = 0
                runtable.loc[i,'Error'] = 1
                runtable.loc[i,'ErrorTime'] = datetime.datetime.now()
                runtable.loc[i,'pgid'] = 0 
                runtable.loc[i,'Comments'] = runtable.loc[i,'Comments'] + "; Table indicated process is running, but could not find the pgid"
                save_runinfo(runtable,'RunTable.csv')

    logging.info('   ...'+str(run_counter)+' jobs running')

    # Launching new job if capacity allows
    logging.info('> Launching new job if capacity allows')

    if run_counter<=5:
    
    	# first check RunTable to see if there are any jobs that are not running,
        # but are not finished yet and did not throw any errors
        ind=runtable.index[(runtable['Running']==0) & (runtable['Error']==0) & (runtable['Finished']==0)]

        if len(ind)>0:
        
            options = ""
            command = "nohup ./submit_nohup.sh &"
            pgid = submit_batch(options, command)
            
        else:
            
        # get table with new runs and append to RunTable, then launch job
            newruns = read_runinfo('NewRuns.csv',runinfo)
        
            if newruns.shape[0]>0:

                runtable = pd.concat([runtable, newruns.iloc[[0]]], ignore_index=True, sort=False)
                newruns = newruns.drop(labels=0)

                save_runinfo(runtable, 'RunTable.csv')
                save_runinfo(newruns, 'NewRuns.csv')

                options = ""
                command = "nohup ./submit_nohup.sh &"
                pgid = submit_batch(options, command)

            else:
            
                logging.info('   ...No new runs in NewRuns.scv: Nothing to do.')

