from utils import submit_batch, active_jobs, read_runinfo, save_runinfo
import datetime
import logging
import pandas as pd
import sys

# Top-level function
if __name__ == "__main__":

    # Create the file
    # and output every level since 'DEBUG' is used
    # and remove all headers in the output
    # using empty format=''
    logging.basicConfig(filename='jobs_master.log', level=logging.INFO, encoding='utf-8', format='%(message)s')
 
    logging.info('============================')
    logging.info(str(datetime.datetime.now()))
    logging.info('============================')
    logging.info('> Reading table with run info')

    runtable = read_runinfo('RunTable.csv')

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

    logging.info('> Launching new job if capacity allows')

    if run_counter<=15:
    
    	# first check RunTable to see if any runs require a restart
        ind1=runtable.index[(runtable['Restart']==1) & (runtable['Finished']==0)]
        # check RunTable to see if any runs require a fresh start
        ind2=runtable.index[(runtable['Restart']==0) & (runtable['Submitted']==0) & (runtable['Finished']==0) & (runtable['Error']==0)]
        
        if len(ind1)>0 or len(ind2)>0:
            
            options = ""
            command = "nohup ./submit_nohup.sh &"
            pgid = submit_batch(options, command)
            
        else:
            
        # get table with new runs and append to RunTable, then launch job
            newruns = read_runinfo('NewRuns.csv')
        
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




