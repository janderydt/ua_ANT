#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ANT_Trans
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=30
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

#SBATCH --partition=highmem
#SBATCH --qos=highmem

##############################################################################################################
# Run multiple instances of UA
# sbatch --export=ALL,UA_CONFIG=<path to UA config file>,ACC=n02-xxxxx -N xxx -A n02-xxxxx ./run_ua_archer2.sh 
#
# -N xxx is the required number of nodes
# -A n02-xxxxx is the budget code
###########################################################################################'##################

# Add top directory to python path (this makes the utils.py file visible to this script)
CASEDIR=$WORK/ua/cases/ANT
export PYTHONPATH=$PWD:$CASEDIR:$PYTHONPATH

# Make MCR available
MCR=$WORK/MCR_2023b/R2023b/

# Shorter variable name
JOBID=${SLURM_JOB_ID}

# Make sure MCR cache (as defined in Ua_MCR.sh) exists
# If you want the cache in a different location, modify it here AND in ua_run/Ua_MCR.sh
if [ ! -d $WORK/mcr_cache ]; then
  mkdir $WORK/mcr_cache
fi

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

# Propagate the cpus-per-task setting from script to srun commands
#    By default, Slurm does not propagate this setting from the sbatch
#    options to srun commands in the job script. If this is not done,
#    process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Get a list of the nodes assigned to this job in a format we can use.
#   scontrol converts the condensed node IDs in the sbatch environment
#   variable into a list of full node IDs that we can use with srun to
#   ensure the subjobs are placed on the correct node. e.g. this converts
#   "nid[001234,002345]" to "nid001234 nid002345"
nodelist=$(scontrol show hostnames $SLURM_JOB_NODELIST)

# Write information to jobs_master_ARCHER2.log
# lock file for writing
while [ -f global_log_active ]; do
    sleep 1
done 
touch global_log_active
currenttime=$(date +"%d-%b-%Y %H:%M:%S")
jobname=$(sacct -j ${JOBID} --format=Jobname%15 2>&1 | sed -n 3p) 
echo "${currenttime} || STARTING ${jobname} (Config file ${UA_CONFIG}, JobID ${JOBID})" >> jobs_master_ARCHER2.log

# start timer
timestart=$(date +"%s")

# make local copy of runtable
python ../copy_runtable.py $UA_CONFIG "" ".tmp"

# get list of experiments to start (read python output as string)
# first obtain the row numbers of experiments that need to be started. NOTE: the index of the first row is 1 (MATLAB syntax)
RowNb_all="`python ../get_runs_to_submit.py $UA_CONFIG all`"
# now obtain the unique ids for each experiment
ExpID_all="`python ../get_runs_to_submit.py $UA_CONFIG expid`"
# convert strings to arrays
RowNb=($RowNb_all)
ExpID=($ExpID_all)
# how many experiments to start?
Nb_experiments_to_start=`python ../get_runs_to_submit.py $UA_CONFIG count`

if [ $Nb_experiments_to_start -gt 0 ]
then
    # define job counter
    jobs_submitted=0
    pids=()
    
    # Loop over the nodes assigned to the job
    for nodeid in $nodelist
    do
        # Loop over 30 subjobs on each node pinning each to different cores
        for i in $(seq 1 30)
        do
            # update remaining walltime in config file
	    timenow=$(date +"%s")
	    seconds_expired=$(expr $timenow - $timestart)
	    python ../update_walltime.py $UA_CONFIG ${seconds_expired}
            
	    # check that job counter <= Nb_experiments_to_start
            if [ $jobs_submitted -lt $Nb_experiments_to_start ]
	    then
                # Launch subjob overriding job settings as required and in the background
                # Make sure to change the amount specified by the `--mem=` flag to the amount
                # of memory required. The amount of memory is given in MiB by default but other
                # units can be specified. If you do not know how much memory to specify, we
                # recommend that you specify `--mem=1500M` (1,500 MiB).
                srun --nodelist=${nodeid} --nodes=1 --ntasks=1 --ntasks-per-node=1 \
                --exact --mem-per-cpu=3000M --output /dev/null \
                --error stderr_jobid${JOBID}_expid${ExpID[$jobs_submitted]}.out ./Ua_MCR.sh \
		$MCR $UA_CONFIG ${SLURM_JOB_ID} "" ${RowNb[$jobs_submitted]} ${ExpID[$jobs_submitted]} \
		& pids+=($!) 

		# add expid to temporary copy of global runtable (note that RowNb follows MATLAB syntax: number of first row is 1 \
		# so we offset the indices by 1 to comply with python numbering convetions)
		python ../tmpruntable_add_expid.py $UA_CONFIG $((${RowNb[$jobs_submitted]}-1)) ${ExpID[$jobs_submitted]}

		# advance the counter
                ((jobs_submitted++))
            fi
        done
    done

    # Write information to jobs_master_ARCHER2.log
    currenttime=$(date +"%d-%b-%Y %H:%M:%S")
    jobname=$(sacct -j ${JOBID} --format=Jobname%15 2>&1 | sed -n 3p)
    echo " > Total number of jobs submitted: ${jobs_submitted} out of ${Nb_experiments_to_start}" >> jobs_master_ARCHER2.log
    echo "---------------------------------------------------" >> jobs_master_ARCHER2.log
    rm global_log_active

    # Wait for all jobsteps to finish and gather exit codes
    rets=()
    for pid in ${pids[*]}; do
        wait $pid
	rets+=($?)
    done

    sleep 5

    # gather information about runtime
    currenttime=$(date +"%d-%b-%Y %H:%M:%S")
    timeelapsed=$(sacct -j ${JOBID}  --format=Elapsed 2>&1 | sed -n 3p) # elapsed time
    exitflag=0

    # Find non-empty error files
    nonempty_error_files=$(find . -maxdepth 1 -name "stderr_jobid${JOBID}*.out" -type f ! -size 0 | wc -l)

    # Write information to jobs_master_ARCHER2.log
    while [ -f global_log_active ]; do
    	sleep 1
    done 
    
    touch global_log_active
       
    echo "${currenttime} || ENDING ${jobname} (Config file ${UA_CONFIG}, JobID ${JOBID})" >> jobs_master_ARCHER2.log
    echo " > Energy Consumption and maximum memory usage:" >> jobs_master_ARCHER2.log
    for i in $(seq 0 $(( ${#rets[*]}-1 ))); do
        EJ=$(sacct --jobs=${JOBID}.${i} --format=ConsumedEnergy 2>&1 | sed -n 3p) # energy usage
        MaxRSS=$(sacct --jobs=${JOBID}.${i} --format=MaxRSS 2>&1 | sed -n 3p) # max memory usage  
	echo "      ExpID ${ExpID[$i]}        Energy usage [J]: ${EJ} || Maximum memory usage [bytes]: ${MaxRSS}" >> jobs_master_ARCHER2.log
    done 
    echo " > Exit codes (if different from zero):" >> jobs_master_ARCHER2.log 
    for i in $(seq 0 $(( ${#rets[*]}-1 ))); do
        if [ ${rets[$i]} -ne 0 ]; then
	    echo "      ExpID ${ExpID[$i]}: exit code ${rets[$i]}" >> jobs_master_ARCHER2.log
	    exitflag=1
	fi
    done
    echo " > Non-empty error log files: ${nonempty_error_files}" >> jobs_master_ARCHER2.log
    if [ ${nonempty_error_files} -gt 0 ]; then
        find . -maxdepth 1 -name "stderr_jobid${JOBID}*.out" -type f ! -size 0 -exec echo "      {}" \; >> jobs_master_ARCHER2.log
        exitflag=1
    fi
    echo " > Time elapsed: ${timeelapsed}" >> jobs_master_ARCHER2.log

    # Copy temporary RunTable
    python ../copy_runtable.py $UA_CONFIG ".tmp" ""

    # Update global RunTable
    python ../update_runtable.py $UA_CONFIG

    # Clean up directory
    find . -maxdepth 1 -name "stderr_jobid${JOBID}*.out" -size 0 | xargs rm -rf

    # Relaunch script if all jobs finished ok
    if [ $exitflag -eq 0 ]; then
        # calculate number of required nodes for resubmission
        Nb_experiments_to_start=`python ../get_runs_to_submit.py $UA_CONFIG count`
        Nodes_required=$(( ($Nb_experiments_to_start/30) + ($Nb_experiments_to_start%30>0) )) 
        if [ ${Nodes_required} -gt 6 ]; then
	    Nodes_required=6
	fi	    
	if [ $Nodes_required -gt 0 ]; then
	    THISSCRIPT=$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/{print $2}')
            currenttime=$(date +"%d-%b-%Y %H:%M:%S")
            echo "${currenttime} || RESUBMITTING ${jobname} on ${Nodes_required} nodes" >> jobs_master_ARCHER2.log 
            python -c "from utils import submit_job; submit_job(budget_code='$ACC',Nnodes='$Nodes_required',sbatch_script='$THISSCRIPT',input_var=['UA_CONFIG=$UA_CONFIG','ACC=$ACC'])"
        else
	    currenttime=$(date +"%d-%b-%Y %H:%M:%S")
            echo "${currenttime} || JOB ${jobname} HAS SUCCESSFULLY FINISHED" >> jobs_master_ARCHER2.log
	fi
    else
        currenttime=$(date +"%d-%b-%Y %H:%M:%S")
        echo "${currenttime} || ERROR in ${jobname} - Aborting" >> jobs_master_ARCHER2.log
    fi
    echo "---------------------------------------------------" >> jobs_master_ARCHER2.log
    
    rm global_log_active
fi
