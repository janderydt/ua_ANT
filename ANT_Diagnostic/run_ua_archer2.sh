#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ANT_Diag
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=30
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --partition=standard
#SBATCH --qos=standard

##############################################################################################################
# Run multiple instances of UA
# sbatch --export=ALL,UA_CONFIG=<path to UA config file>,ACC=n02-xxxxx -N xxx -A n02-xxxxx ./run_ua_archer2.sh 
#
# -N xxx is the required number of nodes
# n02-xxxxx is the budget code
# 
# for runs with a higher memory requirement:
# use SBATCH --partition=highmem, SBATCH --qos=standard and set --mem-per-cpu=3000M in the srun command below
#
# for testing:
# use SBATCH --partition=standard, SBATCH --qos=short
##############################################################################################################

# Add top directory to python path (this makes the utils.py file visible to this script)
CASEDIR=$WORK/ua/cases/ANT
export PYTHONPATH=$PWD:$CASEDIR:$PYTHONPATH

# Make MCR available
MCR=$WORK/MCR_2024b/R2024b/

# Shorter variable name
JOBID=${SLURM_JOB_ID}

# Make sure MCR cache (as defined in Ua_MCR.sh) exists
# If you want the cache in a different location, modify it here AND in Ua_MCR.sh
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

# Start a timer to calculate the total runtime at the end of the job
timestart=$(date +"%s")

# make local copy of runtable
python ../copy_runtable.py $UA_CONFIG "" ".tmp"

# update all experiment runtables based on global runtable
python ../update_runtable.py $UA_CONFIG "global" "local"

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

# We kill any jobsteps but not the batch job 5 minutes before the walltime. This is done
# with the --time flag in the srun command below. Here we calculate the remaining walltime
# and subtract 5 minutes to create the TIME_LIMIT variable, which is then passed on to the srun commands.
# 1. Record time left in job, this will be a string of the format d-hh:mm:ss,
# but zero values are truncated e.g. if d=0, hh=0 then the returned string will be mm:ss
TIME_LIMIT="$(squeue -j $SLURM_JOB_ID -h --Format TimeLeft)"
# 2. Trim white spaces
TIME_LIMIT="$(echo ${TIME_LIMIT} | xargs)"
length=${#TIME_LIMIT}
# 3. Recover full date format (hh:mm:ss)
if [ "$length" -eq "4" ]; then
    TIME_LIMIT="00:0$TIME_LIMIT"
elif [ "$length" -eq "5" ]; then
    TIME_LIMIT="00:$TIME_LIMIT"
elif [ "$length" -eq "7" ]; then
    TIME_LIMIT="0$TIME_LIMIT"
fi
# 4. subtract 5 minutes
TIME_LIMIT_SECS=$(( $(date -d "$TIME_LIMIT" "+%s") - $(date -d "00:05:00" "+%s") ))
# 5. convert back to hh:mm:ss format
h=$(( $TIME_LIMIT_SECS / 3600 ))
m=$(( $(($TIME_LIMIT_SECS - $h * 3600)) / 60 ))
s=$(($TIME_LIMIT_SECS - $h * 3600 - $m * 60))
if [ $h -le 9 ];then h=0$h;fi
if [ $m -le 9 ];then m=0$m;fi
if [ $s -le 9 ];then s=0$s;fi
TIME_LIMIT="$h:$m:$s"

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
                --exact --mem-per-cpu=1500M --output /dev/null --time=${TIME_LIMIT} \
                --error stderr_jobid${JOBID}_expid${ExpID[$jobs_submitted]}.out ./Ua_MCR.sh \
		$MCR $UA_CONFIG ${SLURM_JOB_ID} "" ${RowNb[$jobs_submitted]} ${ExpID[$jobs_submitted]} \
		& pids+=($!) 

		# add expid to temporary copy of global runtable (note that RowNb follows MATLAB syntax: number of first row is 1 \
		# so we offset the indices by 1 to comply with python numbering convetions)
		python ../tmpruntable_add_expid.py $UA_CONFIG $((${RowNb[$jobs_submitted]}-1)) ${ExpID[$jobs_submitted]}

	        # advance the counter
                ((jobs_submitted++))

		# stagger job submissions / not sure if this should be needed but it seems to resolve Matlab MRC error:
                # MCL:ComponentCache
                # Error: Could not access the MATLAB Runtime component cache. Details: fl:filesystem:PathNotFound
	        # Component cache root: '/work/n02/n02/janryd69/mcr_cache/R2024b'
                # Component name: 'Ua'	
		sleep 1
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

    # Sleep for 5secs - this seems to allow some time for energy and memory statistics to be generated for individual jobsteps
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
    EJ=$(sacct --jobs=${JOBID}.0 --format=ConsumedEnergy 2>&1 | sed -n 3p) # energy usage of batch job 
    MaxRSS=$(sacct --jobs=${JOBID}.${i} --format=MaxRSS 2>&1 | sed -n 3p) # max memory usage of batch job 
    echo " > Energy Consumption: ${EJ}, Maximum memory usage: ${MaxRSS}" >> jobs_master_ARCHER2.log
    for i in $(seq 0 $(( ${#rets[*]}-1 ))); do
        EJ=$(sacct --jobs=${JOBID}.${i} --format=ConsumedEnergy 2>&1 | sed -n 3p) # energy usage for jobstep
        MaxRSS=$(sacct --jobs=${JOBID}.${i} --format=MaxRSS 2>&1 | sed -n 3p) # max memory usage for jobstep 
	echo "      ExpID ${ExpID[$i]}        Energy usage [J]: ${EJ} || Maximum memory usage [bytes]: ${MaxRSS}" >> jobsteps_master_ARCHER2.log
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
    python ../update_runtable.py $UA_CONFIG "local" "global"

    # Clean up directory
    find . -maxdepth 1 -name "stderr_jobid${JOBID}*.out" -size 0 | xargs rm -rf

    # Relaunch script if all jobs finished ok
    if [ $exitflag -eq 0 ]; then
        # calculate number of required nodes for resubmission
        Nb_experiments_to_start=`python ../get_runs_to_submit.py $UA_CONFIG count`
        Nodes_required=$(( ($Nb_experiments_to_start/30) + ($Nb_experiments_to_start%30>0) )) 
        if [ ${Nodes_required} -gt 9 ]; then
	    Nodes_required=9
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
