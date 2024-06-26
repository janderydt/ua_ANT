#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ANT_uaconfig
#SBATCH --time=24:00:00
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=30
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

#SBATCH --partition=standard
#SBATCH --qos=standard

#######################################################################################################
# Run multiple instances of UA
# sbatch --export=ALL,UA_CONFIG=<path to UA config file>,ACC=n02-xxxxx -A n02-xxxxx ./run_ua_archer2.sh 
#######################################################################################################

# Add top directory to python path (this makes the utils.py file visible to this script)
CASEDIR=$WORK/ua/cases/ANT
export PYTHONPATH=$PWD:$CASEDIR:$PYTHONPATH

# Make MCR available
MCR=$WORK/MCR_2023b/R2023b/

# Shorter variable name
JOBID=$SLURM_JOB_ID

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
currenttime=$(date +"%d-%m-%Y %H:%M:%S")
jobname=$(sacct -j ${SLURM_JOB_ID}  --format=Jobname%15 2>&1 | sed -n 3p) 
echo ------------------------------------------------------------------------------ >> jobs_master_ARCHER2.log
echo "$currenttime || STARTING $jobname (Config file $UA_CONFIG, JobID $SLURM_JOB_ID)" >> jobs_master_ARCHER2.log
echo ------------------------------------------------------------------------------ >> jobs_master_ARCHER2.log

# start timer
timestart=$(date +"%s")

# make local copy of runtable
python copy_runtable.py $UA_CONFIG

# Loop over the nodes assigned to the job
for nodeid in $nodelist
do
    # Loop over 30 subjobs on each node pinning each to different cores
    for i in $(seq 1 30)
    do
        # update remaining walltime in config file
	timenow=$(date +"%s")
	seconds_expired=$(expr $timenow - $timestart)
	python update_walltime.py $UA_CONFIG ${seconds_expired}

        # Launch subjob overriding job settings as required and in the background
        # Make sure to change the amount specified by the `--mem=` flag to the amount
        # of memory required. The amount of memory is given in MiB by default but other
        # units can be specified. If you do not know how much memory to specify, we
        # recommend that you specify `--mem=1500M` (1,500 MiB).
        srun --nodelist=${nodeid} --nodes=1 --ntasks=1 --ntasks-per-node=1 \
        --exact --mem-per-cpu=1500M --output /dev/null \
        --error stderr_jobid${JOBID}_node${nodeid}_job${i}.out ./Ua_MCR.sh $MCR $UA_CONFIG & 

        # pause until ua job has been submitted
        submitted=0
        while [ $submitted -eq 0 ]
	do 
	    if [ -e ${JOBID}_job_submitted ] ; then
                submitted=1
		rm -f ${JOBID}_job_submitted
	        sleep 5 
	    else
 		sleep 1
            fi
        done
    done
done

# Wait for all subjobs to finish
wait

# gather information about job
currenttime=$(date +"%d-%m-%Y %H:%M:%S")
timeelapsed=$(sacct -j ${JOBID}  --format=Elapsed 2>&1 | sed -n 3p) # elapsed time
EJ=$(sacct -j ${JOBID}  --format=ConsumedEnergy 2>&1 | sed -n 3p) # energy usage
exit=$(sacct -j ${JOBID}  --format=Exitcode 2>&1 | sed -n 3p) # exit code
ex
# Write information to jobs_master_ARCHER2.log
echo -------------------------------------------------------------------------------------------------- >> jobs_master_ARCHER2.log
echo "$currenttime || ENDING $jobname (Config file $UA_CONFIG, JobID $SLURM_JOB_ID) || Exit code $exit" >> jobs_master_ARCHER2.log
echo "Time elapsed: $timeelapsed" >> jobs_master_ARCHER2.log
echo "Energy Consumption [J]: $EJ" >> jobs_master_ARCHER2.log
echo -------------------------------------------------------------------------------------------------- >> jobs_master_ARCHER2.log

# Update global RunTable
python update_runtable.py $UA_CONFIG

# Clean up
rm $(JOB_ID)_job_submitted
find . -maxdepth 1 -name 'stderr_jobid${JOBID}*.out' -size 0 | xargs rm -rf

# Relaunch script
THISSCRIPT=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
python -c "from utils import submit_job; submit_job(budget_code='$ACC',sbatch_script='$THISSCRIPT',input_var=['UA_CONFIG=$UA_CONFIG','ACC=$ACC'])"
