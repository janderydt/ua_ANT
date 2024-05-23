#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ANT_MultiSerial
#SBATCH --time=12:00:0
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block

#SBATCH --account=n02-MRW011816
#SBATCH --partition=standard
#SBATCH --qos=standard

# Make MCR available
MCR=$WORK/MCR_2023b/R2023b/

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

# Loop over the nodes assigned to the job
for nodeid in $nodelist
do
    # Loop over 32 subjobs on each node pinning each to different cores
    for i in $(seq 1 32)
    do
        # Launch subjob overriding job settings as required and in the background
        # Make sure to change the amount specified by the `--mem=` flag to the amount
        # of memory required. The amount of memory is given in MiB by default but other
        # units can be specified. If you do not know how much memory to specify, we
        # recommend that you specify `--mem=1500M` (1,500 MiB).
        srun --nodelist=${nodeid} --nodes=1 --ntasks=1 --ntasks-per-node=1 \
        --exact --mem=8000M --output /dev/null \
        --error stderr${nodeid}_${i}.out ./Ua_MCR.sh $MCR &
        # wait 10 min to make sure first job has started, then 1 min between successive jobs
        if [ $i == 1 ]; then
        sleep 900
        else
        sleep 60
        fi
    done
done

# Wait for all subjobs to finish
wait

