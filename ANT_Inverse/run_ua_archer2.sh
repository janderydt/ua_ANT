#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=ANT_MultiSerial
#SBATCH --time=0:10:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
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

# Loop over 128 subjobs pinning each to a different core
for i in $(seq 1 2)
do
# Launch subjob overriding job settings as required and in the background
# Make sure to change the amount specified by the `--mem=` flag to the amount
# of memory required. The amount of memory is given in MiB by default but other
# units can be specified. If you do not know how much memory to specify, we
# recommend that you specify `--mem=1500M` (1,500 MiB).
srun --nodes=1 --ntasks=1 --ntasks-per-node=1 \
      --exact --mem=5000M ./Ua_MCR.sh $MCR 1>>matlab_std${i}.out 2>>matlab_err${i}.out &
# wait 3 min to make sure previous job has started
sleep 180
done

# Wait for all subjobs to finish
wait
