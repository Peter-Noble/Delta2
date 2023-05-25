#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J Delta2
#SBATCH -o /nobackup/wxmv95/thread_scaling/pairs/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/thread_scaling/pairs/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-65

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set pairs = ( 64 128 512 2048 8192 16384 )
set threads = ( 1 2 4 5 16 17 32 33 64 65 128)

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 0
@ pair_index = $SLURM_ARRAY_TASK_ID % 6
@ pairs = `python -c "print('$pairs'.split()[$pair_index])"`
@ index = $SLURM_ARRAY_TASK_ID / 6
@ number_of_threads = `python -c "print('$threads'.split()[$index])"`

echo Pairs: $pairs
echo Threads: $number_of_threads

time ../../../build/Delta2 -s 3 --print_priority 1 -f 5 --local_ts --scenario_size $pairs --threads $number_of_threads > /nobackup/wxmv95/thread_scaling/pairs/results/pairs_runtime_$pairs.t$number_of_threads.$SLURM_JOB_ID.results
