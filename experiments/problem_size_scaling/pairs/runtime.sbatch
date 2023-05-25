#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J D2-s-p
#SBATCH -o /nobackup/wxmv95/problem_size_scaling/pairs/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/problem_size_scaling/pairs/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-11

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set pairs = ( 1 2 4 8 16 32 64 128 512 2048 8192 16384 )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 0
@ pairs = `python -c "print('$pairs'.split()[$SLURM_ARRAY_TASK_ID])"`
@ number_of_threads = 64

echo Pairs: $pairs

time ../../../build/Delta2 -s 3 --print_priority 1 -f 5 --local_ts --scenario_size $pairs --threads $number_of_threads > /nobackup/wxmv95/problem_size_scaling/pairs/results/pairs_runtime_$pairs.t$number_of_threads.$SLURM_JOB_ID.results
