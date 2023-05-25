#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J Delta2
#SBATCH -o /nobackup/wxmv95/thread_scaling/towers/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/thread_scaling/towers/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-48

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set towers = ( 1 2 4 8 16 32 64 )
set threads = ( 1 2 4 8 16 32 64 )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 0
@ tower_index = $SLURM_ARRAY_TASK_ID % 7
@ towers = `python -c "print('$towers'.split()[$tower_index])"`
@ index = $SLURM_ARRAY_TASK_ID / 7
@ number_of_threads = `python -c "print('$threads'.split()[$index])"`

echo Towers: $towers
echo Threads: $number_of_threads

time ../../../build/Delta2 -s 1 --print_priority 1 -f 5 --local_ts --scenario_size $towers --threads $number_of_threads > /nobackup/wxmv95/thread_scaling/towers/results/towers_runtime_$towers.t$number_of_threads.$SLURM_JOB_ID.results
