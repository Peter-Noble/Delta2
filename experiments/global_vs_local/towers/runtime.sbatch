#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J D2-g-t
#SBATCH -o /nobackup/wxmv95/global_vs_local/towers/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/global_vs_local/towers/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-13

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set towers = ( 1 2 4 8 16 32 64 )
set methods = ( "local_ts" "global_ts" )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 0
@ tower_index = $SLURM_ARRAY_TASK_ID % 7
@ towers = `python -c "print('$towers'.split()[$tower_index])"`
@ index = $SLURM_ARRAY_TASK_ID / 7
set method = `python -c "print('$methods'.split()[$index])"`
set method_flag = `python -c "print('--local_ts' if int($index) == 0 else '')"`
@ number_of_threads = 64

echo Towers: $towers
echo Method: $method

time ../../../build/Delta2 -s 1 --print_priority 1 -f 5 $method_flag --scenario_size $towers --threads $number_of_threads > /nobackup/wxmv95/global_vs_local/towers/results/towers_runtime_$towers.t$number_of_threads.m$method.$SLURM_JOB_ID.results
