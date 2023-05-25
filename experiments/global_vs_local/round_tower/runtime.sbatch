#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J D2-g-rt
#SBATCH -o /nobackup/wxmv95/global_vs_local/round_tower/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/global_vs_local/round_tower/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-15

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set heights = ( 8 16 24 32 40 48 56 64 )
set methods = ( "local_ts" "global_ts" )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 0
# @ tower_height = ($SLURM_ARRAY_TASK_ID % 8) * 8 + 8
@ tower_index = $SLURM_ARRAY_TASK_ID % 8
@ tower_height = `python -c "print('$heights'.split()[$tower_index])"`
# @ number_of_towers = ($SLURM_ARRAY_TASK_ID % 8) * 2 + 2
@ number_of_towers = 50
@ index = $SLURM_ARRAY_TASK_ID / 8
set method = `python -c "print('$methods'.split()[$index])"`
set method_flag = `python -c "print('--local_ts' if int($index) == 0 else '')"`
@ number_of_threads = 64

echo Rows in towers: $tower_height
echo Method: $method

time ../../../build/Delta2 -t 0.1 -s 13 --print_priority 1 -n 1 $method_flag --scenario_size $tower_height --scenario_size2 $number_of_towers --threads $number_of_threads > /nobackup/wxmv95/global_vs_local/round_tower/results/round_tower_runtime_$number_of_towers.$tower_height.t$number_of_threads.m$method.$SLURM_JOB_ID.results
