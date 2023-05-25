#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J D2-g-w
#SBATCH -o /nobackup/wxmv95/global_vs_local/waterfall/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/global_vs_local/waterfall/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-11

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

# Multiples of 5 to match the number of particles used in the hopper experiment
set width = ( 5 10 20 40 80 160 )
set methods = ( "local_ts" "global_ts" )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 5
@ width_index = $SLURM_ARRAY_TASK_ID % 6
@ width = `python -c "print('$width'.split()[$width_index])"`
@ index = $SLURM_ARRAY_TASK_ID / 6
set method = `python -c "print('$methods'.split()[$index])"`
set method_flag = `python -c "print('--local_ts' if int($index) == 0 else '')"`
@ number_of_threads = 64

echo Width: $width
echo Method: $method

time ../../../build/Delta2 -t 0.1 -s 5 --print_priority 1 -f 1 $method_flag --scenario_size $width --threads $number_of_threads > /nobackup/wxmv95/global_vs_local/waterfall/results/waterfall_runtime_$width.t$number_of_threads.m$method.$SLURM_JOB_ID.results
