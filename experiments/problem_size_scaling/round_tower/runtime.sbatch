#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:05:00
#SBATCH -J D2-s-rt
#SBATCH -o /nobackup/wxmv95/problem_size_scaling/round_tower/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/problem_size_scaling/round_tower/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-63

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 9
@ tower_height = ($SLURM_ARRAY_TASK_ID / 8) * 8 + 8
@ number_of_towers = ($SLURM_ARRAY_TASK_ID % 8) * 2 + 2

echo Rows in towers: $tower_height
echo Number of towers: $number_of_towers

time ../../../build/Delta2 -t 0.1 -s 13 --print_priority 1 -n 1 --local_ts --scenario_size $tower_height --scenario_size2 $number_of_towers --threads 64 > /nobackup/wxmv95/problem_size_scaling/round_tower/results/round_tower_runtime_$number_of_towers.$tower_height.$SLURM_JOB_ID.results
