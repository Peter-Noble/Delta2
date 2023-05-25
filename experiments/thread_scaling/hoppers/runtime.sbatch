#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J Delta2
#SBATCH -o /nobackup/wxmv95/thread_scaling/hoppers/stdout/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/thread_scaling/hoppers/stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=0-71

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

set hoppers = ( 1 2 4 8 16 32 )
set threads = ( 1 2 4 8 16 32 48 64 80 96 112 128 )

# @ SLURM_JOB_ID = 0
# @ SLURM_ARRAY_TASK_ID = 5
@ hoppers_index = $SLURM_ARRAY_TASK_ID % 6
@ hoppers = `python -c "print('$hoppers'.split()[$hoppers_index])"`
@ index = $SLURM_ARRAY_TASK_ID / 6
@ number_of_threads = `python -c "print('$threads'.split()[$index])"`

echo Hoppers: $hoppers
echo Threads: $number_of_threads

time ../../../build/Delta2 -t 0.1 -s 15 --print_priority 1 -f 1 --local_ts --scenario_size $hoppers --threads $number_of_threads > /nobackup/wxmv95/thread_scaling/hoppers/results/hoppers_runtime_$hoppers.t$number_of_threads.$SLURM_JOB_ID.results
