#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J Delta2
#SBATCH -o /nobackup/wxmv95/thread_scaling/hoppers/stdout.%J.log
#SBATCH -e /nobackup/wxmv95/thread_scaling/hoppers/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
# #SBATCH --array=32

source /etc/profile.d/modules.sh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

# setenv EXPORT_DIR /nobackup/wxmv95/thread_scaling/hoppers/export

source /apps/developers/tools/vtune/2021.9/1/default/vtune/latest/vtune-vars.sh

@ size = 4
@ threads = 4

vtune -collect hotspots -r /nobackup/wxmv95/thread_scaling/hoppers/size$size.f1.t$threads.$SLURM_JOB_ID.hs -finalization-mode=fast ../../../build/Delta2 -t 0.5 -s 15 --print_priority 1 -f 1 --local_ts --scenario_size $size --threads $threads > /nobackup/wxmv95/thread_scaling/hoppers/results/size$size.f5.t$threads.$SLURM_JOB_ID.results
# vtune -collect hotspots -r /nobackup/wxmv95/problem_size_scaling/round_tower/r$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.hs ../../../build/Delta2 -t 0.1 -s 13 --print_priority 1 -n 1 --local_ts --scenario_size $SLURM_ARRAY_TASK_ID --scenario_size2 2 > /nobackup/wxmv95/problem_size_scaling/round_tower/results/round_tower_$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.results
# vtune -collect hotspots -r r$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.hs -finalization-mode=deferred ../../../build/Delta2a -t 0.1 -s 13 --print_priority 1 -n 1 --local_ts --scenario_size $SLURM_ARRAY_TASK_ID > results/round_tower_$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.results
# vtune -collect threading -r r$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.th ../../../build/Delta2 -t 0.1 -s 13 --print_priority 1 -n 10 --local_ts --scenario_size $SLURM_ARRAY_TASK_ID > results/round_tower_$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.results
