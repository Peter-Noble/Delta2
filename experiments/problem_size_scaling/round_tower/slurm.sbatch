#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:30:00
#SBATCH -J Delta2
#SBATCH -o stdout/stdout.%J.log
#SBATCH -e stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=32

source /etc/profile.d/modules.sh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

source /apps/developers/tools/vtune/2021.9/1/default/vtune/latest/vtune-vars.sh
   
vtune -collect hotspots -r r$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.hs -finalization-mode=deferred ../../../build/Delta2 -t 0.1 -s 13 --print_priority 1 -n 10 --local_ts --scenario_size $SLURM_ARRAY_TASK_ID > results/round_tower_$SLURM_ARRAY_TASK_ID.$SLURM_JOB_ID.results
