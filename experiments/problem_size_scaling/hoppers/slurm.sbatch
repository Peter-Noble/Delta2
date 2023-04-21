#!/bin/bash
#SBATCH --job-name=description
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:10:00
#SBATCH -J Delta2
#SBATCH -o stdout/stdout.%J.log
#SBATCH -e stdout/stderr.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide
#SBATCH --array=1-4

source /etc/profile.d/modules.sh

source ../../../hamilton8.sh
source ../../../build/set_env.sh
   
##./fv-0.05 --threads $SLURM_ARRAY_TASK_ID --threading-model native > ccz4-native-$SLURM_ARRAY_TASK_ID.results

time ../../../build/Delta2 -t 0.1 -s 15 --print_priority 1 -f 1 --local_ts --scenario_size $SLURM_ARRAY_TASK_ID > round_tower_$SLURM_ARRAY_TASK_ID.results
