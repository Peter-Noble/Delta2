#!/bin/csh
#SBATCH --job-name=Delta2
#SBATCH -p multi
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=00:02:00
#SBATCH -J Delta2
#SBATCH -o /nobackup/wxmv95/problem_size_scaling/round_tower/stdout/stdout_maqao.%J.log
#SBATCH -e /nobackup/wxmv95/problem_size_scaling/round_tower/stdout/stderr_maqao.%J.log
# #SBATCH --mail-type=END
# #SBATCH --mail-user=wxmv95@durham.ac.uk
#SBATCH --qos=shortwide

source /etc/profile.d/modules.csh

source ../../../hamilton8.sh
source ../../../build/set_env.sh

<run_command>
