#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=SGPU,GGPU
#SBATCH --gres=gpu:1 --ntasks=1
#SBATCH --time=24:00:00

source activate omm8plumed
module load gcc/9.3.0 openmpi/3.1.6gcc9

export OPENMM_CPU_THREADS=1
export OMP_NUM_THREADS=1

python -u eqm.py > eqm.log
bash migraterst.sh

