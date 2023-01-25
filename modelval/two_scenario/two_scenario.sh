#!/bin/bash
#SBATCH -J two_scenario
#SBATCH -o two_scenario.out
#SBATCH -e two_scenario.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 0-01:30
#SBATCH -p huce_cascade
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2018b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "runtwoscenariotest_23_01_25.m;"
