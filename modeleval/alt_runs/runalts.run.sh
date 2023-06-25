#!/bin/bash
#SBATCH -J runalts
#SBATCH -o runalts.out
#SBATCH -e runalts.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 0-06:00
#SBATCH -p huce_ice
#SBATCH --mem=105G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2021a-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "runalts;"
