#!/bin/bash
#SBATCH -J runalts
#SBATCH -o runalts.out
#SBATCH -e runalts.err
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -t 0-02:30
#SBATCH -p huce_intel
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --mail-user=amdur@g.harvard.edu

module load matlab/R2018b-fasrc01
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -r "runalts;"
