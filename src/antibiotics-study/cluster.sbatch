#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o /home/sankarak/logs/slurm-%j.out  # Write the log in $SCRATCH

cd /home/sankarak/microbiome_plvm/src/antibiotics-study/
module load singularity/3.4;
singularity exec --bind /scratch/sankarak/data/,/scratch/sankarak/checkpoints/ --nv /scratch/sankarak/images/lda.sif Rscript lda.R --K 16