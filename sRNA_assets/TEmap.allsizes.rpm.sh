#!/usr/bin/env bash
#SBATCH --job-name=TEmap
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-31
#SBATCH --time=0:30:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4

cd /groups/nodine/lab/members/Ranjith/mCHH/all24ntsRNAs
ml bedtools/2.27.1-foss-2018b
forTE=$(ls *.bed |  head -n ${SLURM_ARRAY_TASK_ID}| tail -n 1)
bedtools map -a TEs.siRNA.clusterID.sorted -b $forTE -c 5 -o sum -g ara.genome -null 0 > $(basename -s .forTE $forTE).TEmapped
