#!/bin/bash
# Job name:
#SBATCH --job-name=001
#
# Partition:
#SBATCH --partition=savio2
#
#SBATCH --account=fc_bivarsurv
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=yi_li@berkeley.edu

module load r

R CMD BATCH --no-save sim.R sim2_1.Rout
