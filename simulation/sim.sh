#!/bin/bash
# Job name:
#SBATCH --job-name=804
#
# Partition:
#SBATCH --partition=savio3
#
#SBATCH --qos=biostat_savio3_normal
#SBATCH --account=co_biostat
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=72:00:00
#
# Number of nodes for use case:
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=yi_li@berkeley.edu

module load r

R CMD BATCH --no-save sim.R sim804.Rout
