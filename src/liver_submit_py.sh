#!/bin/bash
#SBATCH --job-name tune_coverage_param 
#SBATCH --mail-user yh185@duke.edu
#SBATCH --mail-type ALL
#SBATCH -c 8
#SBATCH -t "3-0"
#SBATCH --output "tune_coverage_param".out
#SBATCH --error "tune_coverage_param".err 
#SBATCH --mem=128G–

file=$1
python $file
