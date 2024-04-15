#!/bin/bash
#SBATCH --job-name=conStruct
#SBATCH --partition=colella
#SBATCH --time=600:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-15
#SBATCH --output=output_conStruct.%j.txt
#SBATCH --reservation=b686w673_43

module load R

file=$(sed -n "$SLURM_ARRAY_TASK_ID"p 01.input.files.txt)
#file=$(head -1 01.input.files.txt)

R --vanilla -f 01.conStruct_analyses.R --args $file