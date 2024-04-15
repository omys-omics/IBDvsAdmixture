#!/bin/bash
#SBATCH --job-name=xVal
#SBATCH --partition=colella
#SBATCH --time=100:00:00
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=output_xVal.%j.txt
#SBATCH --reservation=b686w673_43

module load R

#file=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../01.input.files.txt)
file=$(head -1 ../01.input.files.txt)

prefix=`basename $file | cut -f1,2,3 -d"."`

R --vanilla -f 02.crossValidate.R --args $file $prefix