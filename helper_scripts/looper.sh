#!/bin/bash
#SBATCH --job-name=looper
#SBATCH --partition=colella
#SBATCH --time=1:00:00
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=output_looper.%j.txt
#SBATCH --reservation=b686w673_43

for replicate in 1 2 3 4 5 6 7 8 9 10
do
  for subsample in 1 2 3 4 5 6 7 8 9 10
  do
    vcf=`ls input/*_REPLICATE_"$replicate"_subsample_"$subsample".vcf`
    logfile=`ls input/*logfile_REPLICATE_"$replicate".txt`
    ./split_appended_vcf.sh $vcf $logfile
  done
done
