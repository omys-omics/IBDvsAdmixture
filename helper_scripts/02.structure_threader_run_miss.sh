#!/bin/bash
#SBATCH --job-name=structure
#SBATCH --partition=colella
#SBATCH --time=200:00:00
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --output=output_structure.%j.txt
#SBATCH --reservation=b686w673_43

str_files=`ls str.miss.files/*.str`

for str in $str_files
do
  # format the structure file correctly
  awk '{ for(i=2; i<NF; i++) { printf "%s ", $i } print $NF }' $str > temp.str
  paste inds.txt pops.txt temp.str > $str
  rm temp.str
  # get the generation number from the file name
  gen=`basename $str | cut -f3 -d"."`
  # get the amount of missing data from the file name
  miss=`basename $str | cut -f6 -d"."`
  # count how many loci there are
  num_loci=`head -1 $str | grep -o "_" | wc -l`
  # edit mainparams file for correct number of loci
  sed -i "s/^#define NUMLOCI[[:space:]]\+[0-9]\+/#define NUMLOCI    ${num_loci}/" mainparams
  # run structure
  mkdir str.miss.output/01.output.IBD.gen.$gen.miss.0.$miss.200K.K2
  dos2unix $str
  structure_threader run -t 5 -i $str -o str.miss.output/01.output.IBD.gen.$gen.miss.0.$miss.200K.K2 --params mainparams -Klist 2 -R 5 -st ~/.local/bin/structure --log TRUE
done

# -i infile
# -o output directory
# --params File with run parameters
# -t threads
# -Klist '2 4 6' ['2 4 6' ...] List of Ks to calculate
# -R replicates
# -st path to structure
# --log enable logging