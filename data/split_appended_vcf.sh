#!/bin/bash

# argument 1: appended vcf file
# argument 2: logfile

# split the appended vcf file to separate files
csplit -f vcf $1 /"fileformat"/ '{*}'
rm vcf00 #get rid of empty vcf00

# get generation since contact from log file
newcount=0
gen=`grep generation -v $2 | cut -f2 -d","`
for item in $gen
do
  newcount=$(( $newcount+1 ))
  echo $item > temp.$newcount.txt
done

count=0
ls vcf* | while read vcf
do
  count=$(( $count+1 ))
  name=`cat temp.$count.txt`
  mv $vcf gen.$name.vcf
  rm temp.$count.txt
done

