#!/bin/bash

# this script edits every input file 
	# so the first variant has a sample size equal to 
	# the length of the estimated allele frequency spectrum - 1
	# it also reduces the sample size of anything in excess of the spectrum
	# this is hacky, but it only affects a few sites and is easier than re-writing all the output. 

dim=$(wc -l globalspec/$1\_*spect | cut -f 1 -d ' ')
# echo $dim
dim=$(expr $dim - 1)

for file in $(ls input/*.1_$1_*)
do
echo $1 $file
awk -v dim="$dim" '{OFS="\t"}{if (NR==2) $3=dim}{if ($3>dim) $3=dim}{print}' $file >tmp.$pop && mv tmp.$pop $file
done