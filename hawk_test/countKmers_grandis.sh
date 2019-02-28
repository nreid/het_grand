#!/bin/bash

#number of cores to use 
CORES=24		
KMERSIZE=31

# INSTEAD OF A LOOP OVER ALL FILES
	# it would be easy to parallelize this
	# push the lines that append results from all files
	# to another script and run that when the k-mer
	# counting is finished.

# load samtools
module load samtools/1.3

#directory for read files 
dir=/scratch/nmr15102/grandis_align
# directory for output
outdir=/scratch/nmr15102/grandis_hawk
#directory where hawk is installed
hawkDir=~/bin/HAWK
#directory where jellyfish is installed
jellyfishDir=~/bin/HAWK/supplements/jellyfish-Hawk/bin/
#directory where parallel sort is installed
sortDir=~/bin/coreutils-8.30/src/
# bedtools
BED=~/bin/bedtools2/bin/bedtools


cd ${dir}

# loop through bam files, convert them to fastq, feed to jellyfish
for file in $(find $dir -name "*bam" | sort)
do
	OUTPREFIX=$(echo $file | sed 's/.bam//' | sed 's/.*\///')
	
	mkdir $outdir/${OUTPREFIX}_kmers

	# run jellyfish, converting bam to fastq
	${jellyfishDir}/jellyfish count -C \
	-o $outdir/${OUTPREFIX}_kmers/tmp \
	-m ${KMERSIZE} \
	-t ${CORES} \
	-s 2G \
	<(samtools view -bhF 256 $file | $BED bamtofastq -i /dev/stdin -fq /dev/stdout )

	# merge output files if there are > 1, then move file to another directory and rename
	COUNT=$(ls $outdir/${OUTPREFIX}_kmers/tmp* | wc -l)

	if [ $COUNT -eq 1 ]
	then
 		mv $outdir/${OUTPREFIX}_kmers/tmp_0 $outdir/${OUTPREFIX}_kmers_jellyfish
	else
		${jellyfishDir}/jellyfish merge -o $outdir/${OUTPREFIX}_kmers_jellyfish $outdir/${OUTPREFIX}_kmers/tmp*
	fi

	rm -rf $outdir/${OUTPREFIX}_kmers
	
	COUNT=$(ls $outdir/${OUTPREFIX}_kmers_jellyfish | wc -l)

	# assuming there is only one output file per input file...
	if [ $COUNT -eq 1 ]
	then

		# run jellyfish histo
			# max high count is 10k at default
			# this means # kmers in histogram is not 
			# the same as the # you get by summing 
			# the counts in the sorted.txt files
		${jellyfishDir}/jellyfish histo \
		-f \
		-o $outdir/${OUTPREFIX}.kmers.hist.csv \
		-t ${CORES} $outdir/${OUTPREFIX}_kmers_jellyfish

		# reformat the hist file, only print cols 1 & 2 (dunno if there are more)
		awk '{print $2"\t"$1}' $outdir/${OUTPREFIX}.kmers.hist.csv > $outdir/${OUTPREFIX}_tmp

		# rename
		mv $outdir/${OUTPREFIX}_tmp $outdir/${OUTPREFIX}.kmers.hist.csv

		# count total number of kmers, write to a 
			# the script used here must be modified, it excludes kmers that only appear once
			# THIS FILE CONTAINS INFO FOR ALL SAMPLES AND CAN'T BE PARALLELIZED
		awk -f ${hawkDir}/countTotalKmer.awk $outdir/${OUTPREFIX}.kmers.hist.csv >> ${outdir}/total_kmer_counts.txt

		# don't use this kmer frequency cutoff. edited below. 
		CUTOFF=1 
		echo $CUTOFF > $outdir/${OUTPREFIX}_cutoff.csv

		# run jellyfish dump, which outputs each kmer and its count, I think
		${jellyfishDir}/jellyfish dump \
		-c -L 1 \
		$outdir/${OUTPREFIX}_kmers_jellyfish > $outdir/${OUTPREFIX}_kmers.txt 

		# sort the kmer count file
		${sortDir}/sort --parallel=${CORES} -n -k 1 $outdir/${OUTPREFIX}_kmers.txt > $outdir/${OUTPREFIX}_kmers_sorted.txt
		
		# remove intermediate files	
		rm $outdir/${OUTPREFIX}_kmers_jellyfish	
		rm $outdir/${OUTPREFIX}_kmers.txt		
		
		# write sorted filename to a file
			# THIS FILE CONTAINS INFO FOR ALL SAMPLES AND CAN'T BE PARALLELIZED
		echo "$outdir/${OUTPREFIX}_kmers_sorted.txt" >> ${outdir}/sorted_files.txt
		
	fi

done
