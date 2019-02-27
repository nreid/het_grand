#!/bin/bash

#number of cores to use 
CORES=30		
KMERSIZE=31

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

	# run jellyfish
	${jellyfishDir}/jellyfish count -C \
	-o $outdir/${OUTPREFIX}_kmers/tmp \
	-m ${KMERSIZE} \
	-t ${CORES} \
	-s 20G \
	<(samtools view -bhF 256 $file | $BED bamtofastq -i /dev/stdin -fq /dev/stdout )

	COUNT=$(ls $outdir/${OUTPREFIX}_kmers/tmp* | wc -l)

	if [ $COUNT -eq 1 ]
	then
 		mv $outdir/${OUTPREFIX}_kmers/tmp_0 $outdir/${OUTPREFIX}_kmers_jellyfish
	else
		${jellyfishDir}/jellyfish merge -o $outdir/${OUTPREFIX}_kmers_jellyfish $outdir/${OUTPREFIX}_kmers/tmp*
	fi
	rm -rf $outdir/${OUTPREFIX}_kmers
	
	COUNT=$(ls $outdir/${OUTPREFIX}_kmers_jellyfish | wc -l)

	if [ $COUNT -eq 1 ]
	then

		${jellyfishDir}/jellyfish histo -f -o $outdir/${OUTPREFIX}.kmers.hist.csv -t ${CORES} $outdir/${OUTPREFIX}_kmers_jellyfish
		awk '{print $2"\t"$1}' $outdir/${OUTPREFIX}.kmers.hist.csv > $outdir/${OUTPREFIX}_tmp
		mv $outdir/${OUTPREFIX}_tmp $outdir/${OUTPREFIX}.kmers.hist.csv

		awk -f ${hawkDir}/countTotalKmer.awk ${OUTPREFIX}.kmers.hist.csv >> ${outdir}/total_kmer_counts.txt

		CUTOFF=1 
		echo $CUTOFF > $outdir/${OUTPREFIX}_cutoff.csv


		${jellyfishDir}/jellyfish dump -c -L `expr $CUTOFF + 1` $outdir/${OUTPREFIX}_kmers_jellyfish > $outdir/${OUTPREFIX}_kmers.txt 
		${sortDir}/sort --parallel=${CORES} -n -k 1 $/outdir/${OUTPREFIX}_kmers.txt > $outdir/${OUTPREFIX}_kmers_sorted.txt
	
		rm $outdir/${OUTPREFIX}_kmers_jellyfish	
		rm $outdir/${OUTPREFIX}_kmers.txt		
			
		echo "$outdir/${OUTPREFIX}_kmers_sorted.txt" >> ${dir}/sorted_files.txt
		
	fi

done
