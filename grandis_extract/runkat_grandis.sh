#!/bin/bash



module load samtools/1.3

BED=~/bin/bedtools2/bin/bedtools

# for whatever reason, kat won't run unless I change these paths
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/nmr15102/bin/KAT/deps/boost/build/lib/
LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH

KAT=~/bin/KAT/src/kat

# jellyfish hash used to extract reads
# use hawk kmers
	# $JF count -m 31 -s 150000 -o case_test.jf case_test.fasta
# use all kmers from grandis genomic cnv regions
	# these bam files are in grandis_cnv, and contain all reads mapping to the 9 largest grandis CNV intervals
	# $JF count -m 31 -s 150000 -o grandis_cnv.jf <(samtools merge -fu -b <(ls *_g*bam | cat ) /dev/stdout | $BED bamtofastq -i /dev/stdin -fq /dev/stdout)

# jh=/scratch/nmr15102/hawktest/case_test.jf
jh=/scratch/nmr15102/grandis_cnv/grandis_cnv.jf


# sort bam by query name
# turn it into fastq
# extract reads by k-mer
# test file /scratch/nmr15102/grandis_align/BU000004_VB_B.bam

BAM=$1
BS=$(echo $BAM | sed 's/.*\///' | sed 's/.bam//')

samtools sort -n -O bam $BAM >${BS}_sorted.bam

samtools view -bh -F 2304 ${BS}_sorted.bam | \
$BED bamtofastq -i /dev/stdin -fq ${BS}_1.fastq -fq2 ${BS}_2.fastq

$KAT filter seq \
-o ${BS}_out \
--seq ${BS}_1.fastq \
--seq2 ${BS}_2.fastq \
$jh

rm ${BS}_sorted.bam ${BS}_1.fastq ${BS}_2.fastq

