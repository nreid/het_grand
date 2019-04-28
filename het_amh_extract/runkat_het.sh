#!/bin/bash



module load samtools/1.3

BED=~/bin/bedtools2/bin/bedtools

# for whatever reason, kat won't run unless I change these paths
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/nmr15102/bin/KAT/deps/boost/build/lib/
LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH

KAT=~/bin/KAT/src/kat

# jellyfish hash used to extract reads
# use all kmers from amh region
	# these bam files are in assemble_amh/bam/, and contain all reads mapping to amh for RNA-seq 
	# $JF count -m 31 -s 150000 -o het_amh.jf <(samtools merge -fu -b <(ls amh_het*bam | cat ) /dev/stdout | $BED bamtofastq -i /dev/stdin -fq /dev/stdout)

# jh=/scratch/nmr15102/hawktest/case_test.jf
jh=/scratch/nmr15102/het_amh_extracts/het_amh.jf


# sort bam by query name
# turn it into fastq
# extract reads by k-mer

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

