#!/bin/bash



module load samtools/1.3

BED=~/bin/bedtools2/bin/bedtools

# for whatever reason, kat won't run unless I change these paths
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/nmr15102/bin/KAT/deps/boost/build/lib/
LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH

KAT=~/bin/KAT/src/kat

# jellyfish hash used to extract reads
# use all kmers from amh region
	# $JF count -m 31 -s 150000 -o osmotic_amh.jf <($BAM merge -region NW_012234285.1:160000..164000 -list <(cat /scratch/nmr15102/variants/meta/rnabams.list  | grep osmotic) | $BED bamtofastq -i /dev/stdin -fq /dev/stdout)

# jh=/scratch/nmr15102/hawktest/case_test.jf
jh=/scratch/nmr15102/osmotic_amh_extracts/osmotic_amh.jf


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

