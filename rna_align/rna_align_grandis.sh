#!/bin/bash

# align jane park's trimmed merged grandis rna-seq fastq files

fq1=$1


BWA="/home/nmr15102/bin/bwa/bwa mem"
genome=/home/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz


# define input
echo $fq1
fq2=$(echo $fq1 | sed 's/_R1_/_R2_/')
echo $fq2


# set bwa output variables
bamp=$(echo $fq1 | sed 's/.fq.gz/.bam/' | sed 's/_R[12]_merged//' | sed 's/fastqs/alignments/')

# read group info
run=merged
lib=merged
samp=$(echo $fq1 | grep -oP "(?<=/Sample_)[^/]+" | sed 's/AWJRDD00._//')
lane=merged
sub=merged
rg=$(echo \@RG\\tID:$samp\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$samp)


# run bwa
$BWA -R $rg $genome $fq1 $fq2 | \
samtools sort -O bam -T $bamp.temp >$bamp

