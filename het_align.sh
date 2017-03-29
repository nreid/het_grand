#!/bin/bash

#1 fastq1
#2 fastq2
#3 ID table

# need to load parallel and samtools

bwagenind=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz

BWA=~/bin/bwa/bwa
SBL=~/bin/samblaster/samblaster
BED=~/bin/bedtools2/bin/bedtools

#TEM=/scratch/nmr15102/heteroclitus_bam/BP-29.all.bam

sam=$(echo $1 | sed 's/.all.*//' | sed 's/.*\///')

fq1=$sam.1.fastq
fq2=$sam.2.fastq

fqdir=/scratch/nmr15102/heteroclitus_fastq
outdir=/scratch/nmr15102/heteroclitus_align

cmdline1=samtools\ sort\ -T\ $fqdir/$sam.temp\ -n\ $1
cmdline2=$BED\ bamtofastq\ -i\ stdin\ -fq\ $fqdir/$fq1\ -fq2\ $fqdir/$fq2

rg=$(echo \@RG\\tID:$sam\\tPL:Illumina\\tPU:x\\tLB:combined\\tSM:$sam)

echo $sam $fq1 $fq2 $rg

echo $cmdline1
echo $cmdline2

$cmdline1 | $cmdline2

cmdline3=$BWA\ mem\ $bwagenind\ -t\ 2\ -R\ $rg\ $fqdir/$fq1\ $fqdir/$fq2
echo $cmdline3

# execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
$cmdline3 | $SBL -e -d $outdir/$sam.disc.sam -s $outdir/$sam.split.sam | samtools view -S -h -u - | samtools sort -T $outdir/$sam - >$outdir/$sam.bam

