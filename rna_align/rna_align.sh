#!/bin/bash

fq1=$1


# set some variables
trimmo="java -jar /home/nmr15102/bin/trimmomatic-0.36/trimmomatic-0.36.jar"
BWA="/home/nmr15102/bin/bwa/bwa mem"
genome=/scratch/nmr15102/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz


# define input and output for trimmomatic
echo $fq1
fq2=$(echo $fq1 | sed 's/_R1_/_R2_/')
echo $fq2

ofq1=$(echo $fq1 | sed 's/fastq/pt\.fastq/')
ofq1u=$(echo $fq1 | sed 's/fastq/ut\.fastq/')
ofq2=$(echo $fq2 | sed 's/fastq/pt\.fastq/')
ofq2u=$(echo $fq2 | sed 's/fastq/ut\.fastq/')
echo $ofq1
echo $ofq1u
echo $ofq2
echo $ofq2u

# run trimmomatic
$trimmo PE -threads 4 -phred33 $fq1 $fq2 $ofq1 $ofq1u $ofq2 $ofq2u ILLUMINACLIP:~/bin/trimmomatic-036/NEBnextAdapt.fa:2:30:10 LEADING:5 TRAILING:5

# set bwa output variables
bamp=$(echo $ofq1 | sed 's/_R[12].*/.pt.bam/')
bamu2=$(echo $ofq1u | sed 's/_R[12].*/.ut1.bam/')
bamu1=$(echo $ofq1u | sed 's/_R[12].*/.ut2.bam/')

# read group info
run=$(echo $fq1 | grep -oP '(?<=tolerance_rnaseq/)[^/]+')
lib=$(echo $fq1 | grep -oP '(?<=Sample_)[^/]+')
samp=$(echo $fq1 | grep -oP "(?<=/Sample_)[^/]+")
lane=$(echo $fq1 | grep -oP "L00.")
sub=$(echo $fq1 | grep -oP "(?<=L00._R._)...")
rg=$(echo \@RG\\tID:$samp.$lib.$run.$lane.$sub\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$samp)



$BWA -R $rg $genome $ofq1 $ofq2 | \
samtools sort -O bam -T $bamp.temp >$bamp

$BWA -R $rg $genome $ofq1u | \
samtools sort -O bam -T $bamu1.temp >$bamu1

$BWA -R $rg $genome $ofq2u | \
samtools sort -O bam -T $bamu2.temp >$bamu2


rm $ofq1 $ofq2 $ofq1u $ofq2u