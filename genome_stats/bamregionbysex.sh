#!/bin/bash

# a script to output aggregated bam files for a given region
# split across males, females, etc

module load samtools/1.3
BAM=~/bin/bamtools/bin/bamtools

GF=/scratch/nmr15102/variants/meta/grand_females.list
GFS=/scratch/nmr15102/variants/meta/grand_singlecopyF.list
GM=/scratch/nmr15102/variants/meta/grand_males.list

HF=/scratch/nmr15102/variants/meta/het_females.list
HM=/scratch/nmr15102/variants/meta/het_males.list

HRF=/scratch/nmr15102/variants/meta/het_rna_females.list
HRM=/scratch/nmr15102/variants/meta/het_rna_males.list

GR=/scratch/nmr15102/variants/meta/grand_rnaseq.list

region=$1

$BAM merge -list $GF -region $region >${region}_gf.bam
$BAM merge -list $GFS -region $region >${region}_gfs.bam
$BAM merge -list $GM -region $region >${region}_gm.bam

$BAM merge -list $HF -region $region >${region}_hf.bam
$BAM merge -list $HM -region $region >${region}_hm.bam

$BAM merge -list $HRF -region $region >${region}_hrf.bam
$BAM merge -list $HRM -region $region >${region}_hrm.bam

$BAM merge -list $GR -region $region >${region}_gr.bam