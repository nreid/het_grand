#!/bin/bash

#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem=30000

source /etc/profile.d/modules.sh

module load samtools/1.3
module load xz

FB=~/bin/freebayes/bin/freebayes 
BED=~/bin/bedtools2/bin/bedtools
BAM=~/bin/bamtools/bin/bamtools
BGZ=~/bin/htslib/bgzip
REFGEN=~/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta

$FB -f $REFGEN -L <(cat /scratch/nmr15102/variants/meta/hetgrand.all.list | grep -v "KC-9") -t /scratch/nmr15102/variants/hetgrand_sex_cnv.txt -k  | \
$BGZ -c >/scratch/nmr15102/variants/hetgrand.cnv.vcf.gz

