#!/bin/bash

# no idea which step conducts the bonferroni correction. 
# it must be during the initial hawk run which also outputs uncorrected p-values. 
# the R scripts appear to do the correction for population structure. 

# needed to recompile eigenstrat, 
# and recompiled version won't run without these modules loaded

module load blas/openblas-0.2.18
module load gsl/2.4
module load lapack/3.7.1

# script also calls a batch R script
module load r/3.3.3

# set directories
hawkDir=~/bin/HAWK
eigenstratDir=~/bin/HAWK/supplements/EIG6.0.1-Hawk

# set some metadata
isDiploid=1
noInd=$(cat sorted_files.txt | wc -l);

# prepare some input
$hawkDir/preProcess

cat case_total_kmers.txt control_total_kmers.txt > gwas_eigenstratX.total
cat case.ind control.ind > gwas_eigenstratX.ind

caseCount=$(cat case_sorted_files.txt | wc -l);
controlCount=$(cat control_sorted_files.txt | wc -l);

# run hawk
$hawkDir/hawk $caseCount $controlCount

# run eigenstrat smartpca
if [ "$isDiploid" == "0" ]; then
	$eigenstratDir/bin/smartpca -V -p $hawkDir/parfile.txt > log_eigen.txt
else
	$eigenstratDir/bin/smartpca -p $hawkDir/parfile.txt > log_eigen.txt
fi

# something else in eigenstrat
$eigenstratDir/bin/evec2pca.perl 10 gwas_eigenstrat.evec gwas_eigenstratX.ind gwas_eigenstrat.pca

tail -${noInd} gwas_eigenstrat.pca > pcs.evec

# sort by corrected p-values, make a file with top 200k
sort -g  -k 4 -t $'\t' case_out_w_bonf.kmerDiff > case_out_w_bonf_sorted.kmerDiff
head -200000 case_out_w_bonf_sorted.kmerDiff > case_out_w_bonf_top.kmerDiff

sort -g  -k 4 -t $'\t' control_out_w_bonf.kmerDiff > control_out_w_bonf_sorted.kmerDiff
head -200000 control_out_w_bonf_sorted.kmerDiff > control_out_w_bonf_top.kmerDiff

# these R scripts run an anova accounting for PC1 and PC2
# then output a new table e.g. "pvals_case_top.txt"
Rscript $hawkDir/log_reg_case.R
Rscript $hawkDir/log_reg_control.R

# attach new structure-corrected p-values to bonf corrected p-value files
paste pvals_case_top.txt case_out_w_bonf_top.kmerDiff  > pvals_case_top_merged.txt
sort -g -k 1 -t $'\t' pvals_case_top_merged.txt > pvals_case_top_merged_sorted.txt 

paste pvals_control_top.txt control_out_w_bonf_top.kmerDiff  > pvals_control_top_merged.txt
sort -g -k 1 -t $'\t' pvals_control_top_merged.txt > pvals_control_top_merged_sorted.txt

# presumably turns tables into fasta files. not many k-mers output though...
$hawkDir/convertToFasta    

# remove intermediary files. 
# not sure the bonferroni correction is great here. 
# can't do better (e.g. q-values) without full p-value distribution though. 
rm case_out_w_bonf.kmerDiff
rm case_out_wo_bonf.kmerDiff
rm control_out_w_bonf.kmerDiff
rm control_out_wo_bonf.kmerDiff 
