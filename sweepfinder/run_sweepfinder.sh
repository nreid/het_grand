#!/bin/bash

SF2=~/bin/SF2/SweepFinder2

indir=/scratch/nmr15102/sweepfinder/input
griddir=/scratch/nmr15102/sweepfinder/gridfiles
outdir=/scratch/nmr15102/sweepfinder/output

gridfile=$(echo $1 | grep -oP ".*.1").grid
outfile=$(echo $1 | sed 's/geno/out/')

$SF2 -su $griddir/$gridfile $indir/$1 $outdir/$outfile