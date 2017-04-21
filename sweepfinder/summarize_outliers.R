library(magrittr)
library(dplyr)
library(stringr)
library(ape)
# a bioconductor library used for venn diagrams
library(limma)


# convert alpha to approximate selection coefficient
als <- function(alpha){

	10^-8 * log(100000) / alpha

}


# wherever your data are
load("~/Dropbox/Public/het_grand_outliers.RData")

# objects:
	# admixwin - outlier regions in the admixed heteroclitus population
	# northwin - outlier regions in the northern heteroclitus population
	# southwin - outlier regions in the southern heteroclitus population
	# grandwin - outlier regions in grandis

	# wi - merged outlier windows and corresponding max CLR and alpha statistics for each population
	# sweep - whole genome sweepfinder output (sweepa and sweepxy contain autosomes and xy)
	# als - a function to convert alpha values to approximate selection coefficients

	# fstnull - ten thousand random 10kb windows with Fst
	# fst3kb - all six pairwise Fst values for every window in 'wi' with greater than 15 variants
		# first column is the number of variants in the window, the rest are self-explanatory

fst3kb <- as.matrix(fst3kb)
fstnull <- as.matrix(fstnull)

fst3kb[which(fst3kb <= 0)] <- 0.0001
fstnull[which(fstnull <= 0)] <- 0.0001

# make a venn diagram if you can install limma 
(!is.na(wi[,seq(4,11,2)])) %>% vennCounts() %>% vennDiagram()
wi[rowSums(is.na(wi[,seq(4,11,2)]))==0,]

# which regions are outliers in only north and admixed populations?
noad <- !is.na(wi[,4]) & !is.na(wi[,6]) & is.na(wi[,8]) & is.na(wi[,10])

# what does Fst look like in those windows? 

fst3kb[noad,-1] %>% boxplot()

plot(tab[,2],smooth(tab[,3])/smooth(tab[,4]),pch=20,cex=.2)

