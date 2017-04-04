library(magrittr)
library(dplyr)

fs <- list.files("~/projects/het_grand_data",pattern="\\.1kb.*cov",full.names=TRUE)[-7]

fai <- read.table("~/projects/het_grand_data/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai",stringsAsFactors=FALSE)

bc <- read.table("~/projects/het_grand_data/ACGT_1kb.bed",stringsAsFactors=FALSE)
colnames(bc) <- c(
	"scaf",
	"start",
	"end",
	"at",
	"gc",
	"a",
	"c",
	"g",
	"t",
	"n",
	"o",
	"len"
	)

ct <- list()
for(i in 1:length(fs)){

	ct[[i]] <- read.table(fs[i],stringsAsFactors=FALSE)

	}

ct2 <- do.call(cbind,ct)

ct3 <- ct2[,seq(from=6,to=13*6,by=6)]

cn <- gsub(".*/","",fs) %>% gsub("\\..*","",.)

colnames(ct3) <- cn

ct4 <- cbind(bc[,1:10],ct3)

genomestats <- cbind(bc[,c(1:3,6,7,8,9,10,4,5)],bases=rowSums(bc[,6:9]),ct3)

save(genomestats,fai,file="~/Dropbox/Public/coveragedataV2.RData")


