library(tidyverse)

lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

val <- read.table("~/projects/het_grand_data/popsites.1kb_win.bed.gz",stringsAsFactors=FALSE)
val <- val[lord, ]

d <- duplicated(val[,1])
lgscaf <- cbind(scaf=val[!d,1],chr=lift[!d,1])
rownames(lgscaf) <- lgscaf[,1]

sx <- read.table("sex_chisq.txt.gz",stringsAsFactors=FALSE)
fai <- read.table("GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai",stringsAsFactors=FALSE)
rownames(fai) <- fai[,1]

scount <- group_by(sx,V1) %>% summarize(.,hetcount=sum(V5 < 0.00001,na.rm=TRUE),gcount=sum(V6 < 0.00001,na.rm=TRUE)) %>% data.frame()

scount <- cbind(scount,scaflen=fai[scount[,1],2]/1000000,hden=scount[,2]/fai[scount[,1],2]*1000000)

scount[order(scount[,2]/scount[,4]),]

lgscaf[scount[scount[,2] > 0,1],2] %>% table() %>% sort() %>% data.frame()

fai[scount[scount[,2] > 0,1],2] %>% sum()

hden <- scount[,2]/fai[scount[,1],2]
gden <- scount[,3]/fai[scount[,1],2]

hist(hden,breaks=seq(0,0.015,0.0002))
hist(gden,breaks=seq(0,0.015,0.0002))


scount[scount[,5] > 1000 & scount[,4] > 0.02,4] %>% sum()