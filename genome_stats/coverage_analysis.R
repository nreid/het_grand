library(dplyr)
library(magrittr)
library(ggplot2)

load("~/Dropbox/Public/coveragedataV2.RData")

# fai is a table 
	# column 1: scaffold name
	# column 2: scaffold length

# genomestats is a table of 5kb windows and some statistics:
	# columns 1-3: window location
	# columns 4-11: base composition of reference genome
		# "n" is number of N bases
		# "bases" is number of ACGT bases. silly, I know. 
	# columns 12-24: number of 100bp reads mapped with Q>=30
		# captials are grandis, lower case are heteroclitus
		# er and kc are "southern", the rest are northern
		# sh nyc are admixed heteroclitus
		# bp and fp are northern heteroclitus

# the data are in a ridiculous order. first we'll reorder the table
	# by scaffold size


rownames(fai) <- fai[,1]
ord <- order(-fai[genomestats[,1],2],genomestats[,1],genomestats[,2])

genomestats <- genomestats[ord,]
fai <- fai[order(-fai[,2]),]


# what is the median per-base coverage for each population?
subw <- genomestats[,"bases"] > 50

g2 <- genomestats[subw,]

# calculate quantiles for each population
q50 <- c()
for(i in 12:length(g2[1,])){

	q <- g2[,i]/g2[,"bases"]*100
	q <- quantile(q, probs=seq(0,1,0.005))
	q50 <- cbind(q50, q)

}
colnames(q50) <- colnames(g2)[12:length(g2[1,])]


# create index vector to exclude windows in 2% tails of coverage distribution
subwih <- c()
subwil <- c()
for(i in 12:length(g2[1,])){

	pbcov <- g2[,i]/g2[,"bases"]*100
	subwih <- cbind(subwih,(pbcov >= q50["98.5%",i-11]))
	subwil <- cbind(subwil,(pbcov <= q50["1.5%",i-11]))

}

# basic stats on included/excluded windows
subwi <- (!subwih & !subwil)
rowSums(subwi) %>% table()
sum(rowSums(subwi)>=13)/dim(subwi)[1]
rowSums(!subwih) %>% table()
rowSums(!subwil) %>% table()

subwfh <- rowSums(!subwih)>=13
subwfl <- rowSums(!subwil)>=13
subwf <- rowSums(subwi)>=13

apply(g2[subwf,12:24],MAR=2,FUN=function(x){max(x)/median(x)})
apply(g2[subwf,12:24],MAR=2,FUN=function(x){min(x)/median(x)})

hist(log(g2[subwf,"BNP"]/g2[subwf,"bases"]*100,base=10),breaks=100,freq=FALSE,xlim=c(-2,6))
hist(log(g2[!subwfh,"BNP"]/g2[!subwfh,"bases"]*100,base=10),breaks=250,freq=FALSE,add=TRUE,col=rgb(1,0,0,.5))
hist(log(g2[!subwfl,"BNP"]/g2[!subwfl,"bases"]*100,base=10),breaks=250,freq=FALSE,add=TRUE,col=rgb(0,0,1,.5))

write.table(g2[!subwfh,1:3],row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",file="hicov.1kb.bed")
# then merge overlapping/bookended intervals in bedtools:
	# ~/bin/bedtools/bin/bedtools merge -i hicov.1kb.bed -d 0 | sort -k1,1 -k2,2n >hicov.merge.bed
	# then move the NCXXXX line to the end. 

# We evaluated coverage in 1kb windows for each of thirteen populations. 
# We excluded windows that fell in the 1.5% tails of any of the thirteen populations. 
# This resulted in a total of 27581 (2.9%) 1kb windows excluded for high coverage and 32427 (or 3.4% excluded for low coverage). 
# We retained no windows with greater than 1.98x or less than 0.4% the median coverage. 






st <- subwf & g2[,1] == "Scaffold0"

plot(g2[st,2],g2[st,"BP"]/g2[st,"n_bases"],pch=20,cex=.5)
points(g2[st,2],g2[st,"FP"]/g2[st,"n_bases"],pch=20,cex=.5,col="red")

plot(g2[st,"BP"]/g2[st,"n_bases"],g2[st,"FP"]/g2[st,"n_bases"],pch=20,cex=.5)

par(mfrow=c(2,3))
for(i in 12:length(g2[1,])){

	hist(g2[subwf,i]/g2[subwf,"n_bases"],breaks=100)

}




popn <- "ER"
popd <- "KC"
covrat <- g2[,popn]/g2[,popd]
covrat <- log(covrat,base=2) - log(quantile(covrat,prob=0.5),base=2)

qqn <- quantile(g2[,popn],prob=c(0.1,0.95))
qqd <- quantile(g2[,popd],prob=c(0.1,0.95))
subw4 <-  qqd %>% (function(x){g2[,popd] > x[1] & g2[,popd] < x[2]})

par(mfrow=c(1,1))
plot(covrat[subw4],pch=20,cex=.2,col=factor(g2[subw4,1]))
abline(h=c(-1,1))
qq <- quantile(covrat[subw4],prob=c(0.001,0.999))
qq %T>% abline(h=.,lwd=2,lty=2)

g2[subw4,][(covrat[subw4] < qq[1]),1] %>% table() %>% sort()

dev.new()

scaf <- "Scaffold10074"
sss <- g2[,1]==scaf
ssc <- subw4[sss]
ssr <- (covrat[sss] > qq[2] | covrat[sss] < qq[1]) & ssc

par(mfrow=c(3,1))
plot(g2[sss,2],covrat[sss],pch=20,col=ssr+1,cex=0.5)
abline(h=0);abline(h=qq,lty=2)
plot(g2[sss,2],g2[sss,popn],pch=20,cex=0.5,col=(2*ssc)+1)
abline(h=qqn,lty=2)
plot(g2[sss,2],g2[sss,popd],pch=20,cex=0.5,col=(2*ssc)+1)
abline(h=qqd,lty=2)

g2[sss,][ssr,c("scaf","start","end",popn,popd)] %>% cbind(.,covrat=covrat[sss][ssr])


# how much variation is there around that number?
# how well correlated is windowed coverage among populations?


# in review, 

# which windows are likely to contain repetitive elements?

# do any regions show evidence of "interesting" copy number variation?
	# with respect to the reference genome?
	# among populations? 
		# to test your thinking on this, 
		# there is certainly one large region in ER and 
		# one large region in BP that have important structural variants

# which windows should be excluded from analysis?
	# consider base content
	# number of non-N bases
	# relative levels of coverage within and among populations



