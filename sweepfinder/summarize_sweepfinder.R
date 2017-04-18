library(magrittr)
library(dplyr)
library(stringr)
# a bioconductor library used for venn diagrams
library(limma)

# a function to find selective sweep intervals and identify their peaks
createwindows <- function(sw,CLR_col,alpha_col,thresh,mergedist=2001){
	
	out <- c()
	sw <- sw[,c("Scaffold","pos","pos","pos",CLR_col,alpha_col)]
	colnames(sw)[2:4] <- c("start","end","peak_loc")

	keep <- which(sw[,CLR_col] > thresh)
	out <- sw[keep[1],]

	j <- 1
	minheight <- 0.5 * out[1,5]

	for(i in keep[2]:dim(sw)[1]){

		if(sw[i,CLR_col] > thresh){

				if( (sw[i,3]-out[j,3] < mergedist) & sw[i,1]==out[j,1] & (sw[i,5] > minheight) ){
		
					out[j,3] <- sw[i,3]
		
					if( out[j,5] < sw[i,5] ) {out[j,5] <- sw[i,5]; out[j,4] <- sw[i,4]}
		
					out[j,6] <- min(c(out[j,6],sw[i,6]))
					minheight <- out[j,5] * 0.5
					}
		
				
				if( sw[i,1]!=out[j,1] | ((sw[i,3]-out[j,3] > mergedist) & (sw[i,5] > minheight) ) ){
					out <- rbind(out,sw[i,])
					j <- j + 1 
					}
			}else{ minheight <- thresh }

		}
	colnames(out)[2:3] <- c("start","end")
	return(out)

}

# a function to identify co-localization of sweeps. 
mergesweeps <- function(win1,win2,win3,win4,buffer=20000){

	shared <- c()

	win <- rbind(win1[,c(1,4,4)],win2[,c(1,4,4)],win3[,c(1,4,4)],win4[,c(1,4,4)])
	win <- win[order(win[,1],win[,2]),]	

	out <- win[1,]
	j <- 1
	for(i in 2:dim(win)[1]){

		if( out[j,3] >= (win[i,2] - buffer) & out[j,1] == win[i,1] ) {
			out[j,3] <- win[i,3]
		}
		
		if( out[j,3] < (win[i,2] - buffer) | out[j,1] != win[i,1] ) {
			out <- rbind(out,win[i,])
			j <- j + 1
		}
	}

	cnames <- c(colnames(win1)[5:6],colnames(win2)[5:6],colnames(win3)[5:6],colnames(win4)[5:6])
	out <- cbind(out,NA,NA,NA,NA,NA,NA,NA,NA)
	colnames(out)[4:11] <- cnames
	colnames(out)[2:3] <- c("start","end")

	for(k in 1:dim(out)[1]){

		i1 <- win1[,1]==out[k,1] & win1[,4]>=out[k,2] & win1[,4]<=out[k,3,]
		i2 <- win2[,1]==out[k,1] & win2[,4]>=out[k,2] & win2[,4]<=out[k,3,]
		i3 <- win3[,1]==out[k,1] & win3[,4]>=out[k,2] & win3[,4]<=out[k,3,]
		i4 <- win4[,1]==out[k,1] & win4[,4]>=out[k,2] & win4[,4]<=out[k,3,]

		if( any(i1) ){ out[k,4:5]   <- c(max(win1[i1,5]),min(win1[i1,6])) }
		if( any(i2) ){ out[k,6:7]   <- c(max(win2[i2,5]),min(win2[i2,6])) }
		if( any(i3) ){ out[k,8:9]   <- c(max(win3[i3,5]),min(win3[i3,6])) }
		if( any(i4) ){ out[k,10:11] <- c(max(win4[i4,5]),min(win4[i4,6])) }

	}

	return(out)

}

toppeaks <- function(win,col,n,inv=FALSE){

	r <- rank(-win[,col])
	if(inv) {r <- rank(win[,col]) }

	win[r <= n,]

}

# convert alpha to approximate selection coefficient
als <- function(alpha){

	10^-8 * log(100000) / alpha

}

# create intervals for Fst calculation
# x[,1]=scaffold,x[,2]=start,x[,3]=end
intervals <- function(x, width=10000){

	out <- c()
	for(i in 1:dim(x)[1]){
		
		iwid <- x[i,3]-x[i,2]
		if(iwid < width){
			buff <- (width-iwid) / 2
			x[i,2] <- x[i,2] - buff
			if(x[i,2] < 1){ x[i,2] <- 1 }
			x[i,3] <- x[i,3] + buff
			out <- c(out,paste(x[i,1],":",x[i,2],"-",x[i,3],collapse="",sep=""))
		}else{ out <- c(out,paste(x[i,1],":",x[i,2],"-",x[i,3],collapse="",sep="")) }

	}

	return(out)
}

# wherever your data are
load("~/Dropbox/Public/sweepfinder.RData")

colnames(sweep)[2] <- "pos"

xyscafs <- maxy[maxy[,2]=="xy",1]
sweepxy <- sweep[sweep[,1] %in% xyscafs,]

sweepa <- sweep[!(sweep[,1] %in% xyscafs),]

northwin <- createwindows(sweepa,"North_CLR","North_alpha",200,20000)
admixwin <- createwindows(sweepa,"Admix_CLR","Admix_alpha",200,20000)
southwin <- createwindows(sweepa,"South_CLR","South_alpha",200,20000)
grandwin <- createwindows(sweepa,"Grand_CLR","Grand_alpha",200,20000)

wi <- mergesweeps(northwin,admixwin,southwin,grandwin)
(!is.na(wi[,seq(4,11,2)])) %>% vennCounts() %>% vennDiagram()
wi[rowSums(is.na(wi[,seq(4,11,2)]))==0,]

# write out intervals to be analyzed
# write.table(intervals(wi,width=3000),file="het_grand_data/outlier_intervals3kb.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

cnames <- c("count","north.south","north.admix","north.grand","south.admix","south.grand","admix.grand")
fstnull <- read.table("het_grand_data/fst.null.out",stringsAsFactors=FALSE)
colnames(fstnull) <- cnames
fst10kb <- read.table("het_grand_data/fst.outlier.10kb.out",stringsAsFactors=FALSE)
colnames(fst10kb) <- cnames
fst3kb <- read.table("het_grand_data/fst.outlier.3kb.out",stringsAsFactors=FALSE)
colnames(fst3kb) <- cnames

# wi <- mergesweeps(
# 	toppeaks(northwin,5,200),
# 	toppeaks(admixwin,5,200),
# 	toppeaks(southwin,5,200),
# 	toppeaks(grandwin,5,200))

# (!is.na(wi[,seq(4,11,2)])) %>% vennCounts() %>% vennDiagram()
# wi[rowSums(is.na(wi[,seq(4,11,2)]))==0,]

par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(3,3,1,1))
noad <- !is.na(wi[,4]) & !is.na(wi[,6])
boxplot(x=list(fst3kb[noad,"north.admix"],fstnull[,"north.admix"],fst3kb[noad,"north.south"],fstnull[,"north.south"],fst3kb[noad,"south.admix"],fstnull[,"south.admix"]),col=c("red","blue"))

noad <- is.na(wi[,4]) & !is.na(wi[,6])
boxplot(x=list(fst3kb[noad,"north.admix"],fstnull[,"north.admix"],fst3kb[noad,"north.south"],fstnull[,"north.south"],fst3kb[noad,"south.admix"],fstnull[,"south.admix"]),col=c("red","blue"))

so <- !is.na(wi[,8]) & is.na(wi[,4]) & is.na(wi[,6])
boxplot(x=list(fst3kb[so,"north.admix"],fstnull[,"north.admix"],fst3kb[so,"north.south"],fstnull[,"north.south"],fst3kb[so,"south.admix"],fstnull[,"south.admix"]),col=c("red","blue"))


noad <- !is.na(wi[,4]) & !is.na(wi[,6])
hist(fst3kb[noad,3],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank",xlim=c(0,0.4))
hist(fstnull[,3],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

hist(fst3kb[noad,2],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank")
hist(fstnull[,2],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(3,3,1,1))
noad <- !is.na(wi[,4]) & is.na(wi[,6])
hist(fst3kb[noad,3],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank",xlim=c(0,0.4))
hist(fstnull[,3],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

hist(fst3kb[noad,2],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank")
hist(fstnull[,2],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(3,3,1,1))
noad <- is.na(wi[,4]) & !is.na(wi[,6])
hist(fst3kb[noad,3],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank",xlim=c(0,0.4))
hist(fstnull[,3],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

hist(fst3kb[noad,2],freq=FALSE,breaks=br,col=rgb(0,0,1,.2),lty="blank")
hist(fstnull[,2],freq=FALSE,breaks=br,col=rgb(1,0,0,.2),add=TRUE,lty="blank")

noad <- is.na(wi[,4]) & !is.na(wi[,6])



vc <- cbind(
	North=as.numeric(sweepa[,"North_CLR"]>200),
	Admix=as.numeric(sweepa[,"Admix_CLR"]>200),
	South=as.numeric(sweepa[,"South_CLR"]>100),
	Grand=as.numeric(sweepa[,"Grand_CLR"]>100)
	)

par(mfrow=c(1,3))
vc1 <- vennCounts(vc)
vennDiagram(vc1)
vc1 <- vennCounts(vc[,-4])
vennDiagram(vc1)
vc1 <- vennCounts(vc[,-1])
vennDiagram(vc1)



par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[,"North_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
abline(h=200)
plot(sweepa[,"Admix_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
abline(h=200)
plot(sweepa[,"South_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
abline(h=200)
plot(sweepa[,"Grand_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
abline(h=200)

dev.new()
yl=c(0,0.3)
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(als(sweepa[,"F_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]),ylim=yl)
#abline(h=500)
plot(als(sweepa[,"SH_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]),ylim=yl)
#abline(h=500)
plot(als(sweepa[,"KC_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]),ylim=yl)
#abline(h=500)
plot(als(sweepa[,"BB_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]),ylim=yl)
#abline(h=500)

dev.new()
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(als(sweepa[,"BB_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(als(sweepa[,"VB_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(als(sweepa[,"GB_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(als(sweepa[,"BNP_alpha"]),pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)

par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[,"BB_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(sweepa[,"VB_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(sweepa[,"GB_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)
plot(sweepa[,"BNP_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
#abline(h=500)


plot(sweepa[,"North_CLR"],1/sweepa[,"North_alpha"],pch=20,cex=.2,col=factor(sweepa[,1]))


par(mfrow=c(1,1))
plot(sweepa[,"BP_CLR"],sweepa[,"F_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
plot(sweepa[,"BP_CLR"]-sweepa[,"F_CLR"],pch=20,cex=.2,col=factor(sweepa[,1]))
sweepa[(sweepa[,"BP_CLR"]-sweepa[,"F_CLR"])>800,1] %>% table() %>% sort()


par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[,"GB_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=100)
plot(sweepa[,"BNP_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"BB_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"VB_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)

om[sweepa[sweepa[,"North_CLR"]>500,1] %>% table() %>% sort() %>% names(),]

scaf <- "NW_012234474.1"
subc <- sweepa[,1]==scaf
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[subc,2],sweepa[subc,"F_CLR"],pch=20,cex=.2,type="l")
points(northwin[northwin[,1]==scaf,4],northwin[northwin[,1]==scaf,5],col="red")
abline(h=200)
plot(sweepa[subc,2],sweepa[subc,"BP_CLR"],pch=20,cex=.2,type="l")
points(admixwin[admixwin[,1]==scaf,4],admixwin[admixwin[,1]==scaf,5],col="red")
abline(h=200)
plot(sweepa[subc,2],sweepa[subc,"GB_CLR"],pch=20,cex=.2,type="l")
points(southwin[southwin[,1]==scaf,4],southwin[southwin[,1]==scaf,5],col="red")
abline(h=100)
plot(sweepa[subc,2],sweepa[subc,"BB_CLR"],pch=20,cex=.2,type="l")
points(grandwin[grandwin[,1]==scaf,4],grandwin[grandwin[,1]==scaf,5],col="red")
abline(h=100)

par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[subc,2],als((sweepa[subc,"SH_alpha"])),pch=20,cex=.2,type="l")
abline(h=200)
plot(sweepa[subc,2],als((sweepa[subc,"NYC_alpha"])),pch=20,cex=.2,type="l")
abline(h=200)
plot(sweepa[subc,2],als((sweepa[subc,"F_alpha"])),pch=20,cex=.2,type="l")
abline(h=100)
plot(sweepa[subc,2],als((sweepa[subc,"BP_alpha"])),pch=20,cex=.2,type="l")
abline(h=100)


par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(log(sweepa[,"North_alpha"]),pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(log(sweepa[,"Admix_alpha"]),pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(log(sweepa[,"South_alpha"]),pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(log(sweepa[,"Grand_alpha"]),pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)