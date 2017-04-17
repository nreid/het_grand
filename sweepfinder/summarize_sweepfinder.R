library(magrittr)
library(dplyr)
library(stringr)
# a bioconductor library used for venn diagrams
library(limma)

# a function to find selective sweep intervals and identify their peaks
createintervals <- function(sw,CLR_col,alpha_col,thresh,mergedist=2001){
	
	out <- c()
	sw <- sw[sw[,CLR_col] > thresh,c("Scaffold","pos","pos","pos",CLR_col,alpha_col)]
	colnames(sw)[2:4] <- c("start","end","peak_loc")
	out <- sw[1,]

	j <- 1
	for(i in 2:dim(sw)[1]){

		if( (sw[i,3]-out[j,3] < mergedist) & sw[i,1]==out[j,1] ){

			out[j,3] <- sw[i,3]

			if( out[j,5] < sw[i,5] ) {out[j,5] <- sw[i,5];out[j,4] <- sw[i,4]}

			out[j,6] <- min(c(out[j,6],sw[i,6]))
			}

		
		if( (sw[i,3]-out[j,3] > mergedist) | sw[i,1]!=out[j,1] ){
			out <- rbind(out,sw[i,])
			j <- j + 1 
			}

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

		i1 <- win1[,1]==out[k,1] & win1[,3]>=out[k,2] & win1[,2]<out[k,3,]
		i2 <- win2[,1]==out[k,1] & win2[,3]>=out[k,2] & win2[,2]<out[k,3,]
		i3 <- win3[,1]==out[k,1] & win3[,3]>=out[k,2] & win3[,2]<out[k,3,]
		i4 <- win4[,1]==out[k,1] & win4[,3]>=out[k,2] & win4[,2]<out[k,3,]

		if( any(i1) ){ out[k,4:5]   <- c(max(win1[i1,5]),min(win1[i1,6])) }
		if( any(i2) ){ out[k,6:7]   <- c(max(win2[i2,5]),min(win2[i2,6])) }
		if( any(i3) ){ out[k,8:9]   <- c(max(win3[i3,5]),min(win3[i3,6])) }
		if( any(i4) ){ out[k,10:11] <- c(max(win4[i4,5]),min(win4[i4,6])) }

	}

	return(out)

}


# wherever your data are
load("~/Dropbox/Public/sweepfinder.RData")

fai2 <- fai[]

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
plot(sweepa[,"North_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"Admix_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"South_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"Grand_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)

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

subc <- sweepa[,1]=="NW_012224652.1"
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[subc,2],sweepa[subc,"North_CLR"],pch=20,cex=.2)
abline(h=200)
plot(sweepa[subc,2],sweepa[subc,"Admix_CLR"],pch=20,cex=.2)
abline(h=200)
plot(sweepa[subc,2],sweepa[subc,"South_CLR"],pch=20,cex=.2)
abline(h=100)
plot(sweepa[subc,2],sweepa[subc,"Grand_CLR"],pch=20,cex=.2)
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