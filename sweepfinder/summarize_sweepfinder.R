library(magrittr)
library(dplyr)
library(stringr)

# wherever your data are
load("~/Dropbox/Public/sweepfinder.RData")

fai2 <- fai[]

xyscafs <- maxy[maxy[,2]=="xy",1]
sweepxy <- sweep[sweep[,1] %in% xyscafs,]

sweepa <- sweep[!(sweep[,1] %in% xyscafs),]

par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sweepa[,"North_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"Admix_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"South_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)
plot(sweepa[,"Grand_CLR"],pch=20,cex=.2,col=factor(sweepxy[,1]))
abline(h=500)