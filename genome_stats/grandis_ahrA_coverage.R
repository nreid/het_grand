
library(dplyr)
library(stringr)

# this line used to generate per sample depth for samples in grand.list
# NW_012234474.1:550000-750000 is the region encompassing the AHR deletions. 
# samtools depth -a -f ../../grandis_align/meta/GRAND.list -Q 30 -q 20 -r NW_012234474.1:550000-750000 >NW_012234474.1.depth.txt

# grand.list is a list of bam files. I read it in here to get sample names and population assignments. 

# get sample and population name vectors
gl <- read.table("~/projects/het_grand_data/grand.list",stringsAsFactors=FALSE)
gl <- unlist(gl) %>% str_extract(.,regex("BU.*[0-9]_[A-Z]+"))
pop <- gsub(".*_","",gl)

# load depth
dep <- read.table("~/projects/het_grand_data/NW_012234474.1.depth.txt.gz",stringsAsFactors=FALSE)
dep2 <- dep[,3:290]

# reorder samples by population tolerance

ord <- c(which(pop=="VB"),which(pop=="BB"),which(pop=="PB"),which(pop=="BNP"),which(pop=="SJSP"),which(pop=="SP"),which(pop=="GB"))

gl <- gl[ord]
pop <- pop[ord]
dep2 <- dep2[,ord]

# preliminary plots of depth for all grandis and per population
plot(dep[,2],rowSums(dep2),pch=20,cex=.2,ylim=c(0,300))
	
plot(dep[,2],rowSums(dep2[,pop=="BB"]),pch=20,cex=.2,ylim=c(0,40))

plot(dep[,2],rowSums(dep2[,pop=="VB"]),pch=20,cex=.2,ylim=c(0,50))
	plot(dep[,2],rowSums(dep2[,pop=="VB"]),pch=20,cex=.2,ylim=c(0,50),xlim=c(630000,640000))

plot(dep[,2],rowSums(dep2[,grep("VB|BB",pop)]),pch=20,cex=.2,ylim=c(0,40))

plot(dep[,2],rowSums(dep2[,pop=="PB"]),pch=20,cex=.2,ylim=c(0,50))

plot(dep[,2],rowSums(dep2[,pop=="BNP"]),pch=20,cex=.2,ylim=c(0,50))

plot(dep[,2],rowSums(dep2[,pop=="GB"]),pch=20,cex=.2,ylim=c(0,50))

plot(dep[,2],rowSums(dep2[,pop=="SP"]),pch=20,cex=.2,ylim=c(0,50))


# plot log ratios of population-level coverage
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
# VB/GB
plot(dep[,2],log(rowSums(dep2[,pop=="VB"])/rowSums(dep2[,pop=="GB"]),2),pch=20,cex=.2,ylim=c(-5,5))
# BB/GB
plot(dep[,2],log(rowSums(dep2[,pop=="BB"])/rowSums(dep2[,pop=="GB"]),2),pch=20,cex=.2,ylim=c(-5,5))
# PB/GB
plot(dep[,2],log(rowSums(dep2[,pop=="PB"])/rowSums(dep2[,pop=="GB"]),2),pch=20,cex=.2,ylim=c(-5,5))
# GB/SP
plot(dep[,2],log(rowSums(dep2[,pop=="GB"])/rowSums(dep2[,pop=="SP"]),2),pch=20,cex=.2,ylim=c(-5,5))

# keep: sites with "reasonable" coverage in ref pops SP and GB
# del: a subset of the sites in the deletion region
# ref: a subset of the sites outside the region


keep <- rowSums(dep2[,grep("GB|SP",pop)]) < 90 & rowSums(dep2[,grep("GB|SP",pop)]) > 40

del <- ((dep[,2] > 640000 & dep[,2] < 709000)) & keep
ref <- (((dep[,2] > 550000 & dep[,2] < 630000)) | ((dep[,2] > 710000 & dep[,2] < 750000))) & keep

# plot VB/GB ratio, sum(GB+SP), VB, and GB. 
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(dep[,2],log(rowSums(dep2[,pop=="VB"])/rowSums(dep2[,pop=="GB"]),2),pch=20,cex=.2,ylim=c(-6,6),col=del+(2*ref)+1)
abline(v=seq(550000,750000,10000),col="red")
plot(dep[,2],rowSums(dep2[,grep("GB|SP",pop)]),pch=20,cex=.2,ylim=c(0,110),col=del+(2*ref)+1)
abline(v=seq(550000,750000,10000),col="red")
plot(dep[,2],rowSums(dep2[,pop=="VB"]),pch=20,cex=.2,ylim=c(0,80),col=del+(2*ref)+1)
abline(v=seq(550000,750000,10000),col="red")
plot(dep[,2],rowSums(dep2[,pop=="GB"]),pch=20,cex=.2,ylim=c(0,80),col=del+(2*ref)+1)
abline(v=seq(550000,750000,10000),col="red")


# plot individual genotypes. first sorted by relative coverage
par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
(colMeans(dep2[del,]) / colMeans(dep2[ref,])) %>% sort() %>% plot(.,pch=20,cex=0.75)
abline(h=c(0,0.25,0.5,0.75,1))
# then sorted by population, tolerant to sensitive
(colMeans(dep2[del,]) / colMeans(dep2[ref,])) %>% plot(.,pch=20,cex=0.75,col=factor(pop))
abline(h=c(0,0.25,0.5,0.75,1))
# then plot average coverage per individual in case it's a confounder. 
colMeans(dep2) %>% plot(.,pch=20,cex=0.75,col=factor(pop))

# assign individual genotypes
# given above stats, use 0.25 and 0.7 as boundaries between genotypes. 

gtscore <- (colMeans(dep2[del,]) / colMeans(dep2[ref,]))

gts <- data.frame(sample=gl,pop=pop,genotype=NA,score=gtscore)

gts[gts[,4] < 0.25,3] <- 2
gts[gts[,4] > 0.25,3] <- 1
gts[gts[,4] > 0.7,3] <- 0

group_by(gts,pop) %>% summarize(.,mean(genotype)/2)

table(gts[,2:3])



(rowSums(dep2[,gts$genotype==2])/median(rowSums(dep2[1:75000,gts$genotype==2]))) %>% plot(x=dep[,2],y=.,ylim=c(0,3),pch=20,cex=.2)
(rowSums(dep2[,gts$genotype==0])/median(rowSums(dep2[1:75000,gts$genotype==0]))) %>% points(x=dep[,2],y=.,ylim=c(0,3),pch=20,cex=.2,col=rgb(1,0,0,.1))




