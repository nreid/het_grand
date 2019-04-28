library(stringr)
library(dplyr)
library(ape)
library(phytools)
library(magrittr)
library(viridis)
library(beeswarm)
library(tidyverse)

# some useful functions
pbs <- function(t1,t2,c12){
	
	t1 <- -log(1-t1)
	t2 <- -log(1-t2)
	c12 <- -log(1-c12)

	stat <- (t1 + t2 - c12)/2
	return(stat)
	}

pbspop <- function(tar,ca,cb,s,mat){

	tarca <- sort(c(tar,ca)) %>% paste(.,collapse=".")
	tarcb <- sort(c(tar,cb)) %>% paste(.,collapse=".")
	cacb <- sort(c(ca,cb)) %>% paste(.,collapse=".")
	t1 <- 1-(mat[s,tar]+mat[s,ca])/2/mat[s,tarca]
	t2 <- 1-(mat[s,tar]+mat[s,cb])/2/mat[s,tarcb]
	c12 <- 1-(mat[s,ca]+mat[s,cb])/2/mat[s,cacb]
	pbs(t1=t1,t2=t2,c12=c12)

}


smother <- function(x,winsize){
	
	vec <- 1:length(x)
	start <- vec - winsize; start[start < 1] <- 1
	end <- vec + winsize; end[end > length(x)] <- length(x)
	win <- cbind(start,end)
	out <- apply(win, MAR=1,FUN=function(z){sum(x[z[1]:z[2]],na.rm=TRUE)})
	return(out)

}


smotherM <- function(x,winsize){
	
	vec <- 1:length(x)
	start <- vec - winsize; start[start < 1] <- 1
	end <- vec + winsize; end[end > length(x)] <- length(x)
	win <- cbind(start,end)
	out <- apply(win, MAR=1,FUN=function(z){mean(x[z[1]:z[2]],na.rm=TRUE)})
	return(out)

}

mergewin<-function(win, stat, qu, tails=TRUE, buff=5000){

	ord <- order(win[,1],win[,2],win[,3])
	win <- win[ord,]
	stat <- stat[ord]

	stat<-as.numeric(stat)
	stat[is.na(stat)]<-0
	if(tails=="lesser"){
		qu<-qu[1]
		stat<-stat-qu
		stat[stat>0]<-0	
		}
	if(tails=="greater"){
		qu<-qu[length(qu)]
		stat<-stat-qu
		stat[stat<0]<-0
		}
	if(tails=="both"){
		if(length(qu)==1){stop("two quantiles must be provided for option 'both'")}
		stat[stat>qu[1]&stat<qu[2]]<-0
		stat[stat<qu[1]]<-stat[stat<qu[1]]-qu[1]
		stat[stat>qu[2]]<-stat[stat>qu[2]]-qu[2]
		}
	
	#chrom,start,end,length,auc,highpoint,count
	out<-cbind(win[1,],win[1,3]-win[1,2],stat[1],stat[1],1)
	names(out)<-c("chrom","start","end","length","AUC","highpoint","count")
	for(i in 2:length(win[,1])){
		
		ind<-length(out[,1])
		if(out[ind,1]==win[i,1]&(out[ind,3]+buff)>=win[i,2]){
			out[ind,2]<-min(win[i,2],out[ind,2])
			out[ind,3]<-max(win[i,3],out[ind,3])
			out[ind,4]<-out[ind,3]-out[ind,2]
			out[ind,5]<-out[ind,5]+stat[i]
			if(abs(out[ind,6])<abs(stat[i])){out[ind,6]<-stat[i]}
			if(stat[i]!=0){out[ind,7]<-out[ind,7]+1}
			}else{
					out<-rbind(out,setNames(cbind(win[i,],win[i,3]-win[i,2],stat[i],stat[i],1),names(out)))
					}		
		}
	return(out)
	}
	

# functions for extracting pi and dxy for sex chromosomes, given pi and dxy for sexes
xypi <- function(mf,ff){ 2*mf - ff }

yypi <- function(mm,mf1,mf2,ff){

	xy1 <- xypi(mf1,ff)
	xy2 <- xypi(mf2,ff)

	4 * ( mm - 0.25*xy1 - 0.25*xy2 - 0.25*ff)

	}

# equivalent to above function, just simplified. 
yypi <- function(mm,mf1,mf2,ff){

	4 * ( mm - 0.5*mf1 - 0.5*mf2 + 0.25*ff)

	}


getchrpi <- function(p1,p2,xy1,xy2,mat=fst2){

	# this function uses the sex-stratified pi/dxy matrix and pulls out pi/dxy for x and y chromosomes
	# feed it population(s) and sex chromosome(s)

	# first set of statements gets values within populations
	if(p1==p2){

		if(paste(sort(c(xy1,xy2)),collapse="")=="XY"){
			m1 <- paste(p1,"M",sep="")
			f1 <- paste(p1,"F",sep="")
			mf1 <- which(grepl(m1,colnames(mat)) & grepl(f1,colnames(mat)))
			ff1 <- which(colnames(mat)==f1)
			out <- xypi(mat[,mf1],mat[,ff1])
			return(out)
		}
		if(paste(sort(c(xy1,xy2)),collapse="")=="YY"){
			m1 <- paste(p1,"M",sep="")
			f1 <- paste(p1,"F",sep="")
			mm1 <- which(colnames(mat)==m1)
			mf1 <- which(grepl(m1,colnames(mat)) & grepl(f1,colnames(mat)))
			ff1 <- which(colnames(mat)==f1)
			out <- yypi(mat[,mm1],mat[,mf1],mat[,mf1],mat[,ff1])
			return(out)
		}
		if(paste(sort(c(xy1,xy2)),collapse="")=="XX"){	
			f1 <- paste(p1,"F",sep="")
			ff1 <- which(colnames(mat)==f1)
			out <- mat[,ff1]
			return(out)
		}

	}

	# second set of statements gets values between populations
	if(p1!=p2){

		if(paste((c(xy1,xy2)),collapse="")=="XY"){

			z1 <- paste(p1,"F",sep="")
			z2 <- paste(p2,"F",sep="")
			z3 <- paste(p2,"M",sep="")
			mf1 <- which(grepl(z1,colnames(mat)) & grepl(z3,colnames(mat)))
			ff1 <- which(grepl(z1,colnames(mat)) & grepl(z2,colnames(mat)))
			out <- xypi(mat[,mf1],mat[,ff1])			
			return(out)
		}

		if(paste((c(xy1,xy2)),collapse="")=="YX"){

			z1 <- paste(p1,"M",sep="")
			z2 <- paste(p1,"F",sep="")
			z3 <- paste(p2,"F",sep="")
			mf1 <- which(grepl(z1,colnames(mat)) & grepl(z3,colnames(mat)))
			ff1 <- which(grepl(z2,colnames(mat)) & grepl(z3,colnames(mat)))
			out <- xypi(mat[,mf1],mat[,ff1])
			return(out)
		}


		if(paste((c(xy1,xy2)),collapse="")=="YY"){
			m1 <- paste(p1,"M",sep="")
			m2 <- paste(p2,"M",sep="")
			f1 <- paste(p1,"F",sep="")
			f2 <- paste(p2,"F",sep="")
			mm1 <- which(grepl(m1,colnames(mat)) & grepl(m2,colnames(mat)))
			mf1 <- which(grepl(m1,colnames(mat)) & grepl(f2,colnames(mat)))
			mf2 <- which(grepl(m2,colnames(mat)) & grepl(f1,colnames(mat)))
			ff1 <- which(grepl(f1,colnames(mat)) & grepl(f2,colnames(mat)))
			out <- yypi(mat[,mm1],mat[,mf1],mat[,mf2],mat[,ff1])
			return(out)
		}

		if(paste((c(xy1,xy2)),collapse="")=="XX"){	
			f1 <- paste(p1,"F",sep="")
			f2 <- paste(p2,"F",sep="")
			ff1 <- which(grepl(f1,colnames(mat)) & grepl(f2,colnames(mat)))
			out <- mat[,ff1]
			return(out)
		}



	}

}




makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
	}

options(scipen=99)


# start reading in and processing data

# population segments
segs <- c(
	"northM","northF",
	"southM","southF",
	"admixedM","admixedF",
	"grandM","grandF") 

# old scaffold to NCBI scaffold mappings
om <- read.table("~/projects/het_grand/old_scaffold_mappings.txt",stringsAsFactors=FALSE)
rownames(om) <- om[,1]
om2 <- om[-10180,]
rownames(om2) <- om2[,2]

# window liftover from NCBI scaffold to Miller chromosomes
lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

# table of how many sites per window were analyzed
val <- read.table("~/projects/het_grand_data/popsites.1kb_win.bed.gz",stringsAsFactors=FALSE)
val <- val[lord, ]
val2 <- val
val2[,4] <- smother(val2[,4],10)
subv <- val[,4] > 200

# read in dxy,pi table
fst <- read.table("~/projects/het_grand_data/fst_dxy_sex.gz",stringsAsFactors=FALSE)
pnames <- c("scaf","start","end",segs,combn(segs,2) %>% apply(.,MAR=2,FUN=paste,collapse="."))
colnames(fst) <- pnames

for(i in 4:39){ fst[,i] <- as.numeric(fst[,i])}
for(i in 4:39){ fst[which(fst[,i] > 200),i] <- NA}


fst <- fst[lord,]
fst2 <- fst
for(i in 4:39){ fst2[,i] <- smother(fst2[,i],10);print(i)}

# read in reference genome base counts per window
acgt <- read.table("ACGT_1kb.bed",stringsAsFactors=FALSE)
acgt <- acgt[lord,]
acgt2 <- acgt
for(i in 4:12){acgt2[,i] <- smother(acgt2[,i],10);print(i)}

# read in sex linked interval bedfile
snvint <- read.table("heteroclitus_SNVsex_intervals.bed",stringsAsFactors=FALSE)

# find windows overlapping above windows
slwin <- rep(FALSE,dim(fst)[1])
for(i in 1:dim(snvint)[1]){

	slwin <- slwin | ((snvint[i,2] < fst[,3] & snvint[i,3] > fst[,2]) & snvint[i,1]==fst[,1])

}

# read in sex linked variants, table of sex assignments of individuals
vcf <- read.table("sexSNVs.ap.vcf.gz",stringsAsFactors=FALSE)
h <- scan(pipe("gzcat sexSNVs.vcf.gz | head -n 20000 | grep CHROM"),what="character",sep="\t")
h[1] <- "CHROM"
colnames(vcf) <- h

# get only biallelic variants
vcf <- vcf[!grepl(",",vcf[,5]),]

# table of vcf genotypes only
gts <- as.matrix(vcf[,-(1:9)])
class(gts) <- "numeric"

# get sex table
sexes <- read.table("../het_grand/all_sexes.txt",stringsAsFactors=FALSE)
rownames(sexes) <- sexes[,1]
# three individuals mislabeled, fixed here. 
	# mislabeled because I changed them from the
	# morphological designation on the basis of one scaffold
	# that scaffold is not always part of the SDR. 
	# only impact one population, Flax Pond.
sexes[c("F-20", "F-2", "F-47"),3] <- "M"


# population designations
pop <- c("north","south","admixed","grand")
popxy <- c("northX","northY","southX","southY","admixedX","admixedY","grandX","grandY")

# make a distance matrix of average x,y dxy values
xymat <- matrix(nrow=8,ncol=8)
colnames(xymat) <- popxy
rownames(xymat) <- popxy

for(i in pop){

	for(j in pop){

		for(k in c("X","Y")){

			for(l in c("X","Y")){

				cc <- paste(i,k,sep="")
				rr <- paste(j,l,sep="")
				xymat[rr,cc] <- sum(getchrpi(i,j,k,l,mat=fst[slwin,]),na.rm=TRUE)/sum(as.numeric(val[slwin,4]),na.rm=TRUE)

			}
		}
	}
}

# matrix of X-Y fst values for SDR
xyfstmat <- matrix(nrow=8,ncol=8,data=NA)
colnames(xyfstmat) <- popxy
rownames(xyfstmat) <- popxy

for(i in pop){

	for(j in pop){

		for(k in c("X","Y")){

			for(l in c("X","Y")){

				cc <- paste(i,k,sep="")
				rr <- paste(j,l,sep="")
				xyfstmat[rr,cc] <- 1 - ((xymat[cc,cc] + xymat[rr,rr])/(2*xymat[cc,rr]))

			}
		}
	}
}


# make a table of windowed x,y dxy values
sdnames <- c(colnames(fst2)[1:3],popxy,combn(popxy,2) %>% apply(.,MAR=2,FUN=paste,collapse="."))
sexdxy <- fst
sexdxy2 <- fst2
colnames(sexdxy2) <- sdnames
colnames(sexdxy) <- sdnames

for(i in pop){

	for(j in pop){

		for(k in c("X","Y")){

			for(l in c("X","Y")){

				cc <- paste(i,k,sep="")
				rr <- paste(j,l,sep="")

				if(cc==rr){ 
					sexdxy2[,rr] <- getchrpi(i,j,k,l,mat=fst2) 
					sexdxy[,rr] <- getchrpi(i,j,k,l,mat=fst)
				}
				if(cc!=rr){ 
					cn <- paste(cc,rr,sep=".")
					if(cn %in% colnames(sexdxy2)){ 
						sexdxy2[,cn] <- getchrpi(i,j,k,l,mat=fst2) 
						sexdxy[,cn] <- getchrpi(i,j,k,l,mat=fst)
					}

				}

			}
		}
	}
}

# make a table of windowed x,y Fst values
sdnames <- c(colnames(fst2)[1:3],combn(popxy,2) %>% apply(.,MAR=2,FUN=paste,collapse="."))
sexfst <- fst[,1:31]
sexfst2 <- fst2[,1:31]
colnames(sexfst2) <- sdnames

k <- 4
for(i in 1:7){
	for(j in (i+1):8){
		print(c(popxy[i],popxy[j]))

		p1 <- sexdxy2[,popxy[i]]
		p2 <- sexdxy2[,popxy[j]]
		p12 <- sexdxy2[,paste(popxy[i],popxy[j],sep=".")]
		sexfst2[,k] <- 1 - (p1 + p2)/2/p12

		p1 <- sexdxy[,popxy[i]]
		p2 <- sexdxy[,popxy[j]]
		p12 <- sexdxy[,paste(popxy[i],popxy[j],sep=".")]
		sexfst[,k] <- 1 - (p1 + p2)/2/p12
		k <- k+1
	}
}


#################
# START PLOTS ETC
#################

###########
# make plot of Dxy/Pi/Fst
# official figure, probably maybe?
##########

xymat2 <- xymat[c(2,6,4,1,5,3,7),c(2,6,4,1,5,3,7)]
xypimat <- xymat2

xymat2[lower.tri(xymat2,diag=TRUE)] <- NA

xypimat[upper.tri(xypimat)] <- NA
xypimat[lower.tri(xypimat)] <- NA

xyfstmat2 <- xyfstmat[c(2,6,4,1,5,3,7),c(2,6,4,1,5,3,7)]
xyfstmat2[upper.tri(xyfstmat2,diag=TRUE)] <- NA

pdf(width=16, height=4, file="Fig_1_dxy_pi_fst_tree.pdf")

par(mfrow=c(1,4))

image(xymat2,col=viridis(20),
	x=1:7,y=1:7,xaxt='n',yaxt='n',
	ylab=NA,xlab=NA,bty='n')
tcol <- "black"
for(i in 1:7){
	for(j in 1:7){
		try(if(xymat2[i,j]<0.02){print("yes");tcol <- "white"})
		text(i,j,round(xymat2[i,j],digits=3),col=tcol)
		tcol <- "black"
		}
	}

image(xypimat,col=viridis(20),
	x=1:7,y=1:7,xaxt='n',yaxt='n',
	ylab=NA,xlab=NA,bty='n')

tcol <- "black"
for(i in 1:7){
	for(j in 1:7){
		try(if(xypimat[i,j]<0.012){print("yes");tcol <- "white"})
		text(i,j,round(xypimat[i,j],digits=3),col=tcol)
		tcol <- "black"
		}
	}

image(xyfstmat2,col=viridis(20),
	x=1:7,y=1:7,xaxt='n',yaxt='n',
	ylab=NA,xlab=NA,bty='n')

tcol <- "black"
for(i in 1:7){
	for(j in 1:7){
		try(if(xyfstmat2[i,j]<0.3){print("yes");tcol <- "white"})
		text(i,j,round(xyfstmat2[i,j],digits=3),col=tcol)
		tcol <- "black"
		}
	}


as.dist(xymat[-8,-8]) %>% nj() %>% midpoint.root() %>% plot()
add.scale.bar(length=0.005)

dev.off()

########################

# some window subsets based on number of bases with coverage

# all chromosome 5
chr5 <- val2[,4] > 7000 & grepl("chr5",lift[,1])

# slwin is SDR windows
slwin2 <- slwin & val2[,4] > 7000

# non chr5 windows
subw <- val2[,4] > 7000 & grepl("chr",lift[,1]) & !grepl("chr5",lift[,1])

# subw <- val2[,4] > 7000 & grepl("chr",lift[,1])

# subw <- val2[,4] > 7000 & grepl("NW_012224574.1",val[,1])

#################
# simple divergence time measures:
	# T = (Dxy - piA) / 2
#################

# genetic diversity divergence of everything except chr5
div <- (colSums((sexdxy[subw,-c(1:3)]),na.rm=TRUE)/sum(as.numeric(val[subw,4]))) %>% data.frame()

# genetic diversity divergence of everything except chr5
div2 <- (colSums((fst[subw,-c(1:3)]),na.rm=TRUE)/sum(as.numeric(val[subw,4]))) %>% data.frame()

# genetic diversity divergence of SDR
divSDR <- (colSums((sexdxy[slwin,-c(1:3)]),na.rm=TRUE)/sum(as.numeric(val[slwin,4]))) %>% data.frame()

# divergence time = (dxy - pi_anc) / 2
# X-Y divergence time
	# assume South represents ancestral diversity
	(divSDR["southX.southY",1] - divSDR["southX",1])/2
	# 			0.003002158 (1%,99%: 0.002956652 0.003053281)
	# assume North represents ancestral diversity
	(divSDR["northX.northY",1] - divSDR["northX",1])/2
	# 			0.00449272
	# bootstrap replicate
	bootr <- c()
	for(i in 1:1000){
		sam <- sample(which(slwin),replace=TRUE)
		bootr <- c(bootr, sum(sexdxy[sam,"southX.southY"]-sexdxy[sam,"southX"],na.rm=TRUE)/sum(val[sam,4])/2)
		}

# north-south divergence time
	# assume South represents ancestral diversity
	(div["northY.southY",1] - div["southY",1])/2
	# 			0.0007558629
	# assume North represents ancestral diversity
	(div["northY.southY",1] - div["northY",1])/2
	# 			0.003210888
	# bootstrap replicates
	bootr <- c()
	for(i in 1:1000){
		sam <- sample(which(subw),replace=TRUE)
		bootr <- c(bootr, sum(sexdxy[sam,"northY.southY"]-sexdxy[sam,"southY"],na.rm=TRUE)/sum(val[sam,4])/2)
		}

# north-admixed divergence time, negative b/c ??
	(div["northY.admixedY",1] - div["southX",1])/2
	#			-0.0003114985

# admixed-south divergence time, negative b/c ??
	(div["southX.admixedY",1] - div["southX",1])/2
	#			-0.0003114985

# het-grand divergence time
	# assume South represents ancestral diversity
	(div["northX.grandX",1] - div["southX",1])/2
	# 			0.008069279
	# assume North represents ancestral diversity
	(div["northX.grandX",1] - div["northX",1])/2
	# 			0.01045042	
	# assume Grandis represents ancestral diversity
	(div["southX.grandX",1] - div["grandX",1])/2
	# 			0.0134297

# mean time to coalescence of SDR's by pop
	# North:	0.000998055
	# Admixed:	0.000531351
	# South:	0.000275207
	# NY-SY: 

# X/A and Y/A diversity ratios

	# North:
		# X/A:	1.024976
		# Y/A:	0.1479289

	# Admixed:
		# X/A:	0.9516145
		# Y/A:	0.06406621

	# South:
		# X/A:	0.8716144
		# Y/A:	0.03014976

#######################
# plot of genetic diversity on chromosomes, SDR, PAR, etc
#######################

# all chromosome 5
chr5 <- val2[,4] > 7000 & grepl("chr5",lift[,1])

# slwin is SDR windows

# non chr5 windows
subw <- val2[,4] > 7000 & grepl("chr",lift[,1]) & !grepl("chr5",lift[,1])

# get diversity by population per chromosome/SDR/PAR
chrdiv <- group_by(sexdxy,CHROM=lift[,1]) %>% 
	summarize(., 
		northX=sum(northX,na.rm=TRUE),
		northY=sum(northY,na.rm=TRUE),
		admixedX=sum(admixedX,na.rm=TRUE),
		admixedY=sum(admixedY,na.rm=TRUE),
		southX=sum(southX,na.rm=TRUE),
		southY=sum(southY,na.rm=TRUE),
		grandX=sum(grandX,na.rm=TRUE),
		grandY=sum(grandY,na.rm=TRUE)
		) %>% 
	data.frame() %>% head(.,n=24)

# paste SDR and PAR onto table
chrdiv[25,1] <- "SDR"
chrdiv[25,2:9]  <- colSums(sexdxy[slwin,c(4,5,8,9,6,7,10,11)],na.rm=TRUE)
chrdiv[26,1] <- "PAR"
chrdiv[26,2:9]  <- colSums(sexdxy[chr5 & !slwin,c(4,5,8,9,6,7,10,11)],na.rm=TRUE)

# get # of sites to scale by
chrval <- group_by(val,CHROM=lift[,1]) %>% 
	summarize(.,
		sites=sum(V4)) %>% 
	data.frame() %>% head(.,n=24)

chrval[25,1] <- "SDR"
chrval[25,2] <- sum(val[slwin,4])
chrval[26,1] <- "PAR"
chrval[26,2] <- sum(val[chr5 & !slwin,4])

# scale diversity table by sites
for(i in 2:9){

	chrdiv[,i] <- chrdiv[,i]/chrval[,2]

}

# one column per population, necessary because all calculations made using X-Y comparisons, even on autosomes
chrdiv2 <- data.frame(

	north=c(chrdiv[,"northX"],chrdiv[25,3],chrdiv[26,3]),
	admixed=c(chrdiv[,"admixedX"],chrdiv[25,5],chrdiv[26,5]),
	south=c(chrdiv[,"southX"],chrdiv[25,7],chrdiv[26,7]),
	grand=c(chrdiv[,"grandX"],NA,NA)
	)

rownames(chrdiv2) <- c(chrdiv[1:24,1],"SDRx","PARx","SDRy","PARy")

# long format, exclude SDR,PAR
chrdiv3 <- gather(chrdiv2[1:24,][-20,])

# neutral X/Y expectations:

yy <- colSums(sexdxy[subw,c(4,8,6)],na.rm=TRUE)/sum(val[subw,4]) / 4
xx <- colSums(sexdxy[subw,c(4,8,6)],na.rm=TRUE)/sum(val[subw,4]) / 4 * 3

# plot autosomes, add in SDR, PAR
pdf(width=7,height=7,file="Fig_2_XYAdiversity.pdf")

beeswarm(chrdiv3[,2] ~ factor(chrdiv3[,1],levels=c("north","admixed","south","grand")),
	ylim=c(0,0.021),
	pch=20,col="gray",cex=1.5,
	xlab="Population",
	ylab="Genetic Diversity",
	labels=c("North","Admixed","South","F. grandis"))
points(
	x=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)),
	y=unlist(chrdiv2[25:28,]),
	pch=rep(c(21,20,21,20),4),
	col=rep(c("black","black","red","red"),4),
	lwd=2.5,
	cex=rep(c(1,1.5,1,1.5),4))
legend(y=0.020,x=3.7,
	legend=c("autosome","SDR (X)","PAR (X)","SDR-Y","PAR-Y"),
	pch=c(20,21,20,21,20),
	pt.cex=c(1.5,1,1.5,1,1.5),
	pt.lwd=2.5,
	col=c("gray","black","black","red","red")
	)



# lines(x=c(0.85,1.15),y=c(yy[1],yy[1]),lwd=3,col="red")
# lines(x=c(1.85,2.15),y=c(yy[2],yy[2]),lwd=3,col="red")
# lines(x=c(2.85,3.15),y=c(yy[3],yy[3]),lwd=3,col="red")

# lines(x=c(0.85,1.15),y=c(xx[1],xx[1]),lwd=3,col="black")
# lines(x=c(1.85,2.15),y=c(xx[2],xx[2]),lwd=3,col="black")
# lines(x=c(2.85,3.15),y=c(xx[3],xx[3]),lwd=3,col="black")

# instead of lines, use bxplot to get a distribution...
bxplot((chrdiv3[,2]/4 * 3) ~ factor(chrdiv3[,1],levels=c("north","admixed","south")),add=TRUE,probs=c(0.0,0.5,1),width=0.25)
bxplot((chrdiv3[,2]/4) ~ factor(chrdiv3[,1],levels=c("north","admixed","south")),add=TRUE,probs=c(0.0,0.5,1),col="red",width=0.25)

dev.off()

##########################


#######################
# plot of dxy to grandis for autosomes, SDR, PAR, etc
#######################

# all chromosome 5
chr5 <- val2[,4] > 7000 & grepl("chr5",lift[,1])

# slwin is SDR windows

# non chr5 windows
subw <- val2[,4] > 7000 & grepl("chr",lift[,1]) & !grepl("chr5",lift[,1])

# get diversity by population per chromosome/SDR/PAR
chrdxy <- group_by(sexdxy,CHROM=lift[,1]) %>% 
	summarize(., 
		northX.grandX=sum(northX.grandX,na.rm=TRUE),
		northY.grandX=sum(northY.grandX,na.rm=TRUE),
		admixedX.grandX=sum(admixedX.grandX,na.rm=TRUE),
		admixedY.grandX=sum(admixedY.grandX,na.rm=TRUE),
		southX.grandX=sum(southX.grandX,na.rm=TRUE),
		southY.grandX=sum(southY.grandX,na.rm=TRUE)
		) %>% 
	data.frame() %>% head(.,n=24)

# paste SDR and PAR onto table
chrdxy[25,1] <- "SDR"
chrdxy[25,2:7]  <- colSums(sexdxy[slwin,which(colnames(sexdxy) %in% colnames(chrdxy)[-1])[c(1,2,5,6,3,4)]],na.rm=TRUE)
chrdxy[26,1] <- "PAR"
chrdxy[26,2:7]  <- colSums(sexdxy[chr5 & !slwin,which(colnames(sexdxy) %in% colnames(chrdxy)[-1])[c(1,2,5,6,3,4)]],na.rm=TRUE)

# get # of sites to scale by
chrval <- group_by(val,CHROM=lift[,1]) %>% 
	summarize(.,
		sites=sum(V4)) %>% 
	data.frame() %>% head(.,n=24)

chrval[25,1] <- "SDR"
chrval[25,2] <- sum(val[slwin,4])
chrval[26,1] <- "PAR"
chrval[26,2] <- sum(val[chr5 & !slwin,4])

# scale diversity table by sites
for(i in 2:7){

	chrdxy[,i] <- chrdxy[,i]/chrval[,2]

}

# one column per population, necessary because all calculations made using X-Y comparisons, even on autosomes
chrdxy2 <- data.frame(

	north=c(chrdxy[,"northX.grandX"],chrdxy[25,3],chrdxy[26,3]),
	admixed=c(chrdxy[,"admixedX.grandX"],chrdxy[25,5],chrdxy[26,5]),
	south=c(chrdxy[,"southX.grandX"],chrdxy[25,7],chrdxy[26,7])

	)

rownames(chrdxy2) <- c(chrdxy[1:24,1],"SDRx","PARx","SDRy","PARy")

# long format, exclude SDR,PAR
chrdxy3 <- gather(chrdxy2[1:24,][-20,])


# plot autosomes, add in SDR, PAR
beeswarm(chrdxy3[,2] ~ factor(chrdxy3[,1],levels=c("north","admixed","south")),
	ylim=c(0.03,0.04),
	pch=20,col="gray",cex=1.5,
	xlab="Population",
	ylab="Dxy heterolitus-grandis",
	labels=c("North","Admixed","South"))
points(
	x=c(rep(1,4),rep(2,4),rep(3,4)),
	y=unlist(chrdxy2[25:28,]),
	pch=rep(c(21,20,21,20),3),
	col=rep(c("black","black","red","red"),3),
	lwd=2.5,
	cex=rep(c(1,1.5,1,1.5),3))
legend(y=0.040,x=3,
	legend=c("autosome","SDR (X)","PAR (X)","SDR-Y","PAR-Y"),
	pch=c(20,21,20,21,20),
	pt.cex=c(1.5,1,1.5,1,1.5),
	pt.lwd=2.5,
	col=c("gray","black","black","red","red")
	)

lines(x=c(0.85,1.15),y=c(yy[1],yy[1]),lwd=3,col="red")
lines(x=c(1.85,2.15),y=c(yy[2],yy[2]),lwd=3,col="red")
lines(x=c(2.85,3.15),y=c(yy[3],yy[3]),lwd=3,col="red")

lines(x=c(0.85,1.15),y=c(xx[1],xx[1]),lwd=3,col="black")
lines(x=c(1.85,2.15),y=c(xx[2],xx[2]),lwd=3,col="black")
lines(x=c(2.85,3.15),y=c(xx[3],xx[3]),lwd=3,col="black")


# chromosome 23 
	# is the shortest (in this map)
	# has the highest genetic diversity
	# has the largest dxy to grandis
# the X-SDR
	# has the smallest dxy to grandis

# there's a systematic problem with estimating diversity
	# the more diverse you are, the more you underestimate diversity and divergence
	# y-SDR has the most stable estimate of divergence b/c it has very low diversity and doesn't vary much among pops

# for future reference:
	# the SDR contains 595 of 25692 genes, or 2.3% of the annotated NCBI genes. 
	# ~1.05mb of 25.7mb, or 4.1%, of the SDR is coding. 


##########################



##########################

#construct stupid population tree based on divergence times
ptree <- read.tree(text="(((north:215960.8,south:215960.8):445444.6,hetY:661405.4):1644103,grandis:2305508);")
# replace 
ptree$edge.length <- ptree$edge.length * 3.5e-9
plot(ptree)

plot(sexfst2[chr5,"northX.northY"],pch=20,cex=.2,col=slwin[chr5]+1)

plot((sexdxy2[chr5,"southX.southY"]-sexdxy2[chr5,"southX"])/val2[chr5,4],pch=20,cex=.2,col=slwin[chr5]+1)
points((sexdxy2[chr5,"southX.grandX"]-sexdxy2[chr5,"southX"])/val2[chr5,4],pch=20,cex=.2,col="blue")


################
# Score individuals and windows based on maleness
# Try to identify large-scale variation in SDR/recombination events
################

# check some weird individuals
# weird individuals: F-20  F-2 F-47 SH-16


# score heteroclitus individuals by how many strongly male-linked markers they have. 
hetf <- which(sexes[colnames(vcf)[-(1:9)],3]=="F" & !grepl("BU",colnames(vcf)[-(1:9)]))
hetm <- which(sexes[colnames(vcf)[-(1:9)],3]=="M" & !grepl("BU",colnames(vcf)[-(1:9)]))
hetff <- rowMeans(gts[,hetf],na.rm=TRUE)
hetfm <- rowMeans(gts[,hetm],na.rm=TRUE)

# select markers where alt allele is male-linked
sexsnv <- hetff < 0.1 & hetfm > 0.4

# individuals with less than 80% missing data
subi <- colMeans(is.na(gts)) < 0.8

# plot allele frequencies by sex
plot(hetff,hetfm,pch=20,cex=.2,col=rgb(0,0,0,.2))

# score individuals by average number of male-linked markers
sexscore <- colMeans(gts[sexsnv,],na.rm=TRUE)

# plot scores
sexscore %>% plot(.,pch=20,col=factor(sexes[colnames(gts),3]))
# try this as a beeswarm plot. booo. 
beeswarm(sexscore ~ factor(sexes[colnames(gts),3]) + factor(sexes[colnames(gts),2]))

# calculate a "Y" score in windows across SDR
sexsmooth <- gts[sexsnv,]
for(i in 1:dim(gts)[2]){sexsmooth[,i] <- smotherM(x=sexsmooth[,i],winsize=100); print(i)}

# three individuals mislabeled (fixed above) b/c of rare recombination?
	# F-20 	F-2 	F-47 
	# they don't have the full SDR


# figure out which individuals have variant of SDR
	# exclude males with high missing data. 
	# count number of windows with score falling below 0.2
	# rank individuals. 

lowvarM <- hetm[colMeans(is.na(gts[,hetm])) < 0.75]
lowvarF <- hetf[colMeans(is.na(gts[,hetf])) < 0.75]
colSums(sexsmooth[,lowvarM] < 0.2,na.rm=TRUE) %>% sort() %>% data.frame()

lowin <- rowSums(sexsmooth[,lowvarM] < 0.2 ,na.rm=TRUE)


par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))

plot(sexsmooth[,"F-20"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]))
plot(sexsmooth[,"F-2"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]))
plot(sexsmooth[,"F-47"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]))

plot(sexsmooth[,"KC-20"],pch=20,col=factor(vcf[sexsnv,1]),ylim=c(0,1), cex=(lowin > 3)*1 + 0.1)



######
# Example plots of male, recombinant male and female
# Supp fig 2 or 3 for paper
######
fup <- apply(sexsmooth[,lowvarF],MAR=1,FUN=quantile,prob=0.95,na.rm=TRUE)

pdf(width=8, height=3, file="SDR_variant.pdf")
# png(width=12, height=4, file="SDR_variant.png", units="in",res=300)

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))

plot(sexsmooth[,"BP-20"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]),ylim=c(0,0.7),ylab="Male BP-20")
abline(v=which(!duplicated(vcf[sexsnv,1]))[-1],col=rgb(0,0,0,0.25))
plot(sexsmooth[,"F-20"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]),ylim=c(0,0.7),ylab="Male F-20")
# fup %>% points(.,type="l",col=rgb(0,0,0,.5))
abline(v=which(!duplicated(vcf[sexsnv,1]))[-1],col=rgb(0,0,0,0.25))
plot(sexsmooth[,"KC-37"],pch=20,cex=.1,col=factor(vcf[sexsnv,1]),ylim=c(0,0.7),ylab="Female KC-37",xlab="SNP index")
# fup %>% points(.,type="l",col=rgb(0,0,0,.5))
abline(v=which(!duplicated(vcf[sexsnv,1]))[-1],col=rgb(0,0,0,0.25))

dev.off()

###############
# mds of snps in "recombinant" Y regions for females + 2 recombinant Y males
###############
# vcf[,1]=="NW_012224869.1"
mds <- gts[grep("NW_012224869.1|NW_012234316.1|NW_012234324.1",vcf[,1]),lowvarF] %>% t() %>% dist() %>% cmdscale()
mds2 <- gts[grep("NW_012224869.1|NW_012234316.1|NW_012234324.1",vcf[,1]),c(colnames(sexsmooth)[lowvarF],"F-20","F-2")] %>% t() %>% dist() %>% cmdscale()
rytr <- gts[grep("NW_012224869.1|NW_012234316.1|NW_012234324.1",vcf[,1]),c(colnames(sexsmooth)[lowvarF],"F-20","F-2")] %>% t() %>% dist() %>% nj()

par(mfrow=c(1,2))
plot(mds,col=gsub("-.*","",rownames(mds)) %>% factor())
plot(mds2,col=gsub("-.*","",rownames(mds2)) %>% factor())

plot(rytr,"unrooted",cex=0.5)

# nothing particularly exciting results. 

#############


plot(sexsmooth[,1],pch=20,cex=.2,col="blue",ylim=c(0,1))
points(sexsmooth[,2],pch=20,cex=.2,col="black")
points(sexsmooth[,3],pch=20,cex=.2,col="green")




plot(sexsmooth,pch=20,cex=.2,col=factor(vcf[sexsnv,1]),ylim=c(0,1))

# how many SDR windows in F-20, KC-20 have X-like sequence
# mergewin<-function(win, stat, qu, tails=TRUE, buff=5000)

xsnv <- vcf[sexsnv,1:2][sexsmooth[,"KC-20"] < 0.2,]
xsnv <- data.frame(CHROM=xsnv[,1],START=xsnv[,2]-1,END=xsnv[,2],stringsAsFactors=FALSE)
lost <- mergewin(xsnv,stat=rep(0,dim(xsnv)[1]),qu=1,buff=30000)

plot(sexsmooth[,"F-20"],pch=20,col=factor(vcf[sexsnv,1]),ylim=c(0,1), cex=(lowin > 3)*1 + 0.1)
abline(v=which(vcf[sexsnv,2] %in% (lost[,2]+1)),lwd=2)
abline(v=which(vcf[sexsnv,2] %in% (lost[,3])),col="red")

# "female" scaffolds in weird individuals: "NW_012224869.1","NW_012234324.1"
colMeans(gts[sexsnv & grepl("NW_012224869.1|NW_012234324.1",vcf[,1]),],na.rm=TRUE) %>% plot(.,col=factor(sexes[colnames(gts),2]),pch=(sexes[colnames(gts),3]=="M") + 1)

colMeans(gts[sexsnv,],na.rm=TRUE) %>% plot(.,col=factor(sexes[colnames(gts),2]),pch=(sexes[colnames(gts),3]=="M") + 1)


mds <- t(gts[sexsnv & vcf[,1]=="NW_012224869.1",subi & !grepl("BU", colnames(gts))]) %>% dist() %>% cmdscale()
plot(mds,col=factor(sexes[rownames(mds),3]))

subi <- colMeans(is.na(gts)) < 0.85 & !grepl("BU",colnames(gts))
weirdos <- which(colMeans(sexsmooth > 0.1) > 0.3 & colMeans(sexsmooth > 0.2) < 0.99 & subi)

plot(sexsmooth[,weirdos[3]],pch=20,ylim=c(0,1),col=factor(vcf[sexsnv,1]),cex=.2)
points(rowMeans(sexsmooth[,9:11]),pch=20,ylim=c(0,1),cex=.2)

colMeans(sexsmooth > 0.1) %>% plot()


####################
# Plots of genetic diversity and divergence across chr5
####################

# window subsets
	# all chromosome 5
	chr5 <- val2[,4] > 7000 & grepl("chr5",lift[,1])
	# non chr5 windows
	subw <- val2[,4] > 7000 & grepl("chr",lift[,1]) & !grepl("chr5",lift[,1])
	# slwin is SDR windows

#############
# generate plot of Fst X vs SDR for 3 pops
# use RColorBrewer set 1?
# brewer.pal(3,"Set1")

pdf(width=8, height=3, file="3pop_chr5_fst.pdf")
# png(width=12, height=4, file="3pop_chr5_fst.png", units="in",res=300)

cols=c("black","#E41A1C")

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(-0.1,1.1),xaxt="n")
q <- quantile(sexfst2[subw,"northX.northY"],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexfst2[chr5,"northX.northY"],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(-0.1,1.1),xaxt="n")
q <- quantile(sexfst2[subw,"admixedX.admixedY"],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexfst2[chr5,"admixedX.admixedY"],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(-0.1,1.1))
q <- quantile(sexfst2[subw,"southX.southY"],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexfst2[chr5,"southX.southY"],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

dev.off()

############


#############
# generate plot of dXY X vs SDR for 3 pops
# use RColorBrewer set 1?
# brewer.pal(3,"Set1")

pdf(width=8, height=3, file="3pop_chr5_dxy.pdf")
# png(width=12, height=4, file="3pop_chr5_fst.png", units="in",res=300)

cols=c("black","#E41A1C")

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(-0.1,1.1),xaxt="n")
q <- quantile(sexdxy2[subw,"northX.northY"],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexdxy2[chr5,"northX.northY"],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(-0.1,1.1),xaxt="n")
q <- quantile(sexdxy2[subw,"admixedX.admixedY"],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexdxy2[chr5,"admixedX.admixedY"],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(0,0.04))
q <- quantile(sexdxy2[subw,"southX.southY"]/val2[subw,4],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points(sexdxy2[chr5,"southX.southY"]/val2[chr5,4],pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)

dev.off()

plot(x=NULL,xlim=c(1,sum(chr5,na.rm=TRUE)),ylim=c(0,0.011))
q <- quantile((sexdxy2[subw,"southX.southY"]-sexdxy2[subw,"southX"])/val2[subw,4],prob=c(0.01,0.99))
rect(xleft=-5000,xright=40000,ytop=q[2],ybottom=q[1],border=NA,col=rgb(0,0,0,0.1))
abline(v=which(!duplicated(val[chr5,1]))[-1],col=rgb(0,0,0,.1),lwd=1.5)
points((sexdxy2[chr5,"southX.southY"]-sexdxy2[chr5,"southX"])/val2[chr5,4]/2,pch=20,cex=.2,col=cols[(slwin[chr5]+1)])
abline(h=0)



############



par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sexdxy2[chr5,"northY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"southY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"admixedY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sexdxy2[chr5,"northY"]/sexdxy2[chr5,"northY.grandX"],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"southY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"admixedY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot((sexdxy2[chr5,"northX"]-sexdxy2[chr5,"northY"])/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"southY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"admixedY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
plot(sexdxy2[chr5,"northY.southY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)
plot(sexdxy2[chr5,"southX.southY"]/val2[chr5,4],pch=20,cex=.2,col=(slwin[chr5]+1))
abline(h=0)

# get pi/dxy values for each chromosome, for each sex
sexdxychr <- c()
for(i in unique(lift[,1])[1:24]){


	chr <- lift[,1] == i
	sexdxychr <- rbind(sexdxychr,colSums(sexdxy[chr,-(1:3)],na.rm=TRUE)/sum(as.numeric(val[chr,4])))


}

sexdxychr <- rbind(sexdxychr,colSums(sexdxy[chr5 & slwin,-(1:3)],na.rm=TRUE)/sum(as.numeric(val[chr5 & slwin,4])))
sexdxychr <- rbind(sexdxychr,colSums(sexdxy[chr5 & !slwin,-(1:3)],na.rm=TRUE)/sum(as.numeric(val[chr5 & !slwin,4])))
rownames(sexdxychr) <- c(unique(lift[,1])[1:24],"SDR","chr5-SDR")

plot(sexdxychr[,"northX"],sexdxychr[,"southX"],xlim=c(0.01,0.02),ylim=c(0.01,0.02),col=(unique(lift[,1])[1:24]=="chr5")+1)
text(sexdxychr[,"northX"],sexdxychr[,"southX"],rownames(sexdxychr))

plot(sexdxychr[,"admixedX"],sexdxychr[,"southX"],xlim=c(0.01,0.02),ylim=c(0.01,0.02),col=(unique(lift[,1])[1:24]=="chr5")+1)
text(sexdxychr[,"admixedX"],sexdxychr[,"southX"],rownames(sexdxychr))

plot(sexdxychr[,"grandX"],sexdxychr[,"southX"],xlim=c(0.004,0.01),ylim=c(0.014,0.02),col=(unique(lift[,1])[1:24]=="chr5")+1)
text(sexdxychr[,"grandX"],sexdxychr[,"southX"],rownames(sexdxychr))

par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
chr <- lift[,1]=="chr5"
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"southX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
points(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"southX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2,col="red")
plot(lift[subw & chr,2]/1e6,sexfst2[subw & chr,"northX.southX"],pch=20,cex=.2)

par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
chr <- lift[,1]=="chr18"
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"southX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX"]/sexdxy2[subw & chr,"northX.grandX"] - sexdxy2[subw & chr,"southX"]/sexdxy2[subw & chr,"northX.grandX"],pch=20,cex=.2)
abline(h=0)

chr <- lift[,1]=="chr18"
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX"]/val2[subw & chr,4],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"southX"]/val2[subw & chr,4],pch=20,cex=.2)
plot(lift[subw & chr,2]/1e6,sexdxy2[subw & chr,"northX.grandX"]/val2[subw & chr,4],pch=20,cex=.2)
abline(h=0)

### Ne for autosomes, x, y, from Wilson Sayres 2018

naut <- function(nm,nf){ (4 * nm * nf)/(nm+nf)}

nchrx <- function(nm,nf){(9 * nm * nf)/((4*nm) + (2*nf))}

nchry <- function(nm){nm/2}

xaut <- function(nm,nf){nchrx(nm,nf)/naut(nm,nf)}

yaut <- function(nm,nf){nchry(nm)/naut(nm,nf)}


# plot of dxy southX.southY, etc
plot(sexdxy2[slwin2,"southX.southY"]/val2[slwin2,4],pch=20,cex=.2,ylim=c(0,0.04))
points(sexdxy2[slwin2,"southX"]/val2[slwin2,4],pch=20,cex=.2,col=rgb(1,0,0,.3))
points((sexdxy2[slwin2,"southX.southY"]-sexdxy2[slwin2,"southX"])/val2[slwin2,4],pch=20,cex=.2,col=rgb(0,1,0,.3))
points((sexdxy2[slwin2,"southX.grandX"])/val2[slwin2,4],pch=20,cex=.2,col=rgb(0,0,1,.3))
# points((sexdxy2[slwin2,"northX"])/val2[slwin2,4],pch=20,cex=.2,col=rgb(1,0,1,.3))




