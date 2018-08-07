library(stringr)
library(dplyr)
library(ape)
library(phytools)
library(magrittr)
library(viridis)

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




makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

options(scipen=99)

segs <- c(
	"northM","northF",
	"southM","southF",
	"admixedM","admixedF",
	"grandM","grandF") 

om <- read.table("~/projects/het_grand/old_scaffold_mappings.txt",stringsAsFactors=FALSE)
rownames(om) <- om[,1]
om2 <- om[-10180,]
rownames(om2) <- om2[,2]

lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

val <- read.table("~/projects/het_grand_data/popsites.1kb_win.bed.gz",stringsAsFactors=FALSE)
val <- val[lord, ]
val2 <- val
val2[,4] <- smother(val2[,4],10)
subv <- val[,4] > 200

fst <- read.table("~/projects/het_grand_data/fst_dxy_sex.gz",stringsAsFactors=FALSE)
pnames <- c("scaf","start","end",segs,combn(segs,2) %>% apply(.,MAR=2,FUN=paste,collapse="."))
colnames(fst) <- pnames

for(i in 4:39){ fst[,i] <- as.numeric(fst[,i])}
for(i in 4:39){ fst[which(fst[,i] > 200),i] <- NA}


fst <- fst[lord,]
fst2 <- fst
for(i in 4:39){ fst2[,i] <- smother(fst2[,i],10);print(i)}

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
colnames(vcf) <- h

# get only biallelic variants
vcf <- vcf[!grepl(",",vcf[,5]),]

gts <- as.matrix(vcf[,-(1:9)])
class(gts) <- "numeric"

# get sex table
sexes <- read.table("../het_grand/all_sexes.txt",stringsAsFactors=FALSE)
rownames(sexes) <- sexes[,1]

# polarize by female heteroclitus
hetf <- which(sexes[colnames(vcf)[-(1:9)],3]=="F" & !grepl("BU",colnames(vcf)[-(1:9)]))
hetm <- which(sexes[colnames(vcf)[-(1:9)],3]=="M" & !grepl("BU",colnames(vcf)[-(1:9)]))
hetff <- rowMeans(gts[,hetf],na.rm=TRUE)
hetfm <- rowMeans(gts[,hetm],na.rm=TRUE)

sexsnv <- hetff < 0.1 & hetfm > 0.4
subi <- colMeans(is.na(gts)) < 0.8

plot(hetff,hetfm,pch=20,cex=.2,col=rgb(0,0,0,.2))

colMeans(gts[sexsnv,],na.rm=TRUE) %>% plot(.,pch=20,col=factor(sexes[colnames(gts),3]))

sexscore <- colMeans(gts[sexsnv,],na.rm=TRUE)

sexsmooth <- gts[sexsnv,]
for(i in 1:dim(gts)[2]){sexsmooth[,i] <- smotherM(x=sexsmooth[,i],winsize=100); print(i)}


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

image(xymat[c(1,2,5,6,3,4,7,8),c(1,2,5,6,3,4,7,8)])

image(xymat[c(1,2,5,6,3,4),c(1,2,5,6,3,4)])

as.dist(xymat) %>% nj() %>% midpoint.root() %>% plot()

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



subw <- val2[,4] > 7000 & grepl("chr",lift[,1])

chr5 <- val2[,4] > 7000 & grepl("chr5",lift[,1])

plot(sexfst2[chr5,"northX.northY"],pch=20,cex=.2,col=slwin[chr5]+1)

plot((sexdxy2[chr5,"southX.southY"]-sexdxy2[chr5,"southX"])/val2[chr5,4],pch=20,cex=.2,col=slwin[chr5]+1)
points((sexdxy2[chr5,"southX.grandX"]-sexdxy2[chr5,"southX"])/val2[chr5,4],pch=20,cex=.2,col="blue")


# check some weird individuals
# weird individuals: F-20  F-2 F-47 SH-16


plot(sexsmooth[,1],pch=20,cex=.2,col="blue",ylim=c(0,1))
points(sexsmooth[,2],pch=20,cex=.2,col="black")
points(sexsmooth[,3],pch=20,cex=.2,col="green")


plot(sexsmooth,pch=20,cex=.2,col=factor(vcf[sexsnv,1]),ylim=c(0,1))


# "female" scaffolds in weird individuals: "NW_012224869.1","NW_012234324.1"
colMeans(gts[sexsnv & grepl("NW_012224869.1|NW_012234324.1",vcf[,1]),],na.rm=TRUE) %>% plot(.,col=factor(sexes[colnames(gts),2]))

mds <- t(gts[sexsnv & vcf[,1]=="NW_012224869.1",subi & !grepl("BU", colnames(gts))]) %>% dist() %>% cmdscale()
plot(mds,col=factor(sexes[rownames(mds),3]))

subi <- colMeans(is.na(gts)) < 0.8
weirdos <- which(colMeans(sexsmooth > 0.2) > 0.3 & colMeans(sexsmooth > 0.2) < 0.99 & subi)

plot(sexsmooth[,weirdos[1]],pch=20,ylim=c(0,1),col=factor(vcf[sexsnv,1]))

colMeans(sexsmooth > 0.2) %>% plot()


