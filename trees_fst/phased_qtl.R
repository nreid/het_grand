library(dplyr)
library(stringr)
library(qvalue)
library(viridis)

# some functions for reading in and processing data

# extracts genotype probability for each genotype in a VCF column
	# useful for indentifying individuals or loci with systematically poor genotype calls after phasing
	# didn't really make a difference to analysis
mprob <- function(fmt){

	fmt <- gsub(".*:","",fmt) %>%
		str_split(.,",") %>%
		do.call(rbind,.)
	class(fmt) <- "numeric"
	fmt <- apply(fmt,MAR=1,FUN=max) %>% unlist()	
	return(fmt)

	}

# next four functions are subfunctions of read.phased.somm
splif <- function(x) {
	x <- str_split(x,"\\|") %>% do.call(rbind,.)
	return(x)
	}

makehap <- function(x){
	
	hap <- c()
	for(i in 1:dim(x)[2]){

		temp <- splif(x[,i])
		colnames(temp) <- paste(colnames(x)[i],c("a","b"),sep=".")
		hap <- cbind(hap,temp)
		print(i)
		}
	class(hap) <- "numeric"
	return(hap)
	}

ngt <- function(x){

	x <- gsub("0\\|0","0",x)
	x <- gsub("0\\|1","1",x)
	x <- gsub("1\\|0","1",x)
	x <- gsub("1\\|1","2",x)
	class(x) <- "numeric"
	return(x)
	}

read.phased.somm <- function(file){
	# read in data
	cline <- paste("gzcat ", file, " | head -n 15000 | grep ^#CHROM",sep="",collapse="")
	h <- scan(pipe(cline),what="character",sep="\t")
	vcf <- read.table(file,stringsAsFactors=FALSE)
	colnames(vcf) <- h
	
	#create genotype matrix
	gt <- as.matrix(vcf[,10:dim(vcf)[2]])
	gt <- gsub(":.*","",gt)
	#haplotype matrix
	hap <- makehap(gt)
	#make genotypes numeric
	gt <- ngt(gt)

	# get AR2,DR2
	met <- str_split(vcf[,8],";") %>% do.call(rbind,.) %>% gsub(".*=","",.)
	colnames(met) <- c("ar2","dr2","af")
	class(met) <- "numeric"

	# 
	keepl <- rowMeans(gt/2) > 0.15 & rowMeans(gt/2) < 0.85 & rowMeans(gt==1) < 0.65 & rowMeans(gt==1) > 0.2 & met[,1] > 0.85
	vcf <- vcf[keepl,]
	gt <- gt[keepl,]
	hap <- hap[keepl,]
	met <- met[keepl,]

	return(list(vcf=vcf,gt=gt,hap=hap,met=met))
	}

# calculate chi-square stats
chis <- function(som,ki=NULL){
	if(is.null(ki)){ki <- TRUE}
	som$gt <- som$gt[,ki]
	chib <- c()
	dphen <- pmat[colnames(som$gt),"dphen"]
	for(i in 1:dim(som$gt)[1]){
	
		chib[i] <- tryCatch(chisq.test(som$gt[i,],dphen)$p.value, error=function(err){ NA })
		if((i %% 200)==0){print(i)}
	
		}
	return(chib)
	}



# some metadata
	lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
	win <- read.table("~/Downloads/1kb_win.bed",stringsAsFactors=FALSE)
	win200 <- read.table("~/Downloads/200kb50kb_win.bed",stringsAsFactors=FALSE)
	win500 <- read.table("~/Downloads/500kb50kb_win.bed",stringsAsFactors=FALSE)
	om <- read.table("~/projects/het_grand/old_scaffold_mappings.txt",stringsAsFactors=FALSE)
	rownames(om) <- om[,1]
	om2 <- om[-10180,]
	rownames(om2) <- om2[,2]

# read in phased vcfs
	somm1 <- read.phased.somm("SOMM075.beagle.vcf.gz")
	somm2 <- read.phased.somm("SOMM086.beagle.vcf.gz")
	somm3 <- read.phased.somm("SOMM087.beagle.vcf.gz")
	somm4 <- read.phased.somm("SOMM088.beagle.vcf.gz")

# more metadata
	somm <- read.table("~/Downloads/Library.SOMM.table.txt",stringsAsFactors=FALSE,header=TRUE)
	pmat <- data.frame(sam=somm$library_sample,phen=somm$Phenotype,dphen=as.numeric(somm$Phenotype > 2),pop=somm$pop,stringsAsFactors=FALSE)
	rownames(pmat) <- somm$library_sample

# calculate chi square for g-p association
	chi1 <- chis(somm1)
	chi2 <- chis(somm2)
	chi3 <- chis(somm3)
	chi4 <- chis(somm4)

# plot chi-square
	par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col=factor(somm1$vcf[,1]))
	plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col=factor(somm2$vcf[,1]))
	plot(-log(qvalue(chi3)$qvalues,10),pch=20,cex=.2,col=factor(somm3$vcf[,1]))
	plot(-log(qvalue(chi4)$qvalues,10),pch=20,cex=.2,col=factor(somm4$vcf[,1]))
	
	par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
	plot(-log(chi1,10),pch=20,cex=.2,col=factor(somm1$vcf[,1]))
	plot(-log(chi2,10),pch=20,cex=.2,col=factor(somm2$vcf[,1]))
	plot(-log(chi3,10),pch=20,cex=.2,col=factor(somm3$vcf[,1]))
	plot(-log(chi4,10),pch=20,cex=.2,col=factor(somm4$vcf[,1]))


# windowed haplotype association test. 

	winchi <- function(win,som,pma=pmat){

	win <- unlist(win)
	ind <- som$vcf[,1]==win[1] & som$vcf[,2] > as.numeric(win[2]) & som$vcf[,2] <= as.numeric(win[3])
	if(sum(ind)==0){return(NA)}
	sgt <- som$gt[ind,]
	dphen <- pma[colnames(som$gt),"dphen"]
		
	if(sum(ind)==1){
		pv <- chisq.test(sgt,dphen)$p.value
		return(pv)
	}

	mds <- t(sgt) %>% dist() %>% cmdscale()

	cl <- kmeans(cbind(mds[,1]),3,iter.max=100,nstart=50)
	# print(sum(ind))
	pval <- chisq.test(cl$cluster,dphen)$p.value
	return(pval)
	
	}
	
	
	
	winchi_all <- function(w,s){

	wchi <- c()
	for(i in 1:length(w[,1])){

		wchi[i] <- winchi(w[i,],s)
		if((i %% 500)==0){print(i)}
	}
	return(wchi)
	}
	
	
	wc1 <- winchi_all(win500,somm1)
	wc2 <- winchi_all(win500,somm2)
	wc3 <- winchi_all(win500,somm3)
	wc4 <- winchi_all(win500,somm4)
	
	xl <- c(1,17000)
	par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
	plot(-log(qvalue(wc1)$qvalues,10),pch=20,cex=.2,col=factor(win500[,1]),xlim=xl)
	plot(-log(qvalue(wc2)$qvalues,10),pch=20,cex=.2,col=factor(win500[,1]),xlim=xl)
	plot(-log(qvalue(wc3)$qvalues,10),pch=20,cex=.2,col=factor(win500[,1]),xlim=xl)
	plot(-log(qvalue(wc4)$qvalues,10),pch=20,cex=.2,col=factor(win500[,1]),xlim=xl)
	
	xl <- c(1,20000)
	par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,1,1))
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col=factor(somm1$vcf[,1]),xlim=xl)
	plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col=factor(somm2$vcf[,1]),xlim=xl)
	plot(-log(qvalue(chi3)$qvalues,10),pch=20,cex=.2,col=factor(somm3$vcf[,1]),xlim=xl)
	plot(-log(qvalue(chi4)$qvalues,10),pch=20,cex=.2,col=factor(somm4$vcf[,1]),xlim=xl)





############### SOMM088 - ER cross:
	# AHRb qtl core region: 103853:103898
	# chr13 qtl core region: 74300:74600
		# after removing 
	# chr7 qtl core region: 42400:42600


	# genotype heatmap
	ord <- (t(somm4$gt[103853:103898,]) %>% dist() %>% hclust())$order
	image(somm4$gt[103853:103898,ord])
	
	#mds plot
	mds <- t(somm4$gt[103853:103898,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])
	
	# take only individuals not homozygous for AHRb haplotype, recalculate chi-square
	keepi <- rownames(mds)[mds[,1] > -2 & mds[,1] < 2]
	chi4e <- chis(somm4,keepi)
	
	par(mfrow=c(2,1))
	plot(-log(chi4e,10),pch=20,cex=.4,col=factor(somm4$vcf[,1]),xlim=c(42400,42800))
	plot(-log(qvalue(chi4e)$qvalues,10),pch=20,cex=.4,col=factor(somm4$vcf[,1]))
	
	plot(-log(chi4,10),pch=20,cex=.2,col="black",xlim=c(49200,49700))
	
	ord <- (t(somm4$gt[49200:49700,keepi]) %>% dist() %>% hclust())$order
	image(somm4$gt[49200:49700,keepi][,ord])
	
	ord <- (t(somm4$gt[49200:49700,]) %>% dist() %>% hclust())$order
	image(somm4$gt[49200:49700,][,ord])
	
	
	mds <- t(somm4$gt[49200:49700,keepi]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])
	
	mds2 <- t(somm4$gt[49200:49700,]) %>% dist() %>% cmdscale()
	plot(mds2,col="white")
	text(jitter(mds2[,1],amount=1),jitter(mds2[,2],amount=1),pmat[rownames(mds),"phen"])
	
	
	mds <- t(somm4$gt[103853:103898,]) %>% dist() %>% cmdscale()
	mds2 <- t(somm4$gt[74300:74600,]) %>% dist() %>% cmdscale()
	mds3 <- t(somm4$gt[42400:42600,]) %>% dist() %>% cmdscale()
	mds4 <- t(somm4$gt[49200:49700,]) %>% dist() %>% cmdscale()
	
	par(mfrow=c(1,3))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	plot(mds[,1],mds3[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds3[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	plot(mds[,1],mds4[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds4[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	
	
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),round(mds3[,1]),col=pmat[rownames(mds),"phen"]+1)
	
	
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(mds[,1],mds2[,1],pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	
	cl <- kmeans(cbind(2*mds[,1],mds2[,1]),9,iter.max=50,nstart=50)
	par(mfrow=c(1,3))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),cl$cluster,col=pmat[rownames(mds),"phen"]+1)
	
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	
	plot(jitter(mds3[,1],amount=0.25),jitter(cl$cluster),col="white")
	text(jitter(mds3[,1],amount=0.25),jitter(cl$cluster),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"phen"]+1)
	
	ers <- read.table("~/Dropbox/Public/outliertables/outliersER.txt",stringsAsFactors=FALSE,header=TRUE)
	ert <- read.table("~/Dropbox/Public/outliertables/outliersER_NW.txt",stringsAsFactors=FALSE,header=TRUE)
	erl <- read.table("~/Dropbox/Public/outliertables/outliersER_lift.bed",stringsAsFactors=FALSE)
	eror <- rep(0,length(somm4$gt[,1]))
	for(i in 1:length(erl[,1])){
	
		temp <- somm4$vcf[,1]==erl[i,1] & erl[i,2] < somm4$vcf[,2] & erl[i,3] >= somm4$vcf[,2]
		eror[temp] <- ert[i,"rank"]
		print(i)
		}
	
	par(mfrow=c(2,1))
	plot(-log(qvalue(chi4)$qvalues,10),pch=20,cex=.2,col=factor(somm4$vcf[,1]))
	plot(-log(qvalue(chi4)$qvalues,10),pch=20,cex=.2,col=eror)





############### SOMM075 - NBH

	# chr2: small chunk of AIP region: 7900:8110
	# chr18: AHRb region: 75050:75200 significantly associated after removing AIP homozygotes
	# chr8: 33100:33275 significantly associated after removing AIP homozygotes
	# chr12: 52200:52300 not significant (p=0.06) after controlling or AIP, but suspicious looking
	# chr24: 98800:99000 not significant (p=0.06) after controlling for AIP
	
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col=factor(somm1$vcf[,1]),xlim=c(52200,52300))
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col="black",xlim=c(33000,33500),ylim=c(0,4))
	
	
	# genotype heatmap
	ord <- (t(somm1$gt[7900:8110,]) %>% dist() %>% hclust())$order
	image(somm1$gt[7900:8110,ord])
	
	#mds plot
	mds <- t(somm1$gt[8090:8110,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])
	
	# take only individuals not homozygous for AIP haplotype, recalculate chi-square
	keepi <- rownames(mds)[ mds[,1] < 1]
	chi1e <- chis(somm1,keepi)
	
	par(mfrow=c(3,1))
	plot(-log(chi1,10),pch=20,cex=.4,col=factor(somm1$vcf[,1]),xlim=c(75050,75200))
	plot(-log(chi1e,10),pch=20,cex=.4,col=factor(somm1$vcf[,1]),xlim=c(75050,75200))
	plot(-log(qvalue(chi4e)$qvalues,10),pch=20,cex=.4,col=factor(somm1$vcf[,1]),xlim=c(75050,75200))
	# chr18 remains significantly associated after removing AIP homozygous individuals
		# only 8 out of 56 remaining individuals have healthy phenotype scores in this group
	# chr
	
	
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col=factor(somm1$vcf[,1]),xlim=c(98800,99000))
	plot(-log(qvalue(chi1)$qvalues,10),pch=20,cex=.2,col="black",xlim=c(52200,52300),ylim=c(0,4))
	
	ord <- (t(somm1$gt[98800:99000,keepi]) %>% dist() %>% hclust())$order
	image(somm1$gt[98800:99000,keepi][,ord])
	
	ord <- (t(somm1$gt[98800:99000,]) %>% dist() %>% hclust())$order
	image(somm1$gt[98800:99000,][,ord])
	
	mds <- t(somm1$gt[98800:99000,keepi]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds[,2])))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	mds2 <- t(somm1$gt[98800:99000,]) %>% dist() %>% cmdscale()
	plot(mds2,col="white",,xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds2[,2])))
	text(jitter(mds2[,1],amount=1),jitter(mds2[,2],amount=1),pmat[rownames(mds2),"phen"],col=pmat[rownames(mds2),"dphen"]+1)
	
	mds <- t(somm1$gt[7900:8110,]) %>% dist() %>% cmdscale()
	mds2 <- t(somm1$gt[75050:75200,]) %>% dist() %>% cmdscale()
	mds3 <- t(somm1$gt[33100:33275,]) %>% dist() %>% cmdscale()
	mds4 <- t(somm1$gt[52200:52300,]) %>% dist() %>% cmdscale()
	mds5 <- t(somm1$gt[98800:99000 ,]) %>% dist() %>% cmdscale()
	
	par(mfrow=c(2,2))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds3[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds3[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds4[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds4[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds5[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds5[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds5[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	
	
	cl <- kmeans(mds3[,1],3,nstart=50,iter=100)





############### SOMM086 - NYC 
	
	# chr2: 8975:9016 small chunk of AIP region: 
	# chr18: 76850:76950 AHRb region: significantly associated after removing AIP homozygotes
	# chr11: 50250:50450. heterozygotes have the advantage. 
		# three alleles in parents: AB x AC cross. BC appears to have a lethal interaction with AIPt/AIPt
	   #  cl   0  1  4  5 # column=phenotype, row=cluster
		  # 1  4  0  4 12 # BC genotype
 		 #  2 11  3  4  4 # AB genotype
  		#   3 15  3  1  6 # AC genotype
  		#   4 10  2  6 11 # AA genotype


	# chr8: 38775:38880 not much to say.. 

	par(mfrow=c(2,1))
	plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col=factor(somm2$vcf[,1]))
	plot(-log(qvalue(wc2)$qvalues,10),pch=20,cex=.5,col=factor(win500[,1]),xlim=c(1,27000))
	#plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col="black",ylim=c(0,4))
	
	# genotype heatmap
	ord <- (t(somm2$gt[8975:9016,]) %>% dist() %>% hclust())$order
	image(somm2$gt[8975:9016,ord])
	
	#mds plot
	mds <- t(somm2$gt[8975:9016,]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=1.25*range(mds[,2]))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	# take only individuals not homozygous for AIP haplotype, recalculate chi-square
	keepi <- rownames(mds)[ mds[,1] > -4]
	chi2e <- chis(somm2,keepi)

	par(mfrow=c(3,1))
	plot(-log(chi2,10),pch=20,cex=.4,col=factor(somm2$vcf[,1]))
	plot(-log(chi2e,10),pch=20,cex=.4,col=factor(somm2$vcf[,1]))
	plot(-log(qvalue(chi2e)$qvalues,10),pch=20,cex=.4,col=factor(somm2$vcf[,1]))

	#AHRb retains a very strong association after controlling for AIP.
		# not surprising, it was stronger to begin with

	par(mfrow=c(2,1))
	plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col="black",xlim=c(38775,38880))
	plot(-log(qvalue(wc2)$qvalues,10),pch=20,cex=.2,col="black",xlim=c(5800,5950))
	
	ord <- (t(somm2$gt[38775:38880,keepi]) %>% dist() %>% hclust())$order
	image(somm2$gt[38775:38880,keepi][,ord])
	
	ord <- (t(somm2$gt[38775:38880,]) %>% dist() %>% hclust())$order
	image(somm2$gt[38775:38880,][,ord])
	
	mds <- t(somm2$gt[38775:38880,keepi]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds[,2])))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	mds2 <- t(somm2$gt[38775:38880,]) %>% dist() %>% cmdscale()
	plot(mds2,col="white",,xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds2[,2])))
	text(jitter(mds2[,1],amount=1),jitter(mds2[,2],amount=1),pmat[rownames(mds2),"phen"],col=pmat[rownames(mds2),"dphen"]+1)

	mds <- t(somm2$gt[8975:9016,]) %>% dist() %>% cmdscale()
	mds2 <- t(somm2$gt[76850:76950,]) %>% dist() %>% cmdscale()
	mds3 <- t(somm2$gt[50250:50450,]) %>% dist() %>% cmdscale()
	mds4 <- t(somm2$gt[38775:38880,]) %>% dist() %>% cmdscale()
	# mds5 <- t(somm1$gt[98800:99000 ,]) %>% dist() %>% cmdscale()
	
	par(mfrow=c(2,2))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds3[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds3[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds2[,1],mds[,1],col="white",xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds2[,1],amount=1),jitter(mds3[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds4[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds4[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds2[,1],mds4[,1],col="white",xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds2[,1],amount=1),jitter(mds4[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)



	nycs <- read.table("~/Dropbox/Public/outliertables/outliersNYC.txt",stringsAsFactors=FALSE,header=TRUE)
	nyct <- read.table("~/Dropbox/Public/outliertables/outliersNYC_NW.txt",stringsAsFactors=FALSE,header=TRUE)
	nycl <- read.table("~/Dropbox/Public/outliertables/outliersNYC_lift.bed",stringsAsFactors=FALSE)



############### SOMM087 - BP # half sibs? outcross?
	
	# chr2: 9675:9775 AIP region
	# chr18: 88150:88250 AHRb region
	# chr24: 116450:116600 4 haps in cross. only one associated with tolerance. contributes after removing AHRb
	# chr13: 

	par(mfrow=c(2,1))
	plot(-log(qvalue(chi3)$qvalues,10),pch=20,cex=.2,col=factor(somm3$vcf[,1]))
	plot(-log(qvalue(wc3)$qvalues,10),pch=20,cex=.2,col=factor(win500[,1]))
	#plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col="black",ylim=c(0,4))
	
	# genotype heatmap
	ord <- (t(somm3$gt[9675:9775,]) %>% dist() %>% hclust())$order
	image(somm3$gt[9675:9775,ord])
	
	#mds plot
	mds <- t(somm3$gt[88150:88250,]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=1.25*range(mds[,2]))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	# take only individuals not homozygous for AIP haplotype, recalculate chi-square
	keepi <- rownames(mds)[ mds[,1] > -4]
	chi3e <- chis(somm3,keepi)

	par(mfrow=c(3,1))
	plot(-log(chi3,10),pch=20,cex=.4,col=factor(somm3$vcf[,1]))
	plot(-log(chi3e,10),pch=20,cex=.4,col=factor(somm3$vcf[,1]))
	plot(-log(qvalue(chi3e)$qvalues,10),pch=20,cex=.4,col=factor(somm3$vcf[,1]))

	#AHRb retains a very strong association after controlling for AIP.
		# not surprising
	# chr14 may show an association when controlling for AIP

	par(mfrow=c(2,1))
	plot(-log((chi3),10),pch=20,cex=.2,col=factor(somm3$vcf[,1]),xlim=c(64900,65250)) #116450,116600
	plot(-log((wc3),10),pch=20,cex=.2,col=factor(win500[,1]),xlim=c(15000,18000))
	
	ord <- (t(somm3$gt[64900:65250,keepi]) %>% dist() %>% hclust())$order
	image(somm3$gt[64900:65250,keepi][,ord])
	
	ord <- (t(somm3$gt[64900:65250,]) %>% dist() %>% hclust())$order
	image(somm3$gt[64900:65250,][,ord])
	
	mds <- t(somm3$gt[64900:65250,keepi]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds[,2])))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	mds2 <- t(somm3$gt[64900:65250,]) %>% dist() %>% cmdscale()
	plot(mds2,col="white",,xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds2[,2])))
	text(jitter(mds2[,1],amount=1),jitter(mds2[,2],amount=1),pmat[rownames(mds2),"phen"],col=pmat[rownames(mds2),"dphen"]+1)

	mds <- t(somm3$gt[9675:9775,]) %>% dist() %>% cmdscale()
	mds2 <- t(somm3$gt[88150:88250,]) %>% dist() %>% cmdscale()
	mds3 <- t(somm3$gt[116450:116600,]) %>% dist() %>% cmdscale()
	mds4 <- t(somm3$gt[64900:65250,]) %>% dist() %>% cmdscale()
	
	par(mfrow=c(2,2))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds3[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds3[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds2[,1],mds3[,2],col="white",xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds3[,1])))
	text(jitter(mds2[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds[,1],mds4[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds4[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	plot(mds2[,1],mds4[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	text(jitter(mds2[,1],amount=1),jitter(mds4[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)




############### SOMM088 - ER

	# chr18: AHRb qtl core region: 
		# maybe: 103853:103898
		# 
	# chr13 qtl core region: 74300:74600
		# after removing 
	# chr7 qtl core region: 42400:42600
	

	par(mfrow=c(2,1))
	plot(-log(qvalue(chi4)$qvalues,10),pch=20,cex=.2,col=factor(somm4$vcf[,1]))
	plot(-log((wc4),10),pch=20,cex=.5,col=factor(win500[,1]),xlim=c(1,17000))
	#plot(-log(qvalue(chi2)$qvalues,10),pch=20,cex=.2,col="black",ylim=c(0,4))
	
	# genotype heatmap
	ord <- (t(somm4$gt[102190:102210,]) %>% dist() %>% hclust())$order
	image(somm4$gt[102190:102210,ord])
	
	#mds plot
	mds <- t(somm4$gt[102199:102201,]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=3*range(mds[,2]))
	text(jitter(mds[,1],amount=.1),jitter(mds[,2],amount=0),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	# take only individuals not homozygous for AIP haplotype, recalculate chi-square
	#keepi <- rownames(mds)[ mds[,1] > -4]
	keepi <- names(which(somm4$gt[103648,]!=2))
	chi4e <- chis(somm4,keepi)

	par(mfrow=c(3,1))
	plot(-log(chi4,10),pch=20,cex=.4,col=factor(somm4$vcf[,1]))
	plot(-log(chi4e,10),pch=20,cex=.4,col=factor(somm4$vcf[,1]))
	plot(-log(qvalue(chi4e)$qvalues,10),pch=20,cex=.4,col=factor(somm4$vcf[,1]))

	#AHRb retains a very strong association after controlling for AIP.
		# not surprising
	# chr14 may show an association when controlling for AIP

	par(mfrow=c(2,1))
	plot(-log((chi4),10),pch=20,cex=1,col=factor(somm4$vcf[,1]),xlim=c(102000,104000)) #116450,116600
	plot(-log((wc4),10),pch=20,cex=.5,col=factor(win500[,1]),xlim=c(8000,15000))
	
	ord <- (t(somm4$gt[74300:74600,keepi]) %>% dist() %>% hclust())$order
	image(somm4$gt[74300:74600,keepi][,ord])
	
	ord <- (t(somm4$gt[74300:74600,]) %>% dist() %>% hclust())$order
	image(somm4$gt[74300:74600,][,ord])
	
	mds <- t(somm4$gt[74300:74600,keepi]) %>% dist() %>% cmdscale()
	plot(mds,col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds[,2])))
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	
	mds2 <- t(somm4$gt[74300:74600,]) %>% dist() %>% cmdscale()
	plot(mds2,col="white",,xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds2[,2])))
	text(jitter(mds2[,1],amount=1),jitter(mds2[,2],amount=1),pmat[rownames(mds2),"phen"],col=pmat[rownames(mds2),"dphen"]+1)

	mds <- t(somm4$gt[103853:103898,]) %>% dist() %>% cmdscale()
	mds2 <- t(somm4$gt[74300:74600,]) %>% dist() %>% cmdscale()
	# mds3 <- t(somm4$gt[116450:116600,]) %>% dist() %>% cmdscale()
	# mds4 <- t(somm4$gt[64900:65250,]) %>% dist() %>% cmdscale()
	# mds5 <- t(somm1$gt[98800:99000 ,]) %>% dist() %>% cmdscale()
	
	par(mfrow=c(2,2))
	plot(mds[,1],mds2[,1],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds2[,1])))
	text(jitter(mds[,1],amount=1),jitter(mds2[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds[,1],mds3[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	# text(jitter(mds[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds[,1],mds3[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds3[,1])))
	# text(jitter(mds[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds2[,1],mds3[,2],col="white",xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds3[,1])))
	# text(jitter(mds2[,1],amount=1),jitter(mds3[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds[,1],mds4[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	# text(jitter(mds[,1],amount=1),jitter(mds4[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds2[,1],mds4[,2],col="white",xlim=1.25*range(mds[,1]),ylim=c(1.25*range(mds4[,1])))
	# text(jitter(mds2[,1],amount=1),jitter(mds4[,2],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)
	# plot(mds2[,1],mds4[,1],col="white",xlim=1.25*range(mds2[,1]),ylim=c(1.25*range(mds4[,1])))
	# text(jitter(mds2[,1],amount=1),jitter(mds4[,1],amount=1),pmat[rownames(mds),"phen"],col=pmat[rownames(mds),"dphen"]+1)

	plot(somm4$gt[103859,],somm4$gt[76045,],col="white",xlim=c(-1,3),ylim=c(-1,3))
	text(jitter(somm4$gt[103859,]),jitter(somm4$gt[76045,]),pmat[colnames(somm4$gt),"phen"],col=pmat[colnames(somm4$gt),"dphen"]+1)


nw11 <- c()
for(i in grep("NW_",somm2$vcf[,1])){

	nw11 <- c(nw11, chisq.test(somm2$gt[50433,],somm2$gt[i,])$p.value)
	if((i %% 200)==0){print(i)}

}


get_scafs <- function(chr,start,end){

	scafs <- win[which(lift[,1]==chr & lift[,2] > end & lift[,3] < start),]
	data.frame(NW_scaf=unique(scafs[,1]),old_scaf=om2[unique(scafs[,1]),1],counts=as.vector(table(scafs[,1])))

	}



fai <- read.table("~/projects/het_grand_data/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai",stringsAsFactors=FALSE)

revc <- function(x){

	x <- str_split(x,"")
	x <- sapply(x,FUN=function(y){
		trans <- c(A="T",C="G",G="C",T="A")
		y <- trans[y] %>% rev() %>% paste(.,sep="",collapse="")
		y 
		})
	x

}

# check against pop gen data. 

# AIP NW_012234293.1.vcf.gz
	# orientation is opposite to map assembly
	# 

	scaflen <- fai[fai[,1]=="NW_012234293.1",2]
	offset <- lift[win[,1]=="NW_012234293.1",2] %>% range()
	
	hnames <- scan("~/projects/het_grand/vcf_header.txt",what="character",sep="\t")
	vcf <- read.table("~/projects/het_grand_data/NW_012234293.1.vcf.gz",stringsAsFactors=FALSE)
	colnames(vcf) <- hnames
	vcf <- vcf[!grepl(",",vcf[,5]),]
	vcf[,1] <- "chr2"
	# renumber variants
	vcf[,2] <- vcf[,2] + nchar(vcf[,4]) - 1
	vcf[,2] <- -1 * (vcf[,2] - scaflen - 1) + offset[1]
	vcf[,4] <- revc(vcf[,4])
	vcf[,5] <- revc(vcf[,5])
	vcf <- vcf[order(vcf[,2]),]
	
	pgt <- as.matrix(vcf[,10:297])
		class(pgt) <- "numeric"


	## SOMM1 NBH
	ind1 <- somm1$vcf[,1]=="chr2" & somm1$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm1$vcf[somm1$vcf[,1]=="chr2",2]

	gm <- somm1$hap[ind1,]
	gmg <- somm1$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))

	ntol <- grep("BP|NYC",pop)
	nsen <- grep("F|SH",pop)

	fdif <- (rowMeans(gp[,ntol],na.rm=TRUE) - rowMeans(gp[,nsen],na.rm=TRUE))	
	stat <- colMeans(gmg[fdif > 0.25,]/2)
	# randomize to check
		# stat <- colMeans(gmg[sample(1:60,14),]/2)

	plot(stat)
	text(1:96,stat,pmat[rownames(mds),"phen"])

	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	stat <- colMeans(gm[fdif > 0.25,])
	# randomize to check
		stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm1$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])

	# SOMM2 NYC: AIP haplotype maybe homozygous in cross. 
	ind1 <- somm2$vcf[,1]=="chr2" & somm2$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm2$vcf[somm2$vcf[,1]=="chr2",2]

	gm <- somm2$hap[ind1,]
	gmg <- somm2$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))

	ntol <- grep("BP|NYC",pop)
	nsen <- grep("F|SH",pop)

	fdif <- (rowMeans(gp[,ntol],na.rm=TRUE) - rowMeans(gp[,nsen],na.rm=TRUE))	
	stat <- colMeans(gmg[fdif > 0.25,]/2)
	# randomize to check
		# stat <- colMeans(gmg[sample(1:81,24),]/2)

	plot(stat)
	text(1:96,stat,pmat[rownames(mds),"phen"])

	mean(colMeans(gmg[,stat>0.8]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(.,pch=20,col=ind2+1)

	mds <- t(somm2$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])


	# SOMM3 BP: AIP haplotype probably homozygous in cross. 
	ind1 <- somm3$vcf[,1]=="chr2" & somm3$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm3$vcf[somm3$vcf[,1]=="chr2",2]

	gm <- somm3$hap[ind1,]
	gmg <- somm3$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))

	ntol <- grep("BP|NYC",pop)
	nsen <- grep("F|SH",pop)

	fdif <- (rowMeans(gp[,ntol],na.rm=TRUE) - rowMeans(gp[,nsen],na.rm=TRUE))	
	stat <- colMeans(gmg[fdif > 0.25,]/2)
	# randomize to check
		# stat <- colMeans(gmg[sample(1:101,24),]/2)

	plot(stat)
	text(1:96,stat,pmat[rownames(mds),"phen"])

	mean(colMeans(gmg[,stat>0.8]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	stat <- colMeans(gm[fdif > 0.25,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(.,pch=20,col=ind2+1)

	mds <- t(somm3$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])


	# SOMM4 ER: AIP haplotype not likely homozygous in cross. 
	# high heterozygosity in ER/ER chromosme cluster
	ind1 <- somm4$vcf[,1]=="chr2" & somm4$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm4$vcf[somm4$vcf[,1]=="chr2",2]

	gm <- somm4$hap[ind1,]
	gmg <- somm4$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))

	ntol <- grep("ER",pop)
	nsen <- grep("KC",pop)

	fdif <- (rowMeans(gp[,ntol],na.rm=TRUE) - rowMeans(gp[,nsen],na.rm=TRUE))	
	stat <- colMeans(gmg[fdif > 0.5,]/2)
	# randomize to check
		stat <- colMeans(gmg[sample(1:107,15),]/2)

	plot(stat,ylim=c(0,1))
	text(1:96,stat,pmat[rownames(mds),"phen"])

	mean(colMeans(gmg[,stat<0.2]==1))

	#		plot(colMeans(gmg)/2,colMeans(gmg==1))
	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])


	stat <- colMeans(gm[fdif > 0.5,])
	# randomize to check
		stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])
	

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(.,pch=20,col=ind2+1)

	mds <- t(somm4$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])


# AHR NW_012234474.1.vcf.gz
	# orientation same as map assembly
	# 

	scaflen <- fai[fai[,1]=="NW_012234474.1",2]
	offset <- lift[win[,1]=="NW_012234474.1",2] %>% range()
	
	hnames <- scan("~/projects/het_grand/vcf_header.txt",what="character",sep="\t")
	vcf <- read.table("~/projects/het_grand_data/NW_012234474.1.vcf.gz",stringsAsFactors=FALSE)
	colnames(vcf) <- hnames
	vcf <- vcf[!grepl(",",vcf[,5]),]
	vcf[,1] <- "chr1"
	# renumber variants
	vcf[,2] <- vcf[,2]+offset[1]
	
	pgt <- as.matrix(vcf[,10:297])
		class(pgt) <- "numeric"



	# SOMM4 ER:
	# Not enough information to clear this up. 

	del <- read.table("~/projects/het_grand/heteroclitus_del_gt.txt",stringsAsFactors=FALSE)
	erdel <- del[del[,4]=="ER" & del[,2]==0,1]

	ind1 <- somm4$vcf[,1]=="chr2" & somm4$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm4$vcf[somm4$vcf[,1]=="chr2",2]

	gm <- somm4$hap[ind1,]
	gmg <- somm4$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))

	ntol <- which(colnames(gp) %in% erdel)
	nsen <- grep("KC",pop)

	fdif <- (rowMeans(gp[,ntol],na.rm=TRUE) - rowMeans(gp[,nsen],na.rm=TRUE))	
	stat <- colMeans(gmg[which(fdif > 0.25),]/2)
	# randomize to check
		stat <- colMeans(gmg[sample(1:107,15),]/2)

	plot(stat,ylim=c(0,1))
	text(1:96,stat,pmat[rownames(mds),"phen"])

	mean(colMeans(gmg[,stat<0.2]==1))

	#	plot(colMeans(gmg)/2,colMeans(gmg==1))
	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])


	stat <- colMeans(gm[fdif > 0.5,])
	# randomize to check
		stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])
	

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm4$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"])


	# 





# cyp1a NW_012234324.1:500000-1500000.vcf.gz
	# orientation is opposite to map assembly
	# 

	scaflen <- fai[fai[,1]=="NW_012234324.1",2]
	offset <- lift[win[,1]=="NW_012234324.1",2] %>% range()
	
	hnames <- scan("~/projects/het_grand/vcf_header.txt",what="character",sep="\t")
	vcf <- read.table("~/projects/het_grand_data/NW_012234324.1:500000-1500000.vcf.gz",stringsAsFactors=FALSE)
	colnames(vcf) <- hnames
	vcf <- vcf[!grepl(",",vcf[,5]),]
	vcf[,1] <- "chr5"
	# renumber variants
	vcf[,2] <- vcf[,2]+offset[1]
	
	pgt <- as.matrix(vcf[,10:297])
		class(pgt) <- "numeric"

	
	sexes <- read.table("~/projects/het_grand/all_sexes.txt",stringsAsFactors=FALSE)
	rownames(sexes) <- sexes[,1]

	## SOMM1 NBH: haplotype appears not to be in the cross. 
	ind1 <- somm1$vcf[,1]=="chr5" & somm1$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm1$vcf[somm1$vcf[,1]=="chr5",2]

	gm <- somm1$hap[ind1,]
	gmg <- somm1$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.4,]/2)	

	plot(stat,ylim=c(0,1),col=sex2+1,pch=20)
	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:83,7),]/2)
	points(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.4,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(sstat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm1$gt[ind1,][fdif > 0.25,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)

	# SOMM2 NYC: one copy of CYP1A selected haplotpye in the cross. copy number unknown. 
	ind1 <- somm2$vcf[,1]=="chr5" & somm2$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm2$vcf[somm2$vcf[,1]=="chr5",2]

	gm <- somm2$hap[ind1,]
	gmg <- somm2$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.25,]/2)	

	plot(stat,ylim=c(0,1),col=sex2+1,pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.4,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(sstat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm2$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)



	# SOMM3 BP: cyp1a at best one copy in parents. 
	ind1 <- somm3$vcf[,1]=="chr5" & somm3$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm3$vcf[somm3$vcf[,1]=="chr5",2]

	gm <- somm3$hap[ind1,]
	gmg <- somm3$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.4,]/2)	

	plot(stat,ylim=c(0,1),col=sex2+1,pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.25,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm3$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)




	# SOMM4 ER: two copies in cross. 

	ind1 <- somm4$vcf[,1]=="chr5" & somm4$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm4$vcf[somm4$vcf[,1]=="chr5",2]

	gm <- somm4$hap[ind1,]
	gmg <- somm4$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("ER",pop)
	nsen <- grepl("KC",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.4,]/2)	

	plot(stat,ylim=c(0,1),col=sex2+1,pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.25,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm4$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)



# aip NW_012224618.1:1-500000.vcf.gz
	# orientation is same as map assembly
	# 

	scaflen <- fai[fai[,1]=="NW_012224618.1",2]
	offset <- lift[win[,1]=="NW_012224618.1",2] %>% range()
	
	hnames <- scan("~/projects/het_grand/vcf_header.txt",what="character",sep="\t")
	vcf <- read.table("~/projects/het_grand_data/NW_012224618.1:1-500000.vcf.gz",stringsAsFactors=FALSE)
	colnames(vcf) <- hnames
	vcf <- vcf[!grepl(",",vcf[,5]),]
	vcf[,1] <- "chr18"
	# renumber variants
	vcf[,2] <- vcf[,2]+offset[1]
	
	pgt <- as.matrix(vcf[,10:297])
		class(pgt) <- "numeric"

	
	sexes <- read.table("~/projects/het_grand/all_sexes.txt",stringsAsFactors=FALSE)
	rownames(sexes) <- sexes[,1]

	## SOMM1 NBH: haplotype appears not to be in the cross. 
	ind1 <- somm1$vcf[,1]=="chr18" & somm1$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm1$vcf[somm1$vcf[,1]=="chr18",2]

	gm <- somm1$hap[ind1,]
	gmg <- somm1$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.2,]/2)	

	plot(stat,ylim=c(0,1),pch=20)
	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:83,7),]/2)
	points(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.4,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(sstat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm1$gt[ind1,][fdif > 0.25,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)

	# SOMM2 NYC: one copy of CYP1A selected haplotpye in the cross. copy number unknown. 
	ind1 <- somm2$vcf[,1]=="chr18" & somm2$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm2$vcf[somm2$vcf[,1]=="chr18",2]

	gm <- somm2$hap[ind1,]
	gmg <- somm2$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.25,]/2)	

	plot(stat,ylim=c(0,1),col=sex2+1,pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.4,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(sstat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm2$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)



	# SOMM3 BP: cyp1a at best one copy in parents. 
	ind1 <- somm3$vcf[,1]=="chr18" & somm3$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm3$vcf[somm3$vcf[,1]=="chr18",2]

	gm <- somm3$hap[ind1,]
	gmg <- somm3$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("BP|NYC",pop)
	nsen <- grepl("F|SH",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	sex2 <- sstat > 0.3
	stat <- colMeans(gmg[fdif > 0.2,]/2)	

	plot(stat,ylim=c(0,1),pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.25,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm3$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)




	# SOMM4 ER: two copies in cross. 
	
	ind1 <- somm4$vcf[,1]=="chr18" & somm4$vcf[,2] %in% (vcf[,2])
	ind2 <- vcf[,2] %in% somm4$vcf[somm4$vcf[,1]=="chr18",2]

	gm <- somm4$hap[ind1,]
	gmg <- somm4$gt[ind1,]
	gp <- as.matrix(vcf[ind2,10:297])
		class(gp) <- "numeric"
	pop <- gsub("-.*","",colnames(gp))
	sex <- sexes[colnames(gp),3]

	ntol <- grepl("ER",pop)
	nsen <- grepl("KC",pop)
	male <- sex=="M"
	female <- sex=="F"

	fdif <- (rowMeans(gp[,ntol&female],na.rm=TRUE) - rowMeans(gp[,nsen&female],na.rm=TRUE))	
	sdif <- (rowMeans(gp[,male],na.rm=TRUE) - rowMeans(gp[,female],na.rm=TRUE))	

	sstat <- colMeans(gmg[sdif > 0.4,]/2)
	stat <- colMeans(gmg[fdif > 0.25,]/2)	

	plot(stat,ylim=c(0,1),pch=20)
	rowMeans(gp[,ntol],na.rm=TRUE)[fdif > 0.25]
	rowMeans(gp[,nsen],na.rm=TRUE)[fdif > 0.25]

	# text(1:96,stat,pmat[rownames(mds),"phen"])

	# randomize to check
	stat <- colMeans(gmg[sample(1:64,8),]/2)
	plot(stat,ylim=c(0,1),pch=20,col=rgb(sex2+0,0,0,.5))



	mean(colMeans(gmg[,stat>0.3]==1))

	plot(stat,colMeans(gmg==1))
	image(gmg[,order(stat)])
	image(gmg[fdif>0.25,order(stat)])

	# FOR HAPLOTYPES
	stat <- colMeans(gm[fdif > 0.25,])
	sstat <- colMeans(gm[sdif > 0.4,])
	# randomize to check
		# stat <- colMeans(gmg[sample(1:107,27),]/2)

	plot(stat,ylim=c(0,1))

	plot(stat,colMeans(gm==1))
	image(gm[,order(stat)])
	image(gm[fdif>0.25,order(stat)])

	(rowMeans(pgt[,ntol],na.rm=TRUE) - rowMeans(pgt[,nsen],na.rm=TRUE)) %>% plot(x=vcf[,2],y=.,pch=20,col=ind2+1)

	mds <- t(somm4$gt[ind1,]) %>% dist() %>% cmdscale()
	plot(mds,col="white")
	text(jitter(mds[,1],amount=1),jitter(mds[,2],amount=1),pmat[rownames(mds),"phen"],col=sex2+1)



