library(tidyverse)
library(stringr)
library(qvalue)
library(NMF)
library(beeswarm)

smother <- function(x,winsize){
	
	vec <- 1:length(x)
	start <- vec - winsize; start[start < 1] <- 1
	end <- vec + winsize; end[end > length(x)] <- length(x)
	win <- cbind(start,end)
	out <- apply(win, MAR=1,FUN=function(z){sum(x[z[1]:z[2]],na.rm=TRUE)})
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
	
# old scaffold mappings
om <- read.table("~/projects/het_grand/old_scaffold_mappings.txt",stringsAsFactors=FALSE)
rownames(om) <- om[,1]
om2 <- om[-10180,]
rownames(om2) <- om2[,2]


# scaffold to chromosome translations
lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

val <- read.table("~/projects/het_grand_data/popsites.1kb_win.bed.gz",stringsAsFactors=FALSE)
val <- val[lord, ]

d <- duplicated(val[,1])
lgscaf <- cbind(scaf=val[!d,1],chr=lift[!d,1])
rownames(lgscaf) <- lgscaf[,1]

fai <- read.table("GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.fai",stringsAsFactors=FALSE)
rownames(fai) <- fai[,1]

## read in heteroclitus snp-based sex intervals
sexint <- read.table("heteroclitus_SNVsex_intervals.bed",stringsAsFactors=FALSE)

# read in individual counts in 1kb windows. 
sco <- read.table("1kb_sex_counts.txt.gz",stringsAsFactors=FALSE,header=TRUE)
sco <- sco[lord,]

# fix up names, get population IDs
cname <- colnames(sco) %>% gsub(".*align.","",.) %>% gsub(".bam","",.)
cname[7:294] <- gsub("\\.","-",cname[7:294])
cname <- str_replace(cname,regex("(?<=........)_"),".")
colnames(sco) <- cname

# filter individuals
libsize <- colSums(sco[,7:582])
keep <-  libsize > 500000
sco <- sco[,c(1:6,which(keep)+6)]
libsize <- libsize[keep]

# scale by approximate library size:
for(i in 7:dim(sco)[2]){

	sco[,i] <- sco[,i] / sum(sco[,i]) * 1000000
	if(i %% 10 == 0){print(i)}
}

# smooth values
# sco2 <- sco
# for(i in 7:575){ sco2[,i] <- smother(sco2[,i],5);print(i)}


# filter windows??
# if you change this, you must re-run p-value calcs
rs <- rowSums(sco[,7:dim(sco)[2]])
keep <- rs > 100
sco <- sco[keep,]
# sco2 <- sco2[keep,]
lift <- lift[keep,]
val <- val[keep,]

# get sexes, create metadata table
cname <- colnames(sco)
sexes <- read.table("~/projects/het_grand/all_sexes.txt",stringsAsFactors=FALSE)
rownames(sexes) <- sexes[,1]
sexes <- sexes[which(sexes[,1] %in% colnames(sco)),]
sexes <- sexes[colnames(sco)[7:575],]

covmeta <- cbind(sexes,libsize=libsize)
covmeta <- cbind(covmeta,specsex=paste("h",covmeta[,3],sep="_"),stringsAsFactors=FALSE)
covmeta <- cbind(covmeta,popsex=paste(covmeta[,2],covmeta[,3],sep="_"),stringsAsFactors=FALSE)
covmeta[grep("BU",covmeta[,1]),5] <- gsub("h","g",covmeta[grep("BU",covmeta[,1]),5])

gM <- which(covmeta$specsex=="g_M")
gF <- which(covmeta$specsex=="g_F")

hM <- which(covmeta$specsex=="h_M")
hF <- which(covmeta$specsex=="h_F")


# wilcox test for differences in coverage between the sexes
# run to regenerate p-values

# pvals <- cbind(sco[,1:6],gp=NA,hp=NA)

# for(i in 1:dim(pvals)[1]){

# 	pvals[i,7] <- wilcox.test(x=unlist(sco[i,gM]),y=unlist(sco[i,gF]))$p.value
# 	pvals[i,8] <- wilcox.test(x=unlist(sco[i,hM]),y=unlist(sco[i,hF]))$p.value
# 	if(i %% 1000 == 0){print(i)}
# }

# pvals <- cbind(pvals,gq=qvalue(pvals[,"gp"])$qvalues,hq=qvalue(pvals[,"hp"])$qvalues)

# read in p/q values if already calculated:
# already filtered as above. 
pvals <- read.table("sex_coverage_wilcoxonPQvals.txt.gz",stringsAsFactors=FALSE,header=TRUE)

ord <- order(pvals[,2],pvals[,3])
gwin <- mergewin(win=pvals[which(pvals$gq < 0.05),2:4],stat=pvals$gq[which(pvals$gq < 0.05)],qu=0.05,buff=5000,tail="lesser")
hwin <- mergewin(win=pvals[which(pvals$hq < 0.05),2:4],stat=pvals$hq[which(pvals$hq < 0.05)],qu=0.05,buff=5000,tail="lesser")

owin <- rbind(gwin,hwin)
owin <- owin[order(owin[,1],owin[,2]),]
owin <- mergewin(win=owin[,1:3],stat=rep(0,dim(owin)[1]),qu=0.05,buff=5000,tail="lesser")


# add mean coverage for each interval
# also create individual level means for each interval
gint <- c()
for(i in 1:dim(gwin)[1]){

	subw <- sco[,2]==gwin[i,1] & sco[,3]>=gwin[i,2] & sco[,4]<=gwin[i,3] & pvals$gq < 0.05
	tsum <- rowMeans(sco[subw,gM+6]) %>% mean()
	gwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,gF+6]) %>% mean()
	gwin[i,9] <- tsum
	gint <- rbind(gint,c(gwin[i,1:3],colMeans(sco[subw,c(gM+6,gF+6)])))

}
gint <- data.frame(gint)

hint <- c()
for(i in 1:dim(hwin)[1]){

	subw <- sco[,2]==hwin[i,1] & sco[,3]>=hwin[i,2] & sco[,4]<=hwin[i,3] & pvals$hq < 0.05
	tsum <- rowMeans(sco[subw,hM+6]) %>% mean()
	hwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,hF+6]) %>% mean()
	hwin[i,9] <- tsum
	hint <- rbind(hint,c(hwin[i,1:3],colMeans(sco[subw,c(hM+6,hF+6)])))
}
hint <- data.frame(hint)


# read in SNV p/q values (only those < 0.05)
sx <- read.table("sex_chisq.txt.gz",stringsAsFactors=FALSE)

# read in heteroclitus SDR intervals
hetint <- read.table("heteroclitus_SNVsex_intervals.bed",stringsAsFactors=FALSE)

# counts of significant SNVs per scaffold
scount <- group_by(sx,V1) %>% summarize(.,hetcount=sum(V5 < 0.000001,na.rm=TRUE),gcount=sum(V6 < 0.000001,na.rm=TRUE)) %>% data.frame()
scount <- cbind(scount,scaflen=fai[scount[,1],2]/1000000,hden=scount[,2]/fai[scount[,1],2]*1000000)
scount[order(scount[,2]/scount[,4]),]

lgscaf[scount[scount[,2] > 1,1],2] %>% table() %>% sort() %>% data.frame()

fai[scount[scount[,2] > 1,1],2] %>% sum()


plot(rowMeans(sco[pvals$hq < 0.05,hM+6]),rowMeans(sco[pvals$hq < 0.05,hF+6]),pch=20,cex=.6)
abline(0,1)
abline(h=1,v=1)
abline(h=c(0.5,2),v=c(0.5,2),lty=2)

#plot M/F coverage ratio for grandis and heteroclitus for a scaffold
#subscaf <- grepl("NW_012224610.1|NW_012234400.1|NW_012234431.1",sco[,2])
subscaf <- which(grepl("NW_012234378.1",sco[,2]))[1300:1403]
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM+6])/rowMeans(sco[subscaf,hF+6]),2),pch=20,cex=1,col=(pvals[subscaf,"hq"] < 0.05)+1)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,gM+6])/rowMeans(sco[subscaf,gF+6]),2),pch=20,cex=1,col=(pvals[subscaf,"gq"] < 0.05)+1)


# amh coverage by sex
dev.new()
par(mfrow=c(1,2))
subscaf <- grepl("NW_012234285.1",sco[,2])
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ covmeta[c(gM,gF),5])
points(x=jitter(1+(covmeta[c(gM,gF),5]=="g_M")),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))
boxplot(colMeans(sco[subscaf & pvals$hq < 0.05,c(hM,hF)+6],na.rm=TRUE) ~ covmeta[c(hM,hF),5])
points(x=jitter(1+(covmeta[c(hM,hF),5]=="h_M")),y=colMeans(sco[subscaf & pvals$hq < 0.05,c(hM,hF)+6],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))


dev.new()
par(mfrow=c(1,2))
subscaf <- grepl("NW_012234285.1",sco[,2])
beeswarm(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ covmeta[c(gM,gF),5],pch=20,col=c(rgb(0,0,0,.25),rgb(1,0,0,.75)))
bxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ covmeta[c(gM,gF),5],add=TRUE)
beeswarm(colMeans(sco[subscaf & pvals$hq < 0.05,c(hM,hF)+6],na.rm=TRUE) ~ covmeta[c(hM,hF),5],pch=20,col=c(rgb(0,0,0,.25),rgb(1,0,0,.75)))
bxplot(colMeans(sco[subscaf & pvals$hq < 0.05,c(hM,hF)+6],na.rm=TRUE) ~ covmeta[c(hM,hF),5],add=TRUE)

# plot coverage correlation of two regions (or regional means)
# plot individual coverage, sorted by sex
par(mfrow=c(1,2))
	s1 <- grepl("NW_012228907.1",sco[,2])
	s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# s2 <- grepl("NW_012225189.1",sco[,2])
	plot(
		colMeans(sco[s1 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),
		colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),
		col=factor(covmeta[c(gM,gF),5])
		)

	colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) %>% plot(.,col=factor(covmeta[c(gM,gF),5]))


# plot per individual coverage sorted by coverage:
s1 <- grepl("NW_012228907.1",sco[,2])
s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# single copy females
	hfd <- colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)+6]) < 2
	# add a column to metadata file
	covmeta <- cbind(covmeta,cnvsex=covmeta[,5],stringsAsFactors=FALSE)
	covmeta[names(which(hfd)),7] <- gsub("$","_s",covmeta[names(which(hfd)),7])
	covmeta[,7] <- gsub("g_F$", "g_F_m",covmeta[,7])
	covmeta <- cbind(covmeta,spec=gsub("_.","",covmeta[,5]),stringsAsFactors=FALSE)

cord <- order(colMeans(sco[which(s2 & pvals$gq < 0.05),c(gM,gF)+6]))
plot(colMeans(sco[which(s2 & pvals$gq < 0.05),c(gM,gF)+6])[cord],col=factor(covmeta[c(gM,gF),5][cord]))#,pch=as.numeric(as.factor(pop[cord])))

# how correlated are grandis intervals?
# matrix of p-values for linear models
# need to find a way to explain the redundancy of these variables

gintcor <- matrix(nrow=length(gint[,1]),ncol=length(gint[,1]))
for(i in 1:length(gint[,1])){
	for(j in 1:length(gint[,1])){
		tlm <- lm(unlist(gint[i,-c(1:3)]) ~ unlist(gint[j,-c(1:3)]))
		gintcor[i,j] <- summary(tlm)$coefficient[2,4]

	}
}

which(apply(log(gintcor,10),MAR=1,FUN=median) > -2)
plot(apply(log(gintcor,10),MAR=1,FUN=median))

testmat <- as.matrix(gint[,-c(1:3)])
storage.mode(testmat) <- "numeric"
testcor <- cor(t(testmat))
aheatmap(testmat[-c(15),],annCol=covmeta[c(gM,gF),5],distfun="pearson",scale="row")

dat <- as.matrix(data.frame(gint)[,-c(1:3)])
storage.mode(dat) <- "numeric"

subgwin <- which(gwin$count > 2)
gwin.pca <- (prcomp(dat[subgwin,]))
plot(gwin.pca$rotation[,1],gwin.pca$rotation[,2],col=factor(covmeta[c(gM,gF),5]))

	# PC1 explains 98% of variance

# histogram of PC1
range(gwin.pca$rotation[,1])
hist(gwin.pca$rotation[which(covmeta[c(gM,gF),5]=="g_M"),1],breaks=seq(-0.01,0.25,0.005),xlim=c(-0.01,0.25),ylim=c(0,45),border=NA,col=rgb(0,0,0,0.2))
hist(gwin.pca$rotation[which(covmeta[c(gM,gF),5]=="g_F"),1],,breaks=seq(-0.01,0.25,0.005),border=NA,col=rgb(1,0,0,0.2),add=TRUE)

# how correlated are heteroclitus intervals?

testmat <- as.matrix(hint[,-c(1:3)])
storage.mode(testmat) <- "numeric"
testcor <- cor(t(testmat))
aheatmap(testmat[-c(4,15,28,29,33),],annCol=covmeta[c(hM,hF),5],distfun="pearson",scale="row",Colv=NA)


# plot showing the distribution of CNVs among males and females
# this figure goes in the paper. figure 1? 

par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(3,3,1,1))

subgwin <- which(gwin$count > 2)
subgwin <- subgwin[order(gwin[subgwin,8])]

x1 <- 0.75+c(0,cumsum(gwin[subgwin,"count"][-length(subgwin)]))
x2 <- cumsum(gwin[subgwin,"count"])
y1 <- gwin[subgwin,8]
y2 <- gwin[subgwin,9]

plot(NULL,ylim=c(0,30),xlim=c(0,133))
abline(h=1,lty=2,col="gray",lwd=2)
for(i in 1:length(subgwin)){

	points(x=c(x1[i],x2[i]),y=c(y1[i],y1[i]),type="l",lwd=3)
	points(x=c(x1[i],x2[i]),y=c(y2[i],y2[i]),type="l",lwd=3,col="red")
}

subhwin <- which(hwin$count > 2)
subhwin <- subhwin[order(hwin[subhwin,8])]

x1 <- 0.75+c(0,cumsum(hwin[subhwin,"count"][-length(subhwin)]))
x2 <- cumsum(hwin[subhwin,"count"])
y1 <- hwin[subhwin,8]
y2 <- hwin[subhwin,9]

plot(NULL,ylim=c(0,2.5),xlim=c(0,133))
abline(h=1,lty=2,col="gray",lwd=2)
for(i in 1:length(subhwin)){

	points(x=c(x1[i],x2[i]),y=c(y1[i],y1[i]),type="l",lwd=3)
	points(x=c(x1[i],x2[i]),y=c(y2[i],y2[i]),type="l",lwd=3,col="red")
}


################

cypscaf <- which(sco[,2] =="NW_012234324.1")
cypcnv <- cypscaf[1130:1330]
rowMeans(sco[cypcnv,hM+6]) %>% plot(.,ylim=c(0,5))
rowMeans(sco[cypcnv,hF+6]) %>% points(.,ylim=c(0,5),col="blue")

log(rowMeans(sco[cypcnv,hM+6])/rowMeans(sco[cypcnv,hF+6]),2) %>% plot()

for(i in 1:length(cypcnv)){

	tmp <- wilcox.test(unlist(sco[cypcnv[i],hF+6]),unlist(sco[cypcnv[i],hM+6]))
	print(tmp$p.value)
}


#### plot per individual coverage for grandis regions

par(mfrow=c(1,2))
	s1 <- grepl(gwin[subgwin,1][7],sco[,2])
	# s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# s2 <- grepl("NW_012225189.1",sco[,2])
	s2 <- sco[,2] %in% gwin[subgwin,1][1:9]
	plot(
		colMeans(sco[s1 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),
		colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),
		col=factor(covmeta[c(gM,gF),5])
		)

	colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) %>% plot(.,col=factor(covmeta[c(gM,gF),5]))


# male-biased CNV coverage by sex
# try a beeswarm plot for this

dev.new()
par(mfrow=c(1,2))
subscaf <- sco[,2] %in% gwin[subgwin,1][1:9]
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ factor(covmeta[c(gM,gF),7],levels=c("g_F_s","g_F_m","g_M")))
points(x=jitter(as.numeric(factor(covmeta[c(gM,gF),7],levels=c("g_F_s","g_F_m","g_M")))),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))

dev.new()
par(mfrow=c(1,1))
subscaf <- sco[,2] %in% gwin[subgwin,1][1:9]
beeswarm(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ factor(covmeta[c(gM,gF),7],levels=c("g_F_s","g_F_m","g_M")),corral="random",col=c(rgb(1,0,0,.5),rgb(1,0,0,.5),rgb(0,0,0,.5)),pch=20,ylab="relative coverage",xlab="sex class")
bxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)+6],na.rm=TRUE) ~ factor(covmeta[c(gM,gF),7],levels=c("g_F_s","g_F_m","g_M")),add=TRUE)


# could a focal-gene + copy extension with random termination process explain this?
# simulate. 
# number of copy events for 100 lineages. geometric distribution
nc <- rgeom(1000,.025)

# for each copy event, how far do you copy? also geometric. 
out <- matrix(data=1,nrow=1000,ncol=1000)
for(i in 1:1000){

	if(nc[i]==0){next()}
	#extension distance
	ne <- rgeom(nc[i],0.1)
	for(j in 1:nc[i]){

		if(ne[j]==0){next()}
		out[i,1:ne[j]] <- out[i,1:ne[j]] + 1

	}
	}

k <- 1
l <- 15
plot(out[,k],out[,l])
abline(lm(out[,l] ~ out[,k]))

	# it seems like this process could produce the observed pattern. 


# get combined het grand cnv intervals

combwin <- rbind(hwin[hwin$count > 2,],gwin[gwin$count > 2,])
combwin <- combwin[order(combwin[,1],combwin[,2]),]
combwin <- mergewin(win=combwin[,1:3],stat=rep(0,dim(combwin)[2]),qu=1,buff=2000,tail="lesser")[,1:3]
combwin[,3] <- combwin[,3] + 1000
for(i in 1:dim(combwin)[1]){combwin[i,3] <- min(combwin[i,3],fai[combwin[i,1],2])}
combwin[,2] <- combwin[,2] - 1000
combwin[combwin[,2] < 0,2] <- 0
write.table(combwin,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE,file="hetgrand_sex_cnv.txt")

	#$FB -f ~/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz -L <(cat meta/hetgrand.all.list | grep -v "KC-9") -t hetgrand_sex_cnv.txt -k >hetgrand.cnv.vcf


# mapping rate data
load("~/Dropbox/Public/coveragedataV3.RData")
colnames(tab) <- gsub("(........)_","\\1.",colnames(tab))
tab <- tab[,colnames(tab) %in% covmeta[,1]]

unrawprop <- unlist(tab[1,]/tab[9,])

beeswarm(unrawprop ~ covmeta[,5])

wilcox.test(unrawprop ~ covmeta[,3])


# read in sex-DE gene table

sexde <- read.table("hetgrand_sex_biased_expr.txt",stringsAsFactors=FALSE,header=TRUE)

# read in table of variant calls from cnv regions. 

h <- scan(pipe("gzcat hetgrand.sexcnv.vcf.gz | grep CHROM"),what="character")
vcf <- read.table("hetgrand.sexcnv.vcf.gz",stringsAsFactors=FALSE)
colnames(vcf) <- h

# make a vcf metadata table

# sexes switched m for f
rnasex <- read.table("../het_grand/rna_sexes.txt",stringsAsFactors=FALSE)
rownames(rnasex) <- rnasex[,1]
rnasex[rnasex[,2]=="F",2] <- "Ma"
rnasex[rnasex[,2]=="M",2] <- "F"
rnasex[rnasex[,2]=="Ma",2] <- "M"

rnalibsizes <- read.table("hetgrand_rnaseq_normalized_libsizes.txt",stringsAsFactors=FALSE)
rownames(rnalibsizes) <- rnalibsizes[,1]
rnalibsizes <- cbind(rnalibsizes,sex=NA,stringsAsFactors=FALSE)
rnalibsizes[rnasex[-20,1],3] <- rnasex[-20,2]

rnalibsizes[1:72,1] <- gsub("\\.","-",rnalibsizes[1:72,1])
rnalibsizes[,1] <- gsub("KC-F1","KC-F2",rnalibsizes[,1]) %>% gsub("ER-F2","ER-F1",.) %>% gsub("ER-F1-200-([456])","ER-F1-2000-\\1",.)
rownames(rnalibsizes) <- rnalibsizes[,1]

rnameta <- cbind(rnalibsizes,spec="h",stringsAsFactors=FALSE)
rnameta[grep("ARS|ASR",rnameta[,1]),4] <- "g"
rnameta <- cbind(rnameta,stage=NA)
rnameta[grep("ARS|ASR",rnameta[,1]),5] <- str_extract(rnameta[grep("ARS|ASR",rnameta[,1]),1],regex("..(?=.$)"))

# amh scaffold: NW_012234285.1


# calculate allele depth by sex p-values

# a function to extract allele depths
afunc <- function(v,i,scalar=1,meta=NULL){

	subl <- v[i,-c(1:9)] %>% unlist() %>% str_split(.,":") %>% do.call(rbind,.)
	if(dim(subl)[2]>1){
		subl <- str_split(subl[,3],",") %>% do.call(rbind,.)
		class(subl) <- "numeric"
		rownames(subl) <- meta[,1]
		subl[is.na(subl)] <- 0
		subl <- subl/scalar
		}else{return(NA)}

	return(subl)

}


# subvcf <- vcf[vcf[,1]=="NW_012234285.1",]
subvcf <- vcf[!grepl(",",vcf[,5]),]
subvcf <- subvcf[subvcf[,6] > 30,]
subvcfp <- cbind(subvcf[,1:9],subvcf[,covmeta[,1]])
subvcfr <- cbind(subvcf[,1:9],subvcf[,rnameta[,1]])

snpout <- matrix(nrow=dim(subvcf)[1],ncol=5)
snpcov <- matrix(nrow=dim(subvcf)[1],ncol=16)

M <- covmeta[,3]=="M"
F <- covmeta[,3]=="F"
grand <- covmeta[,8]=="g"
ma <- which(M)
fa <- which(F)
hm <- which(M & !grand)
hf <- which(F & !grand)
gm <- which(M & grand)
gf <- which(F & grand)
gfm <- which(covmeta[,7] == "g_F_m")
gfs <- which(covmeta[,7] == "g_F_s")

gr <- which(rnameta[,4] == "g" & grepl("35.$|HA.$|NH.$",rnameta[,1]))
hrm <- which(rnameta[,4] == "h" & rnameta[,3]=="M")
hrf <- which(rnameta[,4] == "h" & rnameta[,3]=="F")

libscale <- covmeta[,4]/median(covmeta[,4])
libscaler <- rnameta[,2]/median(rnameta[,2])

for(i in 1:dim(subvcf)[1]){
# for(i in 9455:9460){

	# dna
	subl <- subvcfp[i,-c(1:9)] %>% unlist() %>% str_split(.,":") %>% do.call(rbind,.)
	if(dim(subl)[2]>1){
		subl <- str_split(subl[,3],",") %>% do.call(rbind,.)
		class(subl) <- "numeric"
		rownames(subl) <- covmeta[,1]
		subl[is.na(subl)] <- 0
		subl <- subl/libscale
	
		m <- subl[ma,] %>% colSums(.,na.rm=TRUE)
		f <- subl[fa,] %>% colSums(.,na.rm=TRUE)
			if(sum(m,f)>0){
				snpout[i,1] <- chisq.test(cbind(m,f))$p.value
				}else{snpout[i,1] <- NA}
	
		m <- subl[hm,] %>% colSums(.,na.rm=TRUE)
		f <- subl[hf,] %>% colSums(.,na.rm=TRUE)
			if(sum(m,f)>0){
				snpout[i,2] <- chisq.test(cbind(m,f))$p.value
				}else{snpout[i,2] <- NA}
			snpcov[i,1:4] <- c(colMeans(subl[hm,]),colMeans(subl[hf,]))
	
		m <- subl[gm,] %>% colSums(.,na.rm=TRUE)
		f <- subl[gfm,] %>% colSums(.,na.rm=TRUE)
			if(sum(m,f)>0){
				snpout[i,3] <- chisq.test(cbind(m,f))$p.value
				}else{snpout[i,3] <- NA}
			snpcov[i,5:8] <- c(colMeans(subl[gm,]),colMeans(subl[gfm,]))
	
		m <- subl[gm,] %>% colSums(.,na.rm=TRUE)
		f <- subl[gfs,] %>% colSums(.,na.rm=TRUE)
			if(sum(m,f)>0){
				snpout[i,4] <- chisq.test(cbind(m,f))$p.value
				}else{snpout[i,4] <- NA}
			snpcov[i,9:10] <- colMeans(subl[gfs,])
		}else{print("skip")}
	# rna
	subl <- subvcfr[i,-c(1:9)] %>% unlist() %>% str_split(.,":") %>% do.call(rbind,.)
	if(dim(subl)[2]>1){
		subl <- str_split(subl[,3],",") %>% do.call(rbind,.)
		class(subl) <- "numeric"
		rownames(subl) <- rnameta[,1]
		subl[is.na(subl)] <- 0
		subl <- subl/libscaler
	
		m <- subl[hrm,] %>% colSums(.,na.rm=TRUE)
		f <- subl[hrf,] %>% colSums(.,na.rm=TRUE)
			if(sum(m,f)>0){
				snpout[i,5] <- chisq.test(cbind(m,f))$p.value
				}else{snpout[i,5] <- NA}
			snpcov[i,11:16] <- c(colMeans(subl[hrm,]),colMeans(subl[hrf,]),colMeans(subl[gr,]))
	}else{print("skip")}

	if((i %% 200)==0){print(i)}
}

colnames(snpout) <- c("mf_p","hf_hm_p","gfm_gm_p","gfs_gm_p","hf_hm_r")

plot(-log(snpout[,4],10),col=factor(subvcf[,1]))


sl <- which(subvcf[,1]=="NW_012234431.1" & rowSums(snpcov[,11:12]) > 0.5)
# sl <- which(snpout[,3] < 1e-3 & subvcf[,1]=="NW_012224575.1")
# sl <- which(snpout[,4] < 1e-3 & snpout[,2] < 1e-3 & subvcf[,1]=="NW_012234285.1")
# sl <- which(snpout[,4] < 1e-3)[5:100]
dev.new()
par(mfrow=c(5,1),mar=c(0,0,0,0),oma=c(3,3,0,0))
yl <- c(0,range(as.matrix(snpcov[sl,1:10]),na.rm=TRUE)[2]*1.3)
barplot(t(snpcov[sl,1:2]),border=NA,ylim=yl)
barplot(t(snpcov[sl,3:4]),border=NA,ylim=yl)
barplot(t(snpcov[sl,5:6]),border=NA,ylim=yl)
barplot(t(snpcov[sl,7:8]),border=NA,ylim=yl)
barplot(t(snpcov[sl,9:10]),border=NA,ylim=yl)

dev.new()
par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(3,3,0,0))
yl <- c(0,range(as.matrix(snpcov[sl,11:16]),na.rm=TRUE)[2]*1.3)
barplot(t(snpcov[sl,11:12]),border=NA,ylim=yl)
barplot(t(snpcov[sl,13:14]),border=NA,ylim=yl)
barplot(t(snpcov[sl,15:16]),border=NA,ylim=yl)

subvcf[sl,1:5]
#  cnv in grandis, but not heteroclitus. looks sex-DE in heteroclitus.
	# NW_012225185.1	GeneID:105921442

# snx1 is next to a Y deletion on NW_012224575.1. expression is half in males what it is in females. cool!
# kars has a Y duplication just downstream though expression is lower in males. 
# hsp70 on NW_012234431.1 is a shared cnv, but not DE and doesn't have much or any sequence divergence



order(snpout[,3]) %>% head()
vcf[order(snpout[,3]) %>% head(.,n=15),1:5]
vcf[which(snpout[,2] < 1e-4 & snpout[,3] < 1e-4),1:5]

# get one snv
	i <- 10811

	subl <- afunc(subvcfp,i,scalar=libscale,meta=covmeta)

	m1 <- subl[hm,] %>% colSums(.,na.rm=TRUE)
	f1 <- subl[hf,] %>% colSums(.,na.rm=TRUE)
	
	m2 <- subl[gm,] %>% colSums(.,na.rm=TRUE)
	f2 <- subl[gfm,] %>% colSums(.,na.rm=TRUE)
	f3 <- subl[gfs,] %>% colSums(.,na.rm=TRUE)

vcf[i,1:5]
rbind(m1,f1,m2,f2,f3)

	i <- 10811

	subl <- afunc(subvcfr,i,scalar=libscaler,meta=rnameta)
	
	m <- subl[hrm,] %>% colSums(.,na.rm=TRUE)
	f <- subl[hrf,] %>% colSums(.,na.rm=TRUE)
	g <- subl[gr,] %>% colSums(.,na.rm=TRUE)

rbind(m,f,g)



dev.new()
par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(3,3,0,0))
cl <- 1
subl[hrm,][order(rnameta[hrm,5],subl[hrm,1]),] %>% t() %>% barplot()
subl[hrf,][order(rnameta[hrf,5],subl[hrf,1]),] %>% t() %>% barplot()
subl[which(rnameta[,4]=="g"),][order(rnameta[which(rnameta[,4]=="g"),5],subl[which(rnameta[,4]=="g"),cl]),] %>% t() %>% barplot()
subl[gr,][order(rnameta[gr,5],subl[gr,cl]),] %>% t() %>% barplot()

subl[gm,][order(subl[gm,cl]),] %>% t() %>% barplot()


sublg <- 0
for(i in sl[c(25:28,33:35)]){sublg <- sublg + afunc(subvcfr,i,scalar=libscaler,meta=rnameta)}
sublg <- sublg/7
sublg <- afunc(subvcfr,sl[27],scalar=libscaler,meta=rnameta)

suble <- afunc(subvcfr,sl[13],scalar=libscaler,meta=rnameta)

sublep <- afunc(subvcfp,sl[13],scalar=libscale,meta=covmeta)

subl[grep("ARS",rownames(subl)),]
subl[grep("35.$|HA.$|NH.$",rownames(subl)),]
subl[grep("ARS",rownames(subl)),] %>% (function(x){x[,1]/rowSums(x,na.rm=TRUE)}) %>% sort() %>% plot()
subl[grep("35.$|HA.$|NH.$",rownames(subl)),] %>% (function(x){x[,1]/rowSums(x,na.rm=TRUE)}) %>% sort() %>% plot()
subl[grep("35.$|HA.$|NH.$",rownames(subl)),] %>% (function(x){o <- order(rowSums(x)); t(x[o,]) %>% barplot()})

subl[grep("-F",rownames(subl)),]
subl[grep("-F",rownames(subl)),] %>% (function(x){x[,1]/rowSums(x,na.rm=TRUE)}) %>% sort() %>% plot()
subl[grep("-F",rownames(subl)),] %>% (function(x){o <- order(rowSums(x)); t(x[o,]) %>% barplot()})


library(ape)
library(seqinr)
library(phytools)
amhal <- read.alignment("het_grand_amh_aligned_trimmed.fasta",format="fasta")
dist.alignment(amhal,gap=0) %>% nj() %>% midpoint.root() %>% plot(.,show.tip.label=FALSE)

amhalm <- str_split(amhal$seq,"") %>% do.call(rbind,.)
amhalm[amhalm=="-"] <- NA
amhalm <- amhalm[c(1,6,2,3,4,5,7),]
amhd <- matrix(nrow=7,ncol=7)
rownames(amhd) <- amhal$nam[c(1,6,2,3,4,5,7)]
colnames(amhd) <- amhal$nam[c(1,6,2,3,4,5,7)]
amhd <- amhd

# first 380 have a gap for one sequence, wonky nj result with negative branch lengths. still the same topology though. 
subs <- 380:1368
for(i in 1:7){
	for(j in 1:7){

		amhd[i,j] <- 1-mean(amhalm[i,subs]==amhalm[j,subs],na.rm=TRUE)

	}

}
as.dist(amhd) %>% nj() %>% midpoint.root() %>% plot(.,show.tip.label=FALSE)



apply(amhalm,MAR=2,FUN=function(x){table(x) %>% length()}) %>% sort() %>% table()

