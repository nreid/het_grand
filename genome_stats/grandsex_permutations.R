library(tidyverse)
library(stringr)
library(qvalue)
library(NMF)

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

# scale by approximate library size:
for(i in 7:dim(sco)[2]){

	sco[,i] <- sco[,i] / sum(sco[,i]) * 1000000
	if(i %% 10 == 0){print(i)}
}

# smooth values
# sco2 <- sco
# for(i in 7:575){ sco2[,i] <- smother(sco2[,i],5);print(i)}


# filter windows??
rs <- rowSums(sco[,7:dim(sco)[2]])
keep <- rs > 100
sco <- sco[keep,]
sco2 <- sco2[keep,]
lift <- lift[keep,]


# get sexes
cname <- colnames(sco)
sexes <- read.table("~/projects/het_grand/all_sexes.txt",stringsAsFactors=FALSE)
rownames(sexes) <- sexes[,1]

gM <- which(sexes[cname,3] == "M" & grepl("^BU",sexes[cname,1]))
gF <- which(sexes[cname,3] == "F" & grepl("^BU",sexes[cname,1]))

hM <- which(sexes[cname,3] == "M" & !grepl("^BU",sexes[cname,1]))
hF <- which(sexes[cname,3] == "F" & !grepl("^BU",sexes[cname,1]))

specsex <- rep("h",dim(sco)[2]-6)
specsex[which(grepl("BU",cname))-6] <- "g"
specsex <- paste(specsex,sexes[cname[7:length(cname)],3],sep="_")
popsex <- paste(sexes[cname[-c(1:6)],2],sexes[cname[-c(1:6)],3],sep="_")
pop <- paste(sexes[cname[-c(1:6)],2],sep="_")


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
	tsum <- rowMeans(sco[subw,gM]) %>% mean()
	gwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,gF]) %>% mean()
	gwin[i,9] <- tsum
	gint <- rbind(gint,c(gwin[i,1:3],colMeans(sco[subw,c(gM,gF)])))

}
gint <- data.frame(gint)

hint <- c()
for(i in 1:dim(hwin)[1]){

	subw <- sco[,2]==hwin[i,1] & sco[,3]>=hwin[i,2] & sco[,4]<=hwin[i,3] & pvals$hq < 0.05
	tsum <- rowMeans(sco[subw,hM]) %>% mean()
	hwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,hF]) %>% mean()
	hwin[i,9] <- tsum
	hint <- rbind(hint,c(hwin[i,1:3],colMeans(sco[subw,c(hM,hF)])))
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


plot(rowMeans(sco[pvals$hq < 0.05,hM]),rowMeans(sco[pvals$hq < 0.05,hF]),pch=20,cex=.6)
abline(0,1)
abline(h=1,v=1)
abline(h=c(0.5,2),v=c(0.5,2),lty=2)

#plot M/F coverage ratio for grandis and heteroclitus for a scaffold
#subscaf <- grepl("NW_012224610.1|NW_012234400.1|NW_012234431.1",sco[,2])
subscaf <- grepl("NW_012225189.1",sco[,2])
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM])/rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(pvals[subscaf,"hq"] < 0.05)+1)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,gM])/rowMeans(sco[subscaf,gF]),2),pch=20,cex=1,col=(pvals[subscaf,"gq"] < 0.05)+1)


# amh coverage by sex
dev.new()
par(mfrow=c(1,2))
subscaf <- grepl("NW_012234285.1",sco[,2])
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) ~ specsex[c(gM,gF)-6])
points(x=jitter(1+(specsex[c(gM,gF)-6]=="g_M")),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(hM,hF)],na.rm=TRUE) ~ specsex[c(hM,hF)-6])
points(x=jitter(1+(specsex[c(hM,hF)-6]=="h_M")),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(hM,hF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))

# female higher coverage grandis region
dev.new()
par(mfrow=c(1,2))
subscaf <- grepl("NW_012234348.1",sco[,2])
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) ~ specsex[c(gM,gF)-6])
points(x=jitter(1+(specsex[c(gM,gF)-6]=="g_M")),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(hM,hF)],na.rm=TRUE) ~ specsex[c(hM,hF)-6])
points(x=jitter(1+(specsex[c(hM,hF)-6]=="h_M")),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(hM,hF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))


# plot coverage correlation of two regions (or regional means)
# plot individual coverage, sorted by sex
par(mfrow=c(1,2))
	s1 <- grepl("NW_012228907.1",sco[,2])
	s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# s2 <- grepl("NW_012225189.1",sco[,2])
	plot(
		colMeans(sco[s1 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),
		colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),
		col=factor(specsex[c(gM,gF)-6])
		)

	colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) %>% plot(.,col=factor(specsex[c(gM,gF)-6]))


# plot per individual coverage sorted by coverage:
s1 <- grepl("NW_012228907.1",sco[,2])
s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# single copy females
	hfd <- colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)]) < 2

cord <- order(colMeans(sco[which(s2 & pvals$gq < 0.05),c(gM,gF)]))
plot(colMeans(sco[which(s2 & pvals$gq < 0.05),c(gM,gF)])[cord],col=factor(specsex[c(gM,gF)-6][cord]))#,pch=as.numeric(as.factor(pop[cord])))

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
aheatmap(testmat[-c(15),],annCol=specsex[c(gM,gF)-6],distfun="pearson",scale="row")

dat <- as.matrix(data.frame(gint)[,-c(1:3)])
storage.mode(dat) <- "numeric"

gwin.pca <- (prcomp(dat[subgwin,]))
plot(gwin.pca$rotation[,1],gwin.pca$rotation[,2],col=factor(specsex[c(gM,gF)-6]))

	# PC1 explains 98% of variance

# histogram of PC1
range(gwin.pca$rotation[,1])
hist(gwin.pca$rotation[which(specsex[c(gM,gF)-6]=="g_M"),1],breaks=seq(-0.01,0.25,0.005),xlim=c(-0.01,0.25),ylim=c(0,45),border=NA,col=rgb(0,0,0,0.2))
hist(gwin.pca$rotation[which(specsex[c(gM,gF)-6]=="g_F"),1],,breaks=seq(-0.01,0.25,0.005),border=NA,col=rgb(1,0,0,0.2),add=TRUE)

# how correlated are heteroclitus intervals?

testmat <- as.matrix(hint[,-c(1:3)])
storage.mode(testmat) <- "numeric"
testcor <- cor(t(testmat))
aheatmap(testmat[-c(4,15,28,29,33),],annCol=specsex[c(hM,hF)-6],distfun="pearson",scale="row",Colv=NA)


# examine some regions per base per individual
bams <- scan("popgenbams.list",what="character")
bams <- gsub(".*\\/","",bams) %>% gsub(".bam","",.) %>% str_replace(.,"_",".")
dep <- read.table("NW_012234285.1:140000-180000.depth.gz",stringsAsFactors=FALSE)
colnames(dep) <- c("scaf","pos",bams)

rowMeans(dep[,colnames(sco[,gM])]) %>% plot(x=dep[,2],y=.,pch=20,cex=.2)
rowMeans(dep[,colnames(sco[,gF])]) %>% points(x=dep[,2],y=.,col="blue",pch=20,cex=.2)
rowMeans(dep[,names(which(hfd))]) %>% points(x=dep[,2],y=.,col="red",pch=20,cex=.2)
rowMeans(dep[,grep("-",colnames(dep))]) %>% points(x=dep[,2],y=.,col="green",pch=20,cex=.2)


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




cypscaf <- which(sco[,2] =="NW_012234324.1")
cypcnv <- cypscaf[1130:1330]
rowMeans(sco[cypcnv,hM]) %>% plot(.,ylim=c(0,5))
rowMeans(sco[cypcnv,hF]) %>% points(.,ylim=c(0,5),col="blue")

log(rowMeans(sco[cypcnv,hM])/rowMeans(sco[cypcnv,hF]),2) %>% plot()

for(i in 1:length(cypcnv)){

	tmp <- wilcox.test(unlist(sco[cypcnv[i],hF]),unlist(sco[cypcnv[i],hM]))
	print(tmp$p.value)
}


#### plot per individual coverage for grandis regions

par(mfrow=c(1,2))
	s1 <- grepl(gwin[subgwin,1][7],sco[,2])
	# s2 <- grepl("NW_012228907.1|NW_012224610.1|NW_012234400.1|NW_012234431.1|NW_012225850.1|NW_012234285.1",sco[,2])
	# s2 <- grepl("NW_012225189.1",sco[,2])
	s2 <- sco[,2] %in% gwin[subgwin,1][1:9]
	plot(
		colMeans(sco[s1 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),
		colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),
		col=factor(specsex[c(gM,gF)-6])
		)

	colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) %>% plot(.,col=factor(specsex[c(gM,gF)-6]))

# a vector giving single-copy grandis females, cnv females and males

tmp <- colMeans(sco[s2 & pvals$gq < 0.05,gF],na.rm=TRUE)
cnvsex <- specsex
cnvsex[gF-6][tmp < 2] <- "g_F_s"
cnvsex[cnvsex=="g_F"] <- "g_F_m"

# male-biased CNV coverage by sex
dev.new()
par(mfrow=c(1,2))
subscaf <- sco[,2] %in% gwin[subgwin,1][1:9]
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) ~ factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))
points(x=jitter(as.numeric(factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))

subscaf <- sco[,2] %in% gwin[subgwin,1][3]
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) ~ factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))
points(x=jitter(as.numeric(factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))

subscaf <- sco[,2] %in% gwin[subgwin,1][9]
boxplot(colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE) ~ factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))
points(x=jitter(as.numeric(factor(cnvsex[c(gM,gF)-6],levels=c("g_F_s","g_F_m","g_M")))),y=colMeans(sco[subscaf & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE),pch=20,col=rgb(0,0,0,.5))


s1 <- grepl(gwin[subgwin,1][3],sco[,2])
s2 <- grepl(gwin[subgwin,1][9],sco[,2])

x <- colMeans(sco[s1 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE)
y <- colMeans(sco[s2 & pvals$gq < 0.05,c(gM,gF)],na.rm=TRUE)

# fit <- lm(I(y - 0) ~ 0 + x)
fit <- lm(y ~ x)

plot(x, y,xlim=c(0,7),ylim=c(0,90))
abline((fit))
summary(fit)

# could a focal-gene + copy extension with random termination process explain this?

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

	#$FB -f ~/fhet_genome/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fasta.bgz -L <(cat meta/hetgrand.all.list | grep -v "KC-9") -t hetgrand_sex_cnv.txt -k 


# mapping rate data
load("~/Dropbox/Public/coveragedataV3.RData")


sexes <- read.table("~/projects/het_grand/all_sexes.txt",stringsAsFactors=FALSE)
sexes[,1] <- gsub("\\.","_",sexes[,1])
rownames(sexes) <- sexes[,1]

unrawprop <- unlist(tab[1,sexes[colnames(tab),1]]/tab[9,sexes[colnames(tab),1]])

boxplot(unrawprop ~ (grepl("BU",names(unrawprop))+1) + sexes[names(unrawprop),3])

wilcox.test(unrawprop[!grepl("BU",names(unrawprop))] ~ sexes[names(unrawprop),3][!grepl("BU",names(unrawprop))])

