library(tidyverse)
library(stringr)
library(qvalue)

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

# add mean coverage for each interval
for(i in 1:dim(gwin)[1]){

	subw <- sco[,2]==gwin[i,1] & sco[,3]>=gwin[i,2] & sco[,4]<=gwin[i,3]
	tsum <- rowMeans(sco[subw,gM]) %>% mean()
	gwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,gF]) %>% mean()
	gwin[i,9] <- tsum

}

for(i in 1:dim(hwin)[1]){

	subw <- sco[,2]==hwin[i,1] & sco[,3]>=hwin[i,2] & sco[,4]<=hwin[i,3]
	tsum <- rowMeans(sco[subw,hM]) %>% mean()
	hwin[i,8] <- tsum
	tsum <- rowMeans(sco[subw,hF]) %>% mean()
	hwin[i,9] <- tsum

}


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

subscaf <- sco[,2] == "NW_012234431.1"
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM])/rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(pvals[subscaf,"hq"] < 0.05)+1)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,gM])/rowMeans(sco[subscaf,gF]),2),pch=20,cex=1,col=(pvals[subscaf,"gq"] < 0.05)+1)

boxplot(colSums(sco[subscaf & pvals$gq < 0.05,-c(1:6)]) ~ specsex)

plot(
	colMeans(sco[sco[,2]=="NW_012234431.1" & pvals$gq < 0.05,c(gM,gF)]),
	colMeans(sco[sco[,2]=="NW_012234400.1" & pvals$gq < 0.05,c(gM,gF)])
	)
