library(tidyverse)
library(stringr)
library(mscr)

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

data(sexchrs)
sexscaf <- om[sexchrs,2]

# scaffold to chromosome translations
lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

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
sco2 <- sco
for(i in 7:575){ sco2[,i] <- smother(sco2[,i],5);print(i)}


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

# permutations for coverage among sexes
# function to resample columns 
resam <- function(n1,n2,vec){

	sam <- sample(vec)
	return(list(s1=sam[1:n1],s2=sam[(n1+1):length(sam)]))

}

# permutations for grandis:
# iteratively refine low p-values to reduce computational burden

# empirical difference in coverage between sexes. 
	# should matter that I use sum instead of mean b/c permutations
gdiff <- rowSums(sco[,gM]) - rowSums(sco[,gF])


# upper and lower tail p-value initialization
gpvalu <- rep(0,length(gdiff))
gpvall <- rep(0,length(gdiff))

# loop through 100 randomizations
for(i in 1:100){

	sam <- resam(125,161,c(gM,gF))
	rgdiff <- rowSums(sco[,sam[[1]]]) - rowSums(sco[,sam[[2]]])
	gpvalu <- gpvalu + (rgdiff <= gdiff)
	gpvall <- gpvall + (rgdiff >= gdiff)
	if(i %% 10 == 0){print(i)}

}

# turn counts into p-values
gpvalu <- gpvalu/100
gpvall <- gpvall/100


# refine p-values for subset with p < 0.1
	# do 1k loops for these windows
subw <- gpvalu < 0.1 | gpvall < 0.1
scow <- sco[subw,]

gdiffw <- gdiff[subw]
gpvaluw <- rep(0,sum(subw))
gpvallw <- rep(0,sum(subw))

for(i in 1:1000){

	sam <- resam(125,161,c(gM,gF))
	rgdiff <- rowSums(scow[,sam[[1]]]) - rowSums(scow[,sam[[2]]])
	gpvaluw <- gpvaluw + (rgdiff <= gdiffw)
	gpvallw <- gpvallw + (rgdiff >= gdiffw)
	if(i %% 20 == 0){print(i)}

}

# turn counts into p-values
gpvaluw <- gpvaluw/1000
gpvallw <- gpvallw/1000

# place refined p-values into original p-value vectors

gpvalu[subw] <- gpvaluw
gpvall[subw] <- gpvallw

# iterate, this time refine p-values < 0.02
	# do 10k permutations

subw <- gpvalu < 0.02 | gpvall < 0.02
scow <- sco[subw,]

gdiffw <- gdiff[subw]
gpvaluw <- rep(0,sum(subw))
gpvallw <- rep(0,sum(subw))

for(i in 1:10000){

	sam <- resam(125,161,c(gM,gF))
	rgdiff <- rowSums(scow[,sam[[1]]]) - rowSums(scow[,sam[[2]]])
	gpvaluw <- gpvaluw + (rgdiff <= gdiffw)
	gpvallw <- gpvallw + (rgdiff >= gdiffw)
	if(i %% 100 == 0){print(i)}

}

# turn counts into p-values
gpvaluw <- gpvaluw / 10000
gpvallw <- gpvallw / 10000

# place refined p-values into original p-value vectors

gpvalu[subw] <- gpvaluw
gpvall[subw] <- gpvallw


# plot windows with low p-values

plot(log(rowMeans(sco[gpvall < 0.01,gM])/rowMeans(sco[gpvall < 0.01,gF]),2),pch=20,cex=.2)
plot(log(rowMeans(sco[gpvalu < 0.01,gM])/rowMeans(sco[gpvalu < 0.01,gF]),2),pch=20,cex=.2)

plot(log(gpvall,10))


# collapse into outlier windows
upwin <- mergewin(win=sco[gpvall < 0.001,2:4],stat=gpvalu[gpvall < 0.001],qu=0.001,buff=5000)
lowin <- mergewin(win=sco[gpvalu < 0.001,2:4],stat=gpvalu[gpvalu < 0.001],qu=0.001,buff=5000)

upwin <- upwin[order(upwin$count,decreasing=TRUE),]
lowin <- lowin[order(lowin$count,decreasing=TRUE),]

subscaf <- sco[,2] == "NW_012224686.1"
plot(log(rowMeans(sco[subscaf,gM])/rowMeans(sco[subscaf,gF]),2),pch=20,cex=1,col=(gpvalu[subscaf] < 0.001 | gpvall[subscaf] < 0.001)+1)

bigdups <- (gpvall < 0.001) & (sco[,2] %in% upwin[1:7,1])
bigdups <- colSums(sco[bigdups,c(gM,gF)])
MF <- c(rep(1,length(gM)),rep(2,length(gF)))

ord <- order(bigdups)
plot(bigdups[ord],pch=20,col=MF[ord])


#########################################

#########################################

#########################################

# permutations for HETEROCLITUS:
# iteratively refine low p-values to reduce computational burden

# empirical difference in coverage between sexes. 
	# should matter that I use sum instead of mean b/c permutations
hdiff <- rowSums(sco[,hM]) - rowSums(sco[,hF])


# upper and lower tail p-value initialization
hpvalu <- rep(0,length(hdiff))
hpvall <- rep(0,length(hdiff))

# loop through 100 randomizations
for(i in 1:100){

	sam <- resam(128,155,c(hM,hF))
	rhdiff <- rowSums(sco[,sam[[1]]]) - rowSums(sco[,sam[[2]]])
	hpvalu <- hpvalu + (rhdiff <= hdiff)
	hpvall <- hpvall + (rhdiff >= hdiff)
	if(i %% 10 == 0){print(i)}

}

# turn counts into p-values
hpvalu <- hpvalu/100
hpvall <- hpvall/100


# refine p-values for subset with p < 0.1
	# do 1k loops for these windows
subw <- hpvalu < 0.1 | hpvall < 0.1
scow <- sco[subw,]

hdiffw <- hdiff[subw]
hpvaluw <- rep(0,sum(subw))
hpvallw <- rep(0,sum(subw))

for(i in 1:1000){

	sam <- resam(128,155,c(gM,gF))
	rhdiff <- rowSums(scow[,sam[[1]]]) - rowSums(scow[,sam[[2]]])
	hpvaluw <- hpvaluw + (rhdiff <= hdiffw)
	hpvallw <- hpvallw + (rhdiff >= hdiffw)
	if(i %% 20 == 0){print(i)}

}

# turn counts into p-values
hpvaluw <- hpvaluw/1000
hpvallw <- hpvallw/1000

# place refined p-values into original p-value vectors

hpvalu[subw] <- hpvaluw
hpvall[subw] <- hpvallw

# iterate, this time refine p-values < 0.02
	# do 10k permutations

subw <- hpvalu < 0.02 | hpvall < 0.02
scow <- sco[subw,]

hdiffw <- hdiff[subw]
hpvaluw <- rep(0,sum(subw))
hpvallw <- rep(0,sum(subw))

for(i in 1:10000){

	sam <- resam(128,155,c(gM,gF))
	rhdiff <- rowSums(scow[,sam[[1]]]) - rowSums(scow[,sam[[2]]])
	hpvaluw <- hpvaluw + (rhdiff <= hdiffw)
	hpvallw <- hpvallw + (rhdiff >= hdiffw)
	if(i %% 100 == 0){print(i)}

}

# turn counts into p-values
hpvaluw <- hpvaluw / 10000
hpvallw <- hpvallw / 10000

# place refined p-values into original p-value vectors

hpvalu[subw] <- hpvaluw
hpvall[subw] <- hpvallw


# # plot windows with low p-values

# plot(log(rowMeans(sco[gpvall < 0.01,gM])/rowMeans(sco[gpvall < 0.01,gF]),2),pch=20,cex=.2)
# plot(log(rowMeans(sco[gpvalu < 0.01,gM])/rowMeans(sco[gpvalu < 0.01,gF]),2),pch=20,cex=.2)

# plot(log(gpvall,10))


# collapse into outlier windows
hupwin <- mergewin(win=sco[hpvall < 0.001,2:4],stat=hpvalu[hpvall < 0.001],qu=0.001,buff=10000)
hlowin <- mergewin(win=sco[hpvalu < 0.001,2:4],stat=hpvalu[hpvalu < 0.001],qu=0.001,buff=10000)

hupwin <- hupwin[order(hupwin$count,decreasing=TRUE),]
hlowin <- hlowin[order(hlowin$count,decreasing=TRUE),]

subscaf <- sco[,2] == "NW_012224765.1"
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM])/rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,gM])/rowMeans(sco[subscaf,gF]),2),pch=20,cex=1,col=(gpvalu[subscaf] < 0.005 | gpvall[subscaf] < 0.005)+1)

subscaf <- sco[,2] == "NW_012224765.1"
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM])/rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM])/rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)

subscaf <- sco[,2] == "NW_012225537.1"
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)
points(sco[subscaf,3],log(rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+3)
plot(sco[subscaf,3],log(rowMeans(sco[subscaf,hM]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)
points(sco[subscaf,3],log(rowMeans(sco[subscaf,hF]),2),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+3)

plot(sco[subscaf,3],rowMeans(sco[subscaf,hM]) - rowMeans(sco[subscaf,hF]),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)
plot(sco[subscaf,3],rowMeans(sco[subscaf,hM]) - rowMeans(sco[subscaf,hF]),pch=20,cex=1,col=(hpvalu[subscaf] < 0.005 | hpvall[subscaf] < 0.005)+1)


# bigdups <- (gpvall < 0.001) & (sco[,2] %in% upwin[1:7,1])
# bigdups <- colSums(sco[bigdups,c(gM,gF)])
# MF <- c(rep(1,length(gM)),rep(2,length(gF)))

# ord <- order(bigdups)
# plot(bigdups[ord],pch=20,col=MF[ord])

par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
plot(log(rowMeans(sco[,hM])/rowMeans(sco[,hF]),2),pch=20,cex=.3,col=(hpvalu < 0.005 | hpvall < 0.005)+1)
plot(log(rowMeans(sco[,gM])/rowMeans(sco[,gF]),2),pch=20,cex=.3,col=(gpvalu < 0.005 | gpvall < 0.005)+1)
