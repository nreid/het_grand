#!/usr/bin/env Rscript


library(stringr)
library(magrittr)
library(phytools)
library(tidyverse)
library(mscr)

# this is really dumb!

cname <- c(
	"CHROM",
	"POS",
	"ID",
	"REF",
	"ALT",
	"QUAL",
	"FILTER",
	"INFO",
	"FORMAT",
	"BP-10",
	"BP-11",
	"BP-12",
	"BP-13",
	"BP-14",
	"BP-15",
	"BP-16",
	"BP-17",
	"BP-18",
	"BP-19",
	"BP-1",
	"BP-20",
	"BP-21",
	"BP-22",
	"BP-23",
	"BP-24",
	"BP-25",
	"BP-26",
	"BP-27",
	"BP-28",
	"BP-29",
	"BP-2",
	"BP-30",
	"BP-31",
	"BP-32",
	"BP-33",
	"BP-34",
	"BP-38",
	"BP-39",
	"BP-3",
	"BP-40",
	"BP-41",
	"BP-42",
	"BP-43",
	"BP-44",
	"BP-45",
	"BP-46",
	"BP-47",
	"BP-48",
	"BP-49",
	"BP-4",
	"BP-50",
	"BP-51",
	"BP-52",
	"BP-53",
	"BP-5",
	"BP-6",
	"BP-7",
	"BP-8",
	"BP-9",
	"ER-11",
	"ER-12",
	"ER-13",
	"ER-14",
	"ER-15",
	"ER-16",
	"ER-17",
	"ER-18",
	"ER-19",
	"ER-20",
	"ER-21",
	"ER-22",
	"ER-23",
	"ER-24",
	"ER-25",
	"ER-26",
	"ER-27",
	"ER-28",
	"ER-29",
	"ER-30",
	"ER-31",
	"ER-32",
	"ER-33",
	"ER-34",
	"ER-35",
	"ER-36",
	"ER-38",
	"ER-39",
	"ER-40",
	"ER-41",
	"ER-42",
	"ER-43",
	"ER-44",
	"ER-45",
	"ER-46",
	"ER-47",
	"ER-48",
	"ER-49",
	"ER-50",
	"ER-51",
	"ER-52",
	"ER-53",
	"ER-54",
	"ER-55",
	"ER-56",
	"ER-57",
	"ER-58",
	"ER-59",
	"ER-60",
	"F-10",
	"F-11",
	"F-12",
	"F-13",
	"F-14",
	"F-15",
	"F-16",
	"F-17",
	"F-18",
	"F-19",
	"F-1",
	"F-20",
	"F-21",
	"F-22",
	"F-23",
	"F-24",
	"F-25",
	"F-26",
	"F-27",
	"F-28",
	"F-29",
	"F-2",
	"F-30",
	"F-31",
	"F-32",
	"F-33",
	"F-34",
	"F-35",
	"F-39",
	"F-40",
	"F-41",
	"F-42",
	"F-43",
	"F-44",
	"F-45",
	"F-46",
	"F-47",
	"F-49",
	"F-4",
	"F-50",
	"F-51",
	"F-52",
	"F-53",
	"F-54",
	"F-5",
	"F-6",
	"F-7",
	"F-8",
	"F-9",
	"KC-11",
	"KC-13",
	"KC-14",
	"KC-15",
	"KC-16",
	"KC-18",
	"KC-19",
	"KC-1",
	"KC-20",
	"KC-22",
	"KC-23",
	"KC-26",
	"KC-27",
	"KC-28",
	"KC-29",
	"KC-2",
	"KC-30",
	"KC-32",
	"KC-33",
	"KC-34",
	"KC-35",
	"KC-36",
	"KC-37",
	"KC-38",
	"KC-39",
	"KC-3",
	"KC-40",
	"KC-41",
	"KC-42",
	"KC-43",
	"KC-44",
	"KC-45",
	"KC-46",
	"KC-47",
	"KC-48",
	"KC-49",
	"KC-4",
	"KC-50",
	"KC-51",
	"KC-52",
	"KC-54",
	"KC-55",
	"KC-56",
	"KC-5",
	"KC-6",
	"KC-7",
	"KC-9",
	"NYC-10",
	"NYC-11",
	"NYC-12",
	"NYC-13",
	"NYC-14",
	"NYC-15",
	"NYC-16",
	"NYC-17",
	"NYC-18",
	"NYC-19",
	"NYC-20",
	"NYC-21",
	"NYC-22",
	"NYC-23",
	"NYC-24",
	"NYC-25",
	"NYC-26",
	"NYC-27",
	"NYC-28",
	"NYC-29",
	"NYC-30",
	"NYC-31",
	"NYC-32",
	"NYC-33",
	"NYC-34",
	"NYC-40",
	"NYC-41",
	"NYC-42",
	"NYC-43",
	"NYC-44",
	"NYC-45",
	"NYC-46",
	"NYC-47",
	"NYC-48",
	"NYC-49",
	"NYC-50",
	"NYC-51",
	"NYC-52",
	"NYC-53",
	"NYC-54",
	"NYC-55",
	"NYC-8",
	"NYC-9",
	"SH-14",
	"SH-15",
	"SH-16",
	"SH-17",
	"SH-18",
	"SH-19",
	"SH-1",
	"SH-201",
	"SH-202",
	"SH-203",
	"SH-204",
	"SH-205",
	"SH-206",
	"SH-207",
	"SH-208",
	"SH-209",
	"SH-20",
	"SH-210",
	"SH-211",
	"SH-212",
	"SH-213",
	"SH-21",
	"SH-22",
	"SH-23",
	"SH-24",
	"SH-25",
	"SH-26",
	"SH-27",
	"SH-28",
	"SH-29",
	"SH-2",
	"SH-30",
	"SH-31",
	"SH-32",
	"SH-33",
	"SH-34",
	"SH-35",
	"SH-36",
	"SH-37",
	"SH-38",
	"SH-39",
	"SH-3",
	"SH-40",
	"SH-41",
	"SH-42",
	"SH-4",
	"SH-5",
	"SH-6",
	"SH-7",
	"SH-8",
	"BU000004.VB_B",
	"BU000005.VB_B",
	"BU000006.VB_B",
	"BU000007.VB_B",
	"BU000008.VB_B",
	"BU000012.SP",
	"BU000014.SP",
	"BU000017.SP",
	"BU000018.SP",
	"BU000023.SP",
	"BU000024.SP",
	"BU000025.SP",
	"BU000031.SP",
	"BU000032.SP",
	"BU000033.SP",
	"BU000035.SP",
	"BU000036.SP",
	"BU000037.SP",
	"BU000039.SP",
	"BU000041.SP",
	"BU000046.SP",
	"BU000048.SP",
	"BU000049.SP",
	"BU000052.SP",
	"BU000053.SP",
	"BU000054.SP",
	"BU000055.SP",
	"BU000056.SP",
	"BU000057.SP",
	"BU000062.GB",
	"BU000063.GB",
	"BU000064.GB",
	"BU000065.GB",
	"BU000066.GB",
	"BU000067.GB",
	"BU000068.GB",
	"BU000069.GB",
	"BU000070.GB",
	"BU000071.GB",
	"BU000072.GB",
	"BU000073.GB",
	"BU000074.GB",
	"BU000075.GB",
	"BU000076.GB",
	"BU000077.GB",
	"BU000078.GB",
	"BU000081.GB",
	"BU000082.GB",
	"BU000083.GB",
	"BU000084.GB",
	"BU000085.GB",
	"BU000086.GB",
	"BU000087.GB",
	"BU000088.GB",
	"BU000089.GB",
	"BU000090.GB",
	"BU000092.GB",
	"BU000093.GB",
	"BU000094.GB",
	"BU000095.GB",
	"BU000097.GB",
	"BU000100.GB",
	"BU000101.GB",
	"BU000102.GB",
	"BU000103.GB",
	"BU000104.GB",
	"BU000105.GB",
	"BU000106.GB",
	"BU000110.GB",
	"BU000116.GB",
	"BU000120.GB",
	"BU000121.GB",
	"BU000123.GB",
	"BU000124.GB",
	"BU000125.GB",
	"BU000126.GB",
	"BU000127.GB",
	"BU000129.VB_A",
	"BU000130.VB_A",
	"BU000131.VB_A",
	"BU000132.VB_A",
	"BU000133.VB_A",
	"BU000134.VB_A",
	"BU000135.VB_A",
	"BU000136.VB_A",
	"BU000137.VB_A",
	"BU000138.VB_A",
	"BU000139.VB_A",
	"BU000140.VB_A",
	"BU000141.VB_A",
	"BU000142.VB_A",
	"BU000144.VB_A",
	"BU000145.VB_A",
	"BU000148.VB_A",
	"BU000149.VB_A",
	"BU000150.VB_A",
	"BU000151.VB_A",
	"BU000152.VB_A",
	"BU000153.VB_A",
	"BU000155.VB_A",
	"BU000157.VB_A",
	"BU000158.VB_A",
	"BU000160.VB_A",
	"BU000161.VB_A",
	"BU000164.VB_A",
	"BU000165.VB_A",
	"BU000166.VB_A",
	"BU000167.VB_A",
	"BU000168.VB_B",
	"BU000169.VB_B",
	"BU000170.VB_B",
	"BU000171.VB_B",
	"BU000172.VB_B",
	"BU000173.VB_B",
	"BU000174.VB_B",
	"BU000175.VB_B",
	"BU000176.VB_B",
	"BU000177.VB_B",
	"BU000178.VB_B",
	"BU000179.VB_B",
	"BU000180.VB_B",
	"BU000182.PB_A",
	"BU000183.PB_A",
	"BU000184.PB_A",
	"BU000185.PB_A",
	"BU000186.PB_A",
	"BU000187.PB_A",
	"BU000188.PB_A",
	"BU000190.PB_A",
	"BU000191.PB_A",
	"BU000192.PB_A",
	"BU000193.PB_A",
	"BU000194.PB_A",
	"BU000195.PB_A",
	"BU000196.PB_A",
	"BU000197.PB_A",
	"BU000198.PB_A",
	"BU000199.PB_A",
	"BU000200.PB_A",
	"BU000201.PB_A",
	"BU000202.PB_A",
	"BU000203.PB_A",
	"BU000204.PB_A",
	"BU000205.PB_A",
	"BU000206.PB_A",
	"BU000207.PB_B",
	"BU000209.PB_B",
	"BU000210.PB_B",
	"BU000211.PB_B",
	"BU000212.PB_B",
	"BU000213.PB_B",
	"BU000214.PB_B",
	"BU000215.PB_B",
	"BU000217.PB_B",
	"BU000219.PB_B",
	"BU000223.PB_B",
	"BU000225.PB_B",
	"BU000226.PB_B",
	"BU000227.PB_B",
	"BU000228.PB_B",
	"BU000229.PB_B",
	"BU000230.PB_B",
	"BU000231.PB_B",
	"BU000233.PB_B",
	"BU000234.PB_B",
	"BU000235.PB_B",
	"BU000237.PB_B",
	"BU000242.PB_B",
	"BU000244.SP",
	"BU000245.SP",
	"BU000246.SP",
	"BU000248.SP",
	"BU000249.SP",
	"BU000250.SP",
	"BU000252.SP",
	"BU000253.SP",
	"BU000254.SP",
	"BU000255.SP",
	"BU000256.SP",
	"BU000257.SP",
	"BU000259.SP",
	"BU000260.SP",
	"BU000261.SP",
	"BU000262.SP",
	"BU000263.SP",
	"BU000264.SP",
	"BU000265.SP",
	"BU000266.SP",
	"BU000269.SP",
	"BU000270.SP",
	"BU000271.SP",
	"BU000272.SP",
	"BU000318.BNP",
	"BU000319.BNP",
	"BU000320.BNP",
	"BU000321.BNP",
	"BU000322.BNP",
	"BU000323.BNP",
	"BU000324.BNP",
	"BU000325.BNP",
	"BU000326.BNP",
	"BU000327.BNP",
	"BU000328.BNP",
	"BU000329.BNP",
	"BU000330.BNP",
	"BU000331.BNP",
	"BU000332.BNP",
	"BU000333.BNP",
	"BU000334.BNP",
	"BU000335.BNP",
	"BU000336.BNP",
	"BU000337.BNP",
	"BU000338.BNP",
	"BU000339.BNP",
	"BU000340.BNP",
	"BU000341.BNP",
	"BU000343.BNP",
	"BU000344.BNP",
	"BU000345.BNP",
	"BU000346.BNP",
	"BU000347.BNP",
	"BU000348.BNP",
	"BU000349.BNP",
	"BU000350.BNP",
	"BU000351.BNP",
	"BU000352.BNP",
	"BU000354.BNP",
	"BU000355.BNP",
	"BU000356.BNP",
	"BU000357.BNP",
	"BU000358.BNP",
	"BU000359.BNP",
	"BU000361.BNP",
	"BU000362.BNP",
	"BU000364.BNP",
	"BU000366.BNP",
	"BU000367.BNP",
	"BU000372.BNP",
	"BU000375.BNP",
	"BU000382.BB",
	"BU000383.BB",
	"BU000384.BB",
	"BU000386.BB",
	"BU000390.BB",
	"BU000391.BB",
	"BU000392.BB",
	"BU000393.BNP",
	"BU000397.BB",
	"BU000398.BB",
	"BU000400.BB",
	"BU000402.BB",
	"BU000403.BB",
	"BU000405.BB",
	"BU000406.BB",
	"BU000407.BB",
	"BU000408.BB",
	"BU000409.BB",
	"BU000410.BB",
	"BU000411.BB",
	"BU000413.BB",
	"BU000414.BB",
	"BU000415.BB",
	"BU000416.BB",
	"BU000417.BB",
	"BU000418.SJSP",
	"BU000419.SJSP",
	"BU000420.SJSP",
	"BU000421.SJSP",
	"BU000424.SJSP",
	"BU000425.SJSP",
	"BU000426.SJSP",
	"BU000427.SJSP",
	"BU000431.SJSP",
	"BU000432.SJSP",
	"BU000433.SJSP",
	"BU000439.SJSP",
	"BU000441.SJSP",
	"BU000442.SJSP",
	"BU000443.SJSP",
	"BU000444.SJSP",
	"BU000446.SJSP",
	"BU000449.SJSP",
	"BU000451.SJSP",
	"BU000454.SJSP",
	"BU000458.SJSP",
	"BU000460.SJSP",
	"BU000463.SJSP",
	"BU000464.SJSP")

options(scipen=99)

cnamev <- cname


#feed this function a set of windows and corresponding statistics passing a threshold. 
	#it will merge them if they overlap or are within 'buffer' of each other and calculate stats. 
	#modified from a previous version in process_stats.R
	#
mergewin<-function(win, stat, qu, tails=TRUE, buff=5000){
	
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
	
	

# get vector of sample population names
pop <- cname[10:585] %>% gsub("-.*","",.) %>% gsub(".*\\.","",.) %>% gsub("_.*","",.)
popu <- unique(pop)

# get vector of north,admixed,south,grandis designations
sss <- pop
sss[grep("BP|F",sss)] <- "North"
sss[grep("NYC|SH",sss)] <- "Admix"
sss[grep("ER|KC",sss)] <- "South"
sss[grep("GB|BB|SJSP|SP|BNP|PB|VB",sss)] <- "Grand"
sssu <- unique(sss)

# 
lift <- read.table("~/projects/het_grand_data/fst_dxy_allpops_liftover.txt",stringsAsFactors=FALSE)
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

# old scaffold mappings
om <- read.table("~/projects/het_grand/old_scaffold_mappings.txt",stringsAsFactors=FALSE)
rownames(om) <- om[,1]
om2 <- om[-10180,]
rownames(om2) <- om2[,2]

# sexes
sexes <- read.table("../het_grand/all_sexes.txt",stringsAsFactors=FALSE)
rownames(sexes) <- sexes[,1]

# exclude individuals with lower coverage
toss <- c("BP-19","BP-26","BP-38","ER-20","ER-27","ER-36","ER-47","ER-50",
	"ER-53","F-23","F-35","F-39","F-47","F-50","F-52","KC-56",
	"KC-9","NYC-22","NYC-42","SH-32","SH-8","BU000101.GB","BU000136.VB_A","BU000144.VB_A",
	"BU000186.PB_A","BU000211.PB_B","BU000212.PB_B","BU000214.PB_B",
	"BU000217.PB_B","BU000228.PB_B","BU000242.PB_B","BU000349.BNP",
	"BU000415.BB","BU000419.SJSP","BU000421.SJSP","BU000432.SJSP",
	"BU000441.SJSP","BU000442.SJSP","BU000444.SJSP","BU000451.SJSP",
	"BU000458.SJSP","BU000460.SJSP","BU000463.SJSP")

cname <- c("chr","cstart","cend","scaf","start","end",cname[10:585])

keep <- which(!(cname[7:582] %in% toss) & sss=="Grand")
pop <- pop[keep]

cname <- c(cname[1:6],cname[keep+6])

# read in table of windowed distances from heteroclitus to grandis
ghd <- read.table("ghd.txt.gz",stringsAsFactors=FALSE)
colnames(ghd) <- cname
ord <- order(ghd[,1],ghd[,2])
ghd <- ghd[ord,]

# get only data columns
ghd2 <- ghd[,7:272]


# get some derivative stats
# each row scaled by the median
ghdm <- t(apply(ghd2,MAR=1,FUN=function(x){x/median(x,na.rm=TRUE)}))
# each column of above table scaled by the median
ghdmm <- apply(ghdm,MAR=2,FUN=function(x){x/median(x,na.rm=TRUE)})
# coefficient of variation for each row
cv <- apply(ghd2,MAR=1,FUN=function(x){sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)})
# mean for each row
me <- apply(ghd2,MAR=1,FUN=function(x){mean(x,na.rm=TRUE)})
# median for each row
med <- apply(ghd2,MAR=1,FUN=function(x){median(x,na.rm=TRUE)})
# 
q99 <- apply(ghd2,MAR=1,FUN=function(x){quantile(x,prob=0.9,na.rm=TRUE)})
#
q01 <- apply(ghd2,MAR=1,FUN=function(x){quantile(x,prob=0.01,na.rm=TRUE)})

ghdmin <- apply(ghd2,MAR=1,FUN=function(x){min(x,na.rm=TRUE)})


# some subset vectors excluding small windows, or unassigned scaffolds or particular chromosomes
subw <- (ghd[,3] - ghd[,2]) > 50000
subw2 <- ((ghd[,3] - ghd[,2]) > 50000) & (ghd[,1] != "chr1") & (ghd[,1] != "chr10") & (!grepl("^NW",ghd[,1]))
subw3 <- ((ghd[,3] - ghd[,2]) > 50000) & (!grepl("^NW",ghd[,1]))
subw4 <- ((ghd[,3] - ghd[,2]) > 50000) & (!grepl("^NW",ghd[,1])) & med > 30


# some plots
plot(me[subw2],pch=20,cex=.2,col=factor(ghd[subw2,1]))
plot(cv[subw2],pch=20,cex=.2,col=factor(ghd[subw2,1]))
plot(med[subw2],pch=20,cex=.2,col=factor(ghd[subw2,1]))
plot(ghdm[subw2,1],pch=20,cex=.2,col=factor(ghd[subw2,1]))



# plot a single indiviual's stats
par(mfrow=c(2,1),mar=rep(0,4),oma=c(3,3,1,1))
i <- "BU000201.PB_A"
# median scaled
plot(ghdm[subw3,i],pch=20,cex=.5,col=factor(ghd[subw3,1]))
# raw distances
plot(ghd[subw3,colnames(ghd)==i],pch=20,cex=.5,col=factor(ghd[subw3,1]))

thresh <- 0.93
bc <- c()
blist <- list()
for(i in 1:dim(ghdm)[2]){
	#i <- 4
	su <- which(subw3&ghdm[,i]<thresh)
	if(length(su) < 1){
		bc <- c(bc,0)
		blist[[i]] <- data.frame()
		next()
	}
	if(length(su) < 2){
		# double the outlier window just to get it through the function, which expects more than 1 window. 
		block <- mergewin(ghd[c(su,su),c(1:3)],ghdm[c(su,su),i],thresh,tails="lesser",buff=400001)
		# replace the window count to make it accurate. 
		block$count <- 1
		blist[[i]] <- cbind(block,ind=i,sam=colnames(ghdm)[i])
		bc <- c(bc,dim(block)[1])
		next()
	}
	block <- mergewin(ghd[su,c(1:3)],ghdm[su,i],thresh,tails="lesser",buff=400001)
	blist[[i]] <- cbind(block,ind=i,sam=colnames(ghdm)[i])
	bc <- c(bc,dim(block)[1])

	if((i%%10) == 0){print(i)}
}

# n introgression blocks, filtered by min count of windows
bc <- sapply(blist,FUN=function(x){dim(x[x$count > 0,])[1]})

# exclude the weird individual with subp
subp <- (1:length(bc)) != 177
# put all blocks in one table
allblocks <- do.call(rbind,blist[subp])[,-c(5,6)]
ord <- order(allblocks[,1],allblocks[,2],allblocks[,3])
allblocks <- allblocks[ord,]
# are distant blocks being strung together?
hist(allblocks$length/allblocks$count/100000)

# summaries by chromosome
# chr length
chrlen <- lift[grep("chr",lift[,1]),] %>% group_by(.,V1) %>% summarize(.,max(V3)) %>% data.frame()
colnames(chrlen)[1] <- "chrom"
sumblock <- group_by(allblocks[allblocks$count > 4,],chrom) %>% summarize(.,sum(length),length(length)) %>% data.frame()
blockstat <- full_join(chrlen,sumblock)

# does the SDR show signs of introgression?
data(sexchrs)
sexscaf <- om[sexchrs,2]
SDRint <- (ghd[,4] %in% sexscaf) & (rowSums(ghdm < thresh,na.rm=TRUE) > 0) & subw3
rowSums(ghdm[SDRint,subp] < thresh,na.rm=TRUE) %>% plot()

# do sexes show a biased signal?
# consider only BB,VB,PB
subi <- grep("BB|VB|PB",colnames(ghdm))
boxplot(colSums(ghdm[SDRint,subi] < thresh,na.rm=TRUE) ~ sexes[colnames(ghdm[,subi]),3])
barplot(table(colSums(ghdm[SDRint,subi] < thresh,na.rm=TRUE),sexes[colnames(ghdm[,subi]),3]),beside=TRUE)
barplot(table(colSums(ghdm[,subp] < thresh,na.rm=TRUE),sexes[colnames(ghdm[,subp]),3]),beside=TRUE)

# are G scores lower in males in SDR?
gsdr <- ghdm[SDRint,subi]
gsdr[gsdr >= thresh] <- NA
vals <- as.vector(gsdr)
sex <- rep(sexes[colnames(ghdm[,subi]),3],97)

boxplot(as.vector(gsdr)~rep(sexes[colnames(ghdm[,subi]),3],97))

wilcox.test(vals[sex=="M"],vals[sex=="F"])

# merge blocks across individuals
minblocks <- allblocks$count > 4
mblocks <- mergewin(allblocks[minblocks,1:3],stat=rep(0,dim(allblocks[minblocks,])[1]),qu=2,tail="lesser",buff=110)


# visualize blocks per individual per population
# filter by window count first:
bc <- sapply(blist,FUN=function(x){dim(x[x$count > 0,])[1]})

# filter by window count and chromosome:
# bc <- sapply(blist,FUN=function(x){dim(x[x$count > 0 & x$chrom=="chr24",])[1]})

upop <- c("BB","VB","PB","SJSP","BNP","SP","GB")
boxplot(bc[subp] ~ factor(pop[subp],upop))
barplot(table(bc[subp],pop[subp]),beside=TRUE)
table(bc[subp],pop[subp])
(group_by(data.frame(bc=bc[subp],pop=pop[subp]),pop) %>% summarize(.,med=median(bc),m=mean(bc),s=sum(bc),l=length(bc)))[c(1,7,4,5,2,6,3),]

# blocks per individual excluding chr1
cexchr1 <- allblocks[allblocks[,1]!="chr1","ind"] %>% factor(.,levels=1:266) %>% table()
table(pop,cexchr1)
barplot(table((cexchr1),pop),beside=TRUE)

# plot statistic distribution
hist(ghdm[subw3,subp],breaks=seq(0,1.4,0.001),xlim=c(0.7,1),main="G statistic distribution",xlab=NULL)
abline(v=0.93,col="red")

# min median scaled statistic for resistant and sensitive individuals
par(mfrow=c(2,1),mar=c(0,4,0,0),oma=c(3,3,1,1))
sen <- colnames(ghdm)!="BU000078.GB" & (pop=="SP"|pop=="GB")
plot(apply(ghdm[subw4,subp & !sen],MAR=1,FUN=min,na.rm=TRUE),pch=20,col=factor(ghd[subw4,1]),cex=.5,ylim=c(.3,1),ylab="Gmin - contam")
plot(apply(ghdm[subw4,subp & sen],MAR=1,FUN=min,na.rm=TRUE),pch=20,col=factor(ghd[subw4,1]),cex=.5,ylim=c(.3,1),ylab="Gmin - clean")


# estimate admixture time filtering blocks variously

1/mean(allblocks$length[allblocks$count > 4] * 1e-8)
1/mean(allblocks$length[allblocks$count > 0] * 1e-8)

bc <- sapply(blist,FUN=function(x){dim(x[x$count > 0,])[1]})
bcp <- data.frame(bc=bc[subp],pop=pop[subp]) %>% group_by(.,pop) %>% summarize(., mean(bc))
bcp <- bcp[[2]][c(1,7,4,5,2,6,3)]
blp <- c()
for(i in upop){
	blp <- c(blp,do.call(rbind, blist[which(pop==i)])$length %>% mean())

}

# population mean independent block length:
# the idea here is to treat the reference sites as representing "background introgression and error"
# then treat the polluted sites as being an average of recent introgression, background introgression and error
# we can use the count of blocks per individual and the mean length in a reference site to 
# determine the contribution of "recent introgression" in the polluted sites to the mean block length
# then use 1/mean block size in morgans as an estimate of the timing of introgression (from Graham)

PBblocklength <- (blp[3] * bcp[3])-(blp[6] * bcp[6]) / (bcp[3] - bcp[6])

1/(PBblocklength * 1e-8)

BBblocklength <- (blp[1] * bcp[1])-(blp[6] * bcp[6]) / (bcp[1] - bcp[6])

1/(BBblocklength * 1e-8)



# validate a block from individual 99, chr15: NW_012224452.1.vcf.gz


vcf <- read.table("NW_012224452.1.vcf.gz",stringsAsFactors=FALSE)
colnames(vcf) <- cnamev
keep <- !grepl(",",vcf[,5])
vcf <- vcf[keep,]
gt <- as.matrix(vcf[,10:585])
class(gt) <- "numeric"

keep <- which(!colnames(gt) %in% toss)
vcf <- vcf[,c(1:9,keep+9)]
gt <- gt[,keep]

d <- t(gt[1:38000,]) %>% dist(.,method="manhattan")


spec2 <- sss[keep]
tr <- nj(d) #%>% midpoint.root()
mds <- cmdscale(d)

plot(tr,"unrooted",show.tip.label=FALSE)
tiplabels(pch=20,col=factor(spec2),bg="white")
plot(mds,pch=20,col=factor(spec2))
text(mds[,1],mds[,2],rownames(mds),cex=.5)

plot(as.matrix(d)[444,-444]/as.matrix(d)[445,-444])

diffs <- which(rowMeans(gt[,spec2=="Grand"],na.rm=TRUE) - rowMeans(gt[,spec2!="Grand"],na.rm=TRUE) > 0.8)


# plot blocks

chr1 <- allblocks[allblocks[,1]=="chr1",-c(5,6)]
chr1 <- chr1[order(chr1$length,decreasing=TRUE),]


dev.new()
plot(NULL,xlim=c())

i <- 1
vec <- rep(NA,length(chr1[,1]))
vec[1] <- 1
while(any(is.na(vec))){
	for(j in 1:length(chr1[,1])){

		if(!is.na(vec[j])){next()}
		ov <- chr1[j,1] == chr1[which(vec==i),1] & (chr1[j,3] < chr1[which(vec==i),2] | chr1[j,2] > chr1[which(vec==i),3] )
		if(sum(!ov,na.rm=TRUE)){vec[j] <- i}
#		print(c(i,j))
		}
	i <- i + 1

}






