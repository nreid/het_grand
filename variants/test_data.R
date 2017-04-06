# this is really dumb!
cname <- c("CHROM",
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


# read in 400kb of "haploidized" data. 
tab <- read.table("het_grand_data/test.out",stringsAsFactors=FALSE)
colnames(tab) <- cname

# toss any sites with more than 2 alleles
subs <- !grepl(",",tab[,5])
tab2 <- tab[subs,]

# convert haploid genotypes to numeric
tab3 <- as.matrix(tab2[,10:585])
class(tab3) <- "numeric"

# get vector of sample population names
pop <- colnames(tab3) %>% gsub("-.*","",.) %>% gsub(".*\\.","",.) %>% gsub("_.*","",.)

# get vector of north,admixed,south,grandis designations
sss <- pop
sss[grep("BP|F",sss)] <- "North"
sss[grep("NYC|SH",sss)] <- "Admix"
sss[grep("ER|KC",sss)] <- "South"
sss[grep("GB|BB|SJSP|SP|BNP|PB|VB",sss)] <- "Grand"

# plot missing data rate per individual
colMeans(is.na(tab3)) %>% plot(.,pch=20,col=factor(pop))

# toss individuals with greater than 90% missing data
subi <- colMeans(is.na(tab3)) < 0.9

# get ready to make some summary plots
library(ape)
dd <- t(tab3[,subi]) %>% dist() 
md <- cmdscale(dd)
tr <- nj(dd)

# MDS and tree plots
plot(md,col=factor(pop[subi]),pch=20)
plot(tr,"unrooted")

# calculate population allele frequencies
popfreqs <- c()

for(i in 1:13){

	popfreqs <- cbind(popfreqs,rowMeans(tab3[,pop==unique(pop)[i]],na.rm=TRUE))

	}

# calculate sample sizes of haploidized data for north south etc.
sssco <- c()
for(i in 1:4){

	sssco <- cbind(sssco,rowSums(!is.na(tab3[,sss==unique(sss)[i]]),na.rm=TRUE))

	}

# calculate frequencies for north, south, etc
sssfreqs <- c()
for(i in 1:4){

	sssfreqs <- cbind(sssfreqs,rowMeans(tab3[,sss==unique(sss)[i]],na.rm=TRUE))

	}


colnames(sssfreqs) <- unique(sss)

# exclude sites with ANY missing data at north,south etc. level
subs <- rowSums(is.na(sssfreqs))==0

# polarize variants by grandis (assume grandis major allele is ancestral)
sspol <- sssfreqs[subs,]
for(i in 1:dim(sspol)[1]){

	if(sspol[i,4]>0.5){
		sspol[i,] <- 1-sspol[i,]
	}

	}

# get ready to write some sweepfinder output. 
# fake recombination map, in case we need it. 1mb=1cM
basereco <- cbind(
	position=tab2[subs,2],
	rate=c(0,tab2[subs,2][-1]-tab2[subs,2][-sum(subs)])/1000000
	)

# format sweepfinder input tables, exclude sites with no derived alleles
northgeno <- cbind(position=tab2[subs,2],x=round(sspol[,1] * sssco[subs,1]),n=sssco[subs,1],folded=0)
northfgeno <- cbind(position=tab2[subs,2],x=round(sspol[,1] * sssco[subs,1]),n=sssco[subs,1],folded=1)
nn <- northgeno[,"x"] > 0
admixgeno <- cbind(position=tab2[subs,2],x=round(sspol[,2] * sssco[subs,2]),n=sssco[subs,2],folded=0)
an <- admixgeno[,"x"] > 0
southgeno <- cbind(position=tab2[subs,2],x=round(sspol[,3] * sssco[subs,3]),n=sssco[subs,3],folded=0)
sn <- southgeno[,"x"] > 0
grandgeno <- cbind(position=tab2[subs,2],x=round(sspol[,4] * sssco[subs,4]),n=sssco[subs,4],folded=1)
gn <- grandgeno[,"x"] > 0


# write out sweepfinder input files
options(scipen=999)

write.table(northgeno[nn,],file="het_grand_data/north.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(northfgeno[nn,],file="het_grand_data/north.fold.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(basereco[nn,],file="het_grand_data/north.reco",sep="\t",row.names=FALSE,quote=FALSE)

write.table(northgeno[nn,],file="het_grand_data/north.fold.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(basereco[nn,],file="het_grand_data/north.reco",sep="\t",row.names=FALSE,quote=FALSE)


write.table(admixgeno[an,],file="het_grand_data/admix.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(basereco[an,],file="het_grand_data/admix.reco",sep="\t",row.names=FALSE,quote=FALSE)

write.table(southgeno[sn,],file="het_grand_data/south.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(basereco[nn,],file="het_grand_data/south.reco",sep="\t",row.names=FALSE,quote=FALSE)

write.table(grandgeno[gn,],file="het_grand_data/grand.geno",sep="\t",row.names=FALSE,quote=FALSE)
write.table(basereco[nn,],file="het_grand_data/grand.reco",sep="\t",row.names=FALSE,quote=FALSE)


# run sweepfinder e.g. 
# ./SweepFinder2 -s 400 ~/projects/het_grand_data/grand.geno ~/projects/het_grand_data/grand.geno.out

# read in some sweepfinder results
ngof <- read.table("het_grand_data/north.fold.geno.out",stringsAsFactors=FALSE,header=TRUE)
ngo <- read.table("het_grand_data/north.geno.out",stringsAsFactors=FALSE,header=TRUE)
ngo2 <- read.table("het_grand_data/north.geno.out.3",stringsAsFactors=FALSE,header=TRUE)
ngo3 <- read.table("het_grand_data/north.geno.out.2",stringsAsFactors=FALSE,header=TRUE)

sgo <- read.table("het_grand_data/south.geno.out",stringsAsFactors=FALSE,header=TRUE)
sgo2 <- read.table("het_grand_data/south.geno.out.3",stringsAsFactors=FALSE,header=TRUE)
sgo3 <- read.table("het_grand_data/south.geno.out.2",stringsAsFactors=FALSE,header=TRUE)

ggo <- read.table("het_grand_data/grand.geno.out",stringsAsFactors=FALSE,header=TRUE)
ggo2 <- read.table("het_grand_data/grand.geno.out.3",stringsAsFactors=FALSE,header=TRUE)
ggo3 <- read.table("het_grand_data/grand.geno.out.2",stringsAsFactors=FALSE,header=TRUE)


# plot some sweepfinder results
par(mfrow=c(2,1))
plot(ngo[,1],ngo[,2],pch=20,ylim=c(0,23))
points(ngof[,1],ngof[,2],pch=20,col="red")
points(ngo2[,1],ngo2[,2],pch=20,col="green")
points(ngo3[,1],ngo3[,2],pch=20,col="blue")

plot(sgo[,1],sgo[,2],pch=20,ylim=c(0,23))
points(sgo2[,1],sgo2[,2],pch=20,col="red")
points(sgo3[,1],sgo3[,2],pch=20,col="blue")

plot(ggo[,1],ggo[,2],pch=20,ylim=c(0,23))
points(ggo2[,1],ggo2[,2],pch=20,col="red")
points(ggo3[,1],ggo3[,2],pch=20,col="blue")



