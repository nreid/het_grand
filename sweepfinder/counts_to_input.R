library(dplyr)

# srun -J interactive --pty bash
# module load r/3.3.3

chr <- list.files("/scratch/nmr15102/sweepfinder/counts",pattern="^chr.*gz",full.names=TRUE)
ugh <- read.table("~/fhet_genome/fst_dxy_allpops_liftover.txt", stringsAsFactors=FALSE)

for(k in 1:length(chr)){

# genomic windows to use for grid location selection
lift <- ugh
lord <- order(lift[,1],lift[,2])
lift <- lift[lord,]

# read in a chromosome
tab <- read.table(chr[k],stringsAsFactors=FALSE)
tab <- tab[order(tab[,2]),]

# make sure first and last variant in each dataset is the same. 
ss <- which(rowSums(tab[,seq(8,32,2)] > 5)==13 & rowSums(tab[,seq(7,32,2)] > 0)==13 & rowSums(tab[,seq(7,32,2)] < tab[,seq(8,32,2)])==13)
tab <- tab[min(ss):max(ss),]

# windows for grid points
subw <- lift[,1] == tab[1,1] & (lift[,2]+1) > tab[1,2] & lift[,3] < tab[dim(tab)[1],3]
lift <- lift[subw,]
gridp <- lift[,3]
# write grid file
cat(gridp,sep="\n",file=paste(lift[1,1],"grid",sep="."))

popcols <- c("chr","cstart","cend","scaf","start","end","BB.x","BB.n","BNP.x","BNP.n","BP.x","BP.n","ER.x","ER.n","F.x","F.n","GB.x","GB.n","KC.x","KC.n","NYC.x","NYC.n","PB.x","PB.n","SH.x","SH.n","SJSP.x","SJSP.n","SP.x","SP.n","VB.x","VB.n")
colnames(tab) <- popcols

grand <- c(7,9,17,23,27,29,31)
het <- c(11,13,15,19,21,25)

# for grandis populations
for(i in grand){
	#i <- 7
	s1 <- tab[,c(2,i,i+1)]
	ss <- s1[,3]>5 & s1[,2] > 0 & s1[,2] < s1[,3] 
	s1 <- s1[ss,]
	s1 <- cbind(s1,folded=0)
	colnames(s1) <- c("position","x","n","folded")
	rec1 <- cbind(pos=s1[,1],rate=(s1[,1]-s1[1,1])/1e6)
	
	# file names
	ffile <- paste(lift[1,1],gsub(".x","",colnames(tab)[i]),"freq",sep=".")
	# rfile <- paste(lift[1,1],gsub(".x","",colnames(tab)[i]),"rec",sep=".")
	
	# write input files
	write.table(s1,file=ffile,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
	# write.table(rec1,file=rfile,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}

# polarize spectrum assuming grandis major allele is ancestral
gf <- rowSums(tab[,grand])/rowSums(tab[,grand+1])
pol <- which(gf > 0.5)

# for heteroclitus populations
for(i in het){
	#i <- 7
	s1 <- tab[,c(2,i,i+1)]
	s1[pol,2] <- s1[pol,3] - s1[pol,2]
	ss <- s1[,3]>5 & s1[,2] > 0 & s1[,2] < s1[,3] 
	s1 <- s1[ss,]
	s1 <- cbind(s1,folded=0)
	colnames(s1) <- c("position","x","n","folded")
	rec1 <- cbind(pos=s1[,1],rate=(s1[,1]-s1[1,1])/1e6)
	
	# file names
	ffile <- paste(lift[1,1],gsub(".x","",colnames(tab)[i]),"freq",sep=".")
	# rfile <- paste(lift[1,1],gsub(".x","",colnames(tab)[i]),"rec",sep=".")
	
	# write input files
	write.table(s1,file=ffile,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
	# write.table(rec1,file=rfile,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}

print(k)

}



