cc <- read.table("chunkeval_out.txt",stringsAsFactors=FALSE,fill=TRUE,header=FALSE,row.names=NULL)

cc[,1] <- gsub(".vcf.gz","",cc[,1])

cc[,1] <- gsub(".*subvars/","",cc[,1]) %>% gsub("NW_","NW-",.) %>% gsub("NC_","NC-",.)

cc2 <- str_split(cc[,1],regex("_")) %>% do.call(rbind,.)
cc2[,1] <- gsub("-","_",cc2[,1])

cc3 <- cbind(cc2,cc[,3:4],stringsAsFactors=FALSE)
cc3[,2] <- as.numeric(cc3[,2])
cc3[,3] <- as.numeric(cc3[,3])

colnames(cc3) <- c("scaf","start","end","first","last")

cc3 <- cbind(cc3,diff=cc3[,"end"]-cc3[,"last"]+cc3[,"first"]-cc3[,"start"])
cc3 <- cbind(cc3,front=cc3[,"first"]-cc3[,"start"])
cc3 <- cbind(cc3,back=cc3[,"end"]-cc3[,"last"])

cc3[is.na(cc3[,"diff"]),] <- cc3[is.na(cc3[,"diff"]),"end"] - cc3[is.na(cc3[,"diff"]),"start"]

orda <- order(cc3[,"diff"])
ordf <- order(cc3[,"front"])
ordb <- order(cc3[,"back"])

val <- tail(cc3[ordb,c(1,5,3)],n=242)
val <- val[grep("NW",val[,1]),]

write.table(val,file="check_for_N.bed",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

# check base composition using BEDTOOLS

# BED=~/bin/bedtools2/bin/bedtools

# $BED nuc ...

nn <- read.table("check_for_N_eval.bed",stringsAsFactors=FALSE,header=TRUE,comment.char="")
colnames(nn) <- gsub("X.*[0-9]_","",colnames(nn))


plot(nn[,"seq_len"],nn[,"num_N"])

# problem intervals. recall. 
# NW_012234324.1   1200000   1300000   
# NW_012234324.1   1100000   1200000  
