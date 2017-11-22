library(dplyr)
library(stringr)


options(scipen=99)

f <- file("stdin")
open(f)
ind <- 1

minsam <- c(3,3,5,5,5,3,5,5,3,5,3,3,3)

while(length(line <- readLines(f,n=1)) > 0) {

	# skip comment lines
	if(grepl("^#",line)){next()}
	ind <- ind + 1
	
	line <- str_split(line,"\\t") %>% unlist()
	samsize <- as.numeric(line[3:length(line)])

	if(!all(samsize>=minsam)){next()}

	cat(paste(c(line[1],as.numeric(line[2])-1,line[2],samsize),collapse="\t",sep=""),"\n",sep="")

	if((ind %% 100000) == 0){write(ind,stderr())}

}
