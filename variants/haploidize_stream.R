#!/usr/bin/env Rscript

# 0/0:1:1,0:1:40:0:0:0,-0.30103,-3.99568
# 3/3:1:0,0,0,1,0,0:0:0:0,0,1,0,0:0,0,37,0,0:-3.69783,-3.69783,-3.69783,-3.69783,-3.69783,-3.69783,-0.30103,-0.30103,-0.30103,0,-3.69783,-3.69783,-3.69783,-0.30103,-3.69783,-3.69783,-3.69783,-3.69783,-0.30103,-3.69783,-3.69783

library(stringr)
library(magrittr)


	nal <- function(x){
		x <- str_split(x,",") %>% unlist() %>% length()
		x + 1
		}

	sam<-function(x,len){
		if(grepl("^\\.",x)){return(".")}

		x <- str_split(x,":")[[1]][3] %>% 
			str_split(.,",") %>% 
			unlist() %>% 
			as.numeric() # pull allele depth TAG field
		x <- x!=0 # alleles present/absent
		if(sum(x)==0){return(".")} # no genotype if no alleles
		x <- x/sum(x) # probability vector
		x <- sample(x=len,size=1,prob=x) - 1
		return(x)
		}


f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  if(grepl("^#",line)){write(line,stdout());next()}

  line <- str_split(line,"\\t") %>% unlist()
  ncol <- length(line)
  nl <- nal(line[5])
  line[10:ncol] <- sapply(line[10:ncol],sam,len=nl)
  line <- paste(line,collapse="\t")
  write(line,stdout())
  # write(line, stderr())
  # process line
}
