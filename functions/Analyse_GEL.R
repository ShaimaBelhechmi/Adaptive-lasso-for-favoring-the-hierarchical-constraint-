#-------------------------------------------------------------------------#
#                     Group exponential lasso (gel)                       #
#-------------------------------------------------------------------------#


#### Load library ####
library(parallel) 
library(survival)
library(grpreg) 

set.seed(17920)

## Number of PC-cores used
pCores <- 50 

## Definition of the cluster via 'makeCluster'
cl <- makeCluster(pCores, outfile ="OUTanalyseGEL.txt")
r <- 500

## Change the PATH 
# v <- paste0('~/PATH/datas/scen',rep(1:4, each =r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')

analyse <- function(character){
  load(character)
  print('load')
  source("./functions/GEL.R")
  respath <- gsub('(datas|data-)', 'resultsGEL', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    res <- GEL(data[1:1500,])          
    flush.console()
    save(res, file = respath)
    flush.console()
  } 
}

## Execution of the program
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)




