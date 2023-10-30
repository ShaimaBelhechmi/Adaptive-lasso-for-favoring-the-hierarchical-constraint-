#-------------------------------------------------------------------------#
#                    Adaptive lasso Ridge (Ridge)                         #
#-------------------------------------------------------------------------#

 
#### Load library ####
library(biospear)
library(survival) 
library(dplyr)
library(Matrix)
library(corpcor)   
library(parallel) 


set.seed(17920)

## Number of PC-cores used
pCores <- 60 

## Definition of the cluster via 'makeCluster'
cl <- makeCluster(pCores, outfile ="OUTanalyseALR.txt")
r <- 500 

## Change the PATH 
# v <- paste0('~/PATH/datas/scen',rep(1:4, each =r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')


analyse <- function(character){
  load(character)
  print('load')
  print(character)
  source("./functions/ALRidge.R")
  respath <- gsub('(datas|data-)', 'resultsALR', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    res <- cv.ALRidge(data[1:1500,])
    print(res)
    print('Finish')
    flush.console()
    save(res, file = respath)
    print('Saved')
    flush.console()
  } 
}

## Execution of the program
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
print(paste("t.parallel",t.parallel))
stopCluster(cl)




