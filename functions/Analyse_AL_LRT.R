#-------------------------------------------------------------------------#
#              Adaptive lasso Likelihood ratio test (LRT)                 #
#-------------------------------------------------------------------------#


#### Load library ####
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(corpcor) 
library(parallel) 

thresh <- 1e-16
nfolds <- 5
set.seed(17920)


## Number of PC-cores used
pCores <- 60

## Definition of the cluster via 'makeCluster' 
cl <- makeCluster(pCores, outfile ="OUTanalyseAL.txt")
r <- 500

## Change the PATH 
# v <- paste0('~/PATH/datas/scen',rep(1:4, each =r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')

analyse <- function(character){
  load(character)
  print('load')
  source("./functions/AL_LRT.R")
  respath <- gsub('(datas|data-)', 'resultsAL', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    res <- AL(data[1:1500,], nfolds = 5, thresh = 1e-16)          
    print('Finish')
    flush.console()
    save(res, file = respath)
    print('Saved')
    flush.console() 
  } 
}

#### Execution of the program ####
start <- Sys.time()
z <- clusterApplyLB(cl = cl, v, analyse)
end <- Sys.time()
t.parallel <- difftime(end, start, units = "sec")
stopCluster(cl)



