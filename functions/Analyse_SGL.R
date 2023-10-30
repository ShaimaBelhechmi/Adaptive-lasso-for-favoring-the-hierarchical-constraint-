#-------------------------------------------------------------------------#
#                     Sparse Group Lasso (SGL)                            #
#-------------------------------------------------------------------------#


#### Load library ####
library(dplyr)
library(survival)
library(Matrix)
library(SGL)
library(corpcor)   
library(parallel) 

set.seed(17920)

## Number of PC-cores used 
pCores <- 70

## Definition of the cluster via 'makeCluster'
cl <- makeCluster(pCores, outfile ="OUTanalyseSGL.txt")
r <- 500  

## Change the PATH 
# v <- paste0('~/PATH/datas/scen',rep(1:4, each =r),'/simdata-',sprintf("%03d",c(1:r)),'.Rdata')

index <- c(rbind(c(1:500), c(1:500)),501)

analyse <- function(character){
  load(character)
  source("./functions/SGL.R")
  respath <- gsub('(datas|data-)', 'resultsSGL', character)
  
  if (isTRUE(file.exists(respath))){
    print('Found')
    flush.console()
  }else{
    res <- cv.SGL(data[1:1500,], index = index)          
    # print('Finish')
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
print(paste("t.parallel",t.parallel))
stopCluster(cl)




