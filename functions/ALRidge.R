#-------------------------------------------------------------------------#
#                    Adaptive lasso Ridge (Ridge)                         #
#-------------------------------------------------------------------------#

 
#### Load library ####
library(biospear)
library(survival) 

set.seed(17920)


cv.ALRidge <- function(data) {
  
  biomarkers_active.main = attributes(data)$biomarkers$active.main
  biomarkers_active.inter = attributes(data)$biomarkers$active.inter
  
  pos.bioi <- grep("bm", names(data))
  nbiom <- length(pos.bioi)

  X <- as.matrix(data[, pos.bioi])
  XT <- (X * data[, 'treat'])
  colnames(XT) <- paste0(colnames(XT), ".I")
  XXT <- matrix(rbind(X, XT), nrow = dim(X)[1])
  colnames(XXT) <- c(rbind(colnames(X), colnames(XT)))
  XXTT = cbind(XXT, data$treat)
  colnames(XXTT) <- c(colnames(XXT), "treat")
  
  print('Begin cvALR')
  cv <-
    BMsel(
      data = data,
      x = 2:501,
      y = 502:503,
      tt = 1,
      inter = TRUE,
      std.x = FALSE,
      std.i = FALSE,
      std.tt = FALSE,
      method = c('alassoR'),
      folds = 5,
      showWarn = TRUE, 
      trace = TRUE
    )
## print('End cvALR')
  
  coef.ALR = coef.ALRI = coef.T = "0"
  beta.ALR = beta.ALRI = beta.T = 0
  names(cv$alassoR) <- gsub(":treat", ".I", names(cv$alassoR))

  if (length(cv$alassoR) == 0) {
    coef.ALRI <- "0"
    beta.ALRI <- 0
    coef.ALR <- "0"
    beta.ALR <- 0
    coef.T <- "0"
    beta.T <- 0
    lp <- rep(0, dim(data)[1])
    
  } else {
    allcoef.ALR = names(cv$alassoR)
    allbeta.ALR = as.numeric(cv$alassoR)
    if(length(which(names(cv$alassoR) %in% "treat")) != 0){
      coef.T <- "treat"
      beta.T <- as.numeric(cv$alassoR[which(names(cv$alassoR) == "treat")])
    }
    if(length(grep("I", names(cv$alassoR))) != 0){
      coef.ALRI <- names(cv$alassoR)[grep("I", names(cv$alassoR))]
      coef.ALR <- setdiff(names(cv$alassoR), c("treat",coef.ALRI))
      beta.ALRI <- as.numeric(cv$alassoR[coef.ALRI])
      beta.ALR <- as.numeric(cv$alassoR[coef.ALR])
    }else{
      coef.ALRI <- "0"
      coef.ALR <- setdiff(names(cv$alassoR), c("treat",coef.ALRI))
      beta.ALRI <- 0
      beta.ALR <- as.numeric(cv$alassoR[coef.ALR])
    }
    
## print('End beta')
    
    allcoef.ALR<- gsub(":treat", ".I", allcoef.ALR)
    lp <- XXTT[, allcoef.ALR, drop = FALSE] %*% allbeta.ALR
    
## print('end lp')
  }
  
  mcox <-
    coxph(Surv(data$time, data$status) ~ offset(lp), data = data)
  
## baseline function
  basehaz.ALR = basehaz(mcox)
  allresALR <- list(
    coef = allcoef.ALR,
    beta = allbeta.ALR,
    coef.ALR = allcoef.ALR,
    beta.ALR = allbeta.ALR,
    basehaz.ALR =  basehaz.ALR,
    lp = lp,
    beta.T = beta.T,
    biomarkers_active.main = biomarkers_active.main,
    biomarkers_active.inter = biomarkers_active.inter,
    coef.ALR.I = coef.ALRI,
    beta.ALR.I = beta.ALRI,
    coef.ALR.m = coef.ALR,
    beta.ALR.m = beta.ALR
  )
  return(allresALR)
  
}
