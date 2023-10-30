#-------------------------------------------------------------------------#
#                     Sparse Group Lasso (SGL)                            #
#-------------------------------------------------------------------------#


#### Load library ####
library(survival)
library(SGL)

set.seed(17920)

index <- c(rbind(c(1:500), c(1:500)),501)

cv.SGL <- function(data, index) {
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
  datas <- list(x = XXTT,
                time = data$time,
                status = data$status)
  print('Begin cvSGL')
  cv <-
    cvSGL(
      data = datas,
      index = index,
      type = "cox",
      maxit = 1000,
      thresh = 0.001,
      min.frac = 0.05,
      nlam = 20,
      gamma = 0.8,
      nfold = 5,
      standardize = FALSE,
      verbose = FALSE,
      step = 1,
      reset = 10,
      alpha = 0.95,
      lambdas = NULL
    )
  
  all.coef = coef.SGL = coef.m = coef.I = "0"
  beta.SGL = 0
  lp <- rep(0, dim(data)[1])
  if (length(which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0))) == 0) {

    coef.SGL =  coef.m = coef.I = "0"
    beta.SGL <- 0
    lp <- rep(0, dim(data)[1])
  } else {
    all.coef <-
      which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0))

    coef.SGL <- colnames(XXTT)[all.coef]
    if(sum(grep("I",coef.SGL)) != 0){
      coef.I <- coef.SGL[grep("I",coef.SGL)] 
      coef.m <- coef.SGL[-pmatch(coef.I, coef.SGL)]
    }else{
      coef.I <- "0"
      coef.m <- coef.SGL
    }
    beta.SGL <-
      cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))][which((cv$fit$beta[, which(cv$lldiff == min(cv$lldiff))] != 0))]
    lp <- XXTT[, coef.SGL, drop = FALSE] %*% beta.SGL
  }
  
  mcox <-
    coxph(Surv(data$time, data$status) ~ offset(lp), data = data)
  
## baseline function
  basehaz.SGL = basehaz(mcox)
  allressgl <- list(
    coef = coef.SGL,
    beta = beta.SGL,
    coef.SGL = coef.SGL,
    beta.SGL = beta.SGL,
    basehaz.SGL =  basehaz.SGL,
    biomarkers_active.main = biomarkers_active.main,
    biomarkers_active.inter = biomarkers_active.inter,
    lp,
    coef.SGL.m = coef.m,
    coef.SGL.I = coef.I
  )
  return(allressgl)
  
} 