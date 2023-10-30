#-------------------------------------------------------------------------#
#                     Group exponential lasso (gel)                       #
#-------------------------------------------------------------------------#


#### Load library ####
library(survival)
library(grpreg)

set.seed(17920)
groups = c(rbind(c(1:500), c(1:500)),0)

GEL <- function(data) {
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
  
  y = as.matrix(cbind(time = data$time, status = data$status))
  groups = c(rbind(c(1:500), c(1:500)),0)

  res.GEL <-
    cv.grpsurv(
      XXTT,
      y,
      group = groups,
      nfolds = 5,
      penalty = "gel",
      tau = 1 / 3,
      se = c('quick'),
      returnX = TRUE,
      returnY = TRUE,
      trace = TRUE
    )
  if (length(which(res.GEL$fit$beta[, res.GEL$min] != 0)) == 0) {
    coef.GEL <- "0"
    beta.GEL <- 0
    lp <- rep(0, dim(data)[1])
  } else {
    coef.GEL <- names(which(res.GEL$fit$beta[, res.GEL$min] != 0))
    beta.GEL <- res.GEL$fit$beta[coef.GEL, res.GEL$min]
    lp <- XXTT[, coef.GEL, drop=FALSE] %*% beta.GEL
  }
  
  mcox <-
    coxph(Surv(data$time, data$status) ~ offset(lp), data = data)
  
## baseline function
  basehaz.GEL = basehaz(mcox)
  
  allresGEL <- list(
    coef = coef.GEL,
    beta = beta.GEL,
    coef.GEL = coef.GEL,
    beta.GEL = beta.GEL,
    basehaz.GEL =  basehaz.GEL,
    biomarkers_active.main = biomarkers_active.main,
    biomarkers_active.inter = biomarkers_active.inter,
    lp = lp
  )
  
  return(allresGEL)
}