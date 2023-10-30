#-------------------------------------------------------------------------#
#                Composite minimax concave penalty (cMCP)                 #
#-------------------------------------------------------------------------#


#### Load library ####
library(survival)
library(grpreg) 

set.seed(17920) 
groups = c(rbind(c(1:500), c(1:500)),0)

cMCP <- function(data) {
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

  res.cMCP <-
    cv.grpsurv(
      XXTT,
      y,
      group = groups,
      nfolds = 5,
      penalty = "cMCP",
      se = c('quick'),
      returnX = TRUE,
      returnY = TRUE,
      trace = TRUE
    )
  
  if (length(which(res.cMCP$fit$beta[, res.cMCP$min] != 0)) == 0) {
    coef.cMCP <- "0"
    beta.cMCP <- 0
    lp <- rep(0, dim(data)[1])
  } else {
    coef.cMCP <-
      names(which(res.cMCP$fit$beta[, res.cMCP$min] != 0))
    beta.cMCP <- res.cMCP$fit$beta[coef.cMCP, res.cMCP$min]
    lp <- XXTT[, coef.cMCP, drop=FALSE] %*% beta.cMCP
  }
  
  mcox <-
    coxph(Surv(data$time, data$status) ~ offset(lp), data = data)
  
## baseline function
  basehaz.cMCP = basehaz(mcox)
  
  allrescMCP <- list(
    coef = coef.cMCP,
    beta = beta.cMCP,
    coef.cMCP = coef.cMCP,
    beta.cMCP = beta.cMCP,
    basehaz.cMCP =  basehaz.cMCP,
    biomarkers_active.main = biomarkers_active.main,
    biomarkers_active.inter = biomarkers_active.inter,
    lp = lp
  )
  
  return(allrescMCP)
}