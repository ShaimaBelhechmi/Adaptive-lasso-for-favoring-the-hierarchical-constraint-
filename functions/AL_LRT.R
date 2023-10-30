#-------------------------------------------------------------------------#
#              Adaptive lasso Likelihood ratio test (LRT)                 #
#-------------------------------------------------------------------------#


#### Load library ####

library(survival)
library(glmnet)

thresh <- 1e-16
nfolds <- 5
set.seed(17920)


AL <- function(data, nfolds, thresh) {
  if (missing(nfolds)) 
    nfolds <- 5
  else
    nfolds <- as.numeric(nfolds)
  
  if (missing(thresh))
    thresh <- 1e-16
  else
    thresh <- as.numeric(thresh)
  
  biomarkers_active.main = attributes(data)$biomarkers$active.main
  biomarkers_active.inter = attributes(data)$biomarkers$active.inter
  
  pos.bioi <- grep("bm", names(data))
  nbiom <- length(pos.bioi)
  Y <- cbind(time = data$time, status = data$status)
  X <- as.matrix(data[, pos.bioi])
  XT <- (X * data[, 'treat'] )
  colnames(XT) <- paste0(colnames(XT), ".I")
  XXT <- matrix(rbind(X, XT), nrow = dim(X)[1])
  colnames(XXT) <- c(rbind(colnames(X), colnames(XT)))
  XXTT = cbind(XXT, data$treat)
  colnames(XXTT) <- c(colnames(XXT), "treat")
  
  f.weights <- function(k) {
    m0 <- coxph(Surv(time, status) ~ treat, data = data)
    m1 <- coxph(Surv(time, status) ~ treat + X[, k], data = data)
    m2 <-
      coxph(Surv(time, status) ~ treat + X[, k] + XT[, k], data = data)
    wx <- (anova(m0, m2)$Chisq)[2]
    wxt <- (anova(m1, m2)$Chisq)[2]
    return(list(wx, wxt))
  }
  
  weights <- 1 / as.numeric(sapply(c(1:nbiom), f.weights))
  weights <- c(weights, 0)

  # Cross-validation to estimate the optimal lambda
  cv <- cv.glmnet(
    x = XXTT,
    y = Y,
    family = "cox",
    alpha = 1,
    nfolds = nfolds,
    grouped = TRUE,
    standardize = FALSE,
    thresh = thresh,
    penalty.factor = weights
  )
  
  # Fit of the final model
  fit.AL <- glmnet(
    x = XXTT,
    y = Y,
    family = "cox",
    alpha = 1,
    lambda = cv$lambda.min * ((nfolds - 1) / nfolds),
    standardize = FALSE,
    thresh = thresh,
    penalty.factor = weights
  )
  
  if (length(which(coef(fit.AL) != 0)) == 0) {
    coef.AL <- "0"
    beta.AL <- 0
    lp <- rep(0, dim(data)[1])
  } else{
    coef.AL <- coef(fit.AL)@Dimnames[[1]][which(coef(fit.AL) != 0)]
    beta.AL <- coef(fit.AL)@x
    lp <- XXTT[, coef.AL, drop=FALSE] %*% beta.AL
  }
  
  mcox <-
    coxph(Surv(data$time, data$status) ~ offset(lp), data = data)
  
  #baseline function
  basehaz.AL = basehaz(mcox)
  
  allresAL <- list(
    coef = coef.AL,
    beta = beta.AL,
    coef.AL = coef.AL,
    beta.AL = beta.AL,
    basehaz.AL =  basehaz.AL,
    biomarkers_active.main = biomarkers_active.main,
    biomarkers_active.inter = biomarkers_active.inter,
    lp = lp, 
    Weights = weights
  )
  
  return(allresAL)
  
}
