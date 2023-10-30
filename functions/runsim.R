#-----------------------------------------------------------------------#
# Generation of training sets for the simulation study                  #
#-----------------------------------------------------------------------#
# The simulations should contains three studies:                        #
# (1): simulation with p = 500 biomarkers                               #
# (2): simulation with n = 3000 biomarkers                              #
# (3): more realistic setting (higher censor rates, lower effect sizes) #
#                                                                       #
# For each study, six scenarios have been implemented: three null and   #
# three alternative.                                                    #
# (S1): Complete null scenario                                          #
# (S2): One treatment-effect modifier                                   #
# (S3): Multiple treatment-effect modifiers                             #
# (S4): Multiple treatment-effect modifiers and prognostic biomarkers   #
#                                                                       #
# For more information, please refer to the section Simulation study    #
# of the Belhechmi et al. BMC Bioinformatics.                           #      
# https://doi.org/10.1186/s12859-023-05162-x                            #
#-----------------------------------------------------------------------#

#---------------------------------------------#
####            Data generation            ####
#---------------------------------------------#

library(optimx)
library(MASS)
library(biospear)
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(corpcor)
library(parallel) 
#A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5)
set.seed(12345)
r <- 500

# v <- paste0('~/path',sprintf("%03d",c(1:r)),'.Rdata')

character = v
respath = v
h <- 500
h.active <- 10
b.corr.by <- 25
a <- b1 <- b2 <- log(0.5) 
fu <- 2
recr <- 3

if (isTRUE(file.exists(respath))){
  print('Found')
  flush.console()
}else{
  for(i in 1:r){
    # Scenario 1: Complete null scenario
    data <- simdata(
      n = 3000,              # the sample size = number of patients
      p = h,                 # the number of biomarkers
      q.main = 0,            # change the number of true prognostic biomarkers
      q.inter = 0,           # change the number of true prognostic biomarkers
      prob.tt = 0.5,         # the treatement assignement probability
      m0 = 1,                # the baseline median survival time fixed for each arm
      alpha.tt = 0,          # change the effect of the treatment (in log-scale)
      beta.main = 0,         # change the effect of the prognostic biomarkers (in log-scale)
      beta.inter = 0,        # change the effect of the biomarkers interacting with the treatment (in log-scale)
      b.corr = 0.7,          # the correlation between biomarker blocks
      b.corr.by = b.corr.by, # change the size of the blocks of correlated biomarkers
      wei.shape = 1,         # the shape parameter of the Weibull distribution
      recr = recr,           # the recruitment period duration
      fu = fu,               # the follow-up period duration
      timefactor = 1         # the scale multiplicative factor for times (i.e. 1 = times in years)
    )
    save(data, file = paste0("./datas/scen1/simdata-",sprintf("%03d",i), ".Rdata"))
  }
}


respath <- gsub('scen1', 'scen2', character)
if (isTRUE(file.exists(respath))){
  print('Found')
  flush.console()
}else{
  for(i in 1:r){
    # Scenario 2: One treatment-effect modifier (q.inter = 1, beta.inter = log(0.5))
    data <- simdata(
      n = 3000,
      p = h,
      q.main = 0, 
      q.inter = 1, 
      prob.tt = 0.5,
      m0 = 1,
      alpha.tt = 0, 
      beta.main = 0, 
      beta.inter = b2, 
      b.corr = 0.7,
      b.corr.by = b.corr.by, 
      wei.shape = 1,
      recr = recr,
      fu = fu,
      timefactor = 1
    )
    save(data, file = paste0("./datas/scen2/simdata-",sprintf("%03d",i), ".Rdata"))
  }
}

respath <- gsub('scen1', 'scen3', character)
if (isTRUE(file.exists(respath))){
  print('Found')
  flush.console()
}else{    
  for(i in 1:r){
    # Scenario 3: Multiple treatment-effect modifiers (q.inter = h.active, beta.inter = log(0.5))
    data <- simdata(
      n = 3000,
      p = h,
      q.main = 0, 
      q.inter = h.active, 
      prob.tt = 0.5,
      m0 = 1,
      alpha.tt = 0, 
      beta.main = 0, 
      beta.inter = b2, 
      b.corr = 0.7,
      b.corr.by = b.corr.by, 
      wei.shape = 1,
      recr = recr,
      fu = fu,
      timefactor = 1
    )
    save(data, file = paste0("./datas/scen3/simdata-",sprintf("%03d",i), ".Rdata"))
  }
}

respath <- gsub('scen1', 'scen4', character)
if (isTRUE(file.exists(respath))){
  print('Found')
  flush.console()
}else{  
  for(i in 1:r){
    # Scenario 4: Multiple treatment-effect modifiers and prognostic biomarkers (q.main = q.inter = h.active, beta.main = beta.inter = log(0.5))
    data <- simdata(
      n = 3000,
      p = h,
      q.main = h.active, 
      q.inter = h.active, 
      prob.tt = 0.5,
      m0 = 1,
      alpha.tt = 0, 
      beta.main = b1, 
      beta.inter = b2, 
      b.corr = 0.7,
      b.corr.by = b.corr.by, 
      wei.shape = 1,
      recr = recr,
      fu = fu,
      timefactor = 1
    )
    save(data, file = paste0("./datas/scen4/simdata-",sprintf("%03d",i), ".Rdata"))
  }
}
