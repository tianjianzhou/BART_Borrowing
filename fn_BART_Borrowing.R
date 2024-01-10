# install packages first if not already installed
# install.packages("BART")
# install.packages("gtools")

library(BART)
library(gtools)

########################################################################
# Generate posterior samples of the ATE
########################################################################


calc_ATE_bart = function(data, burnin = 100, niter = 1000) {
  
  # DATA FORMAT
  # data: a list with y, X, s, t
  #    y: outcome vector of length n
  #       where n is number of total patients in current trial and external datasets
  #       y[i] is outcome for patient i (continuous outcome)
  #    X: covariate matrix of dimensions n * Q
  #       X[i, ] contains covariates for patient i
  #    s: data source vector of length n
  #       s[i] is indicator of data source for patient i 
  #       s[i] = 0 if patient i from current trial
  #       s[i] = 1, ..., J if patient i from j-th external dataset
  #    t: treatment vector of length n
  #       t[i] is indicator of treatment received for patient i
  #       t[i] = 0 if patient i received control
  #       t[i] = 1 if patient i received (new) treatment
  
  J = max(data$s)

  # Extract control arm data
  X_ctl = as.matrix(data$X[data$t == 0, ])
  s_ctl = factor(data$s[data$t == 0], levels = 0:J)
  y_ctl = data$y[data$t == 0]
  X_s_ctl = data.frame(X = X_ctl, s = s_ctl)
  
  # Extract treatment arm data
  X_trt = as.matrix(data$X[data$t == 1, ])
  y_trt = data$y[data$t == 1]
  
  # Extract covariates for patients in the current trial
  X_trial = as.matrix(data$X[data$s == 0, ])
  s_trial = factor(data$s[data$s == 0], levels = 0:J)
  X_s_trial = data.frame(X = X_trial, s = s_trial)
  # N_trial is number of patients in current trial
  N_trial = nrow(X_trial)

  # MCMC for control arm
  invisible(capture.output(
    MCMC_spls_ctl <- wbart(x.train = X_s_ctl, 
                           y.train = y_ctl,
                           x.test = X_s_trial,
                           ndpost = niter,
                           nskip = burnin,
                           printevery = 10000)))
  
  # MCMC for treatment arm
  invisible(capture.output(
    MCMC_spls_trt <- wbart(x.train = X_trt, 
                           y.train = y_trt,
                           x.test = X_trial,
                           ndpost = niter,
                           nskip = burnin,
                           printevery = 10000)))
  
  # Predict the outcomes of trial patients
  # CTE_spls: N_trial * niter matrix (after transpose)
  #           CTE_spls[i, ] contains posterior samples of the conditional treatment effect for patient i
  CTE_spls = MCMC_spls_trt$yhat.test - MCMC_spls_ctl$yhat.test
  CTE_spls = t(CTE_spls)
  # CTE_mean: length N_trial vector
  #           CTE_mean[i] is posterior mean of the conditional treatment effect for patient i
  CTE_mean = rowMeans(CTE_spls)

  # Calculating conditional average treatment effect (CATE)
  # CATE_spls: length niter vector containing posterior samples of CATE
  CATE_spls = colMeans(CTE_spls)
  
  # Calculating population average treatment effect (PATE)
  # zeta_spls: N_trial * niter matrix (after transpose)
  #            CTE_spls[i, ] contains posterior samples of the weight of patient i in 
  #            calculating CATE in Bayesian bootstrap
  zeta_spls = rdirichlet(niter, rep(0.01 + 1, N_trial))
  zeta_spls = t(zeta_spls)
  # PATE_spls: length niter vector containing posterior samples of PATE
  PATE_spls = colSums(CTE_spls * zeta_spls)
  
  # OUTPUT 
  return(list(PATE_spls = PATE_spls, 
              CATE_spls = CATE_spls,
              CTE_mean = CTE_mean))

}









