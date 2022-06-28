## Read argument
args<-commandArgs(TRUE)
print("begin")
print(args)
print("end")
N_index <- eval( parse(text=args[1]) )
library(utils, methods)
library(Rcpp, lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")

## Aux Functions
tune_it <- function(n_eval=10, args, args_tune, opt_acc, method = 'hmc_logit'){
  acc_vals <- c()
  for(i in 1:n_eval){
    method_eval <- do.call(method, args)
    acc <- mean(method_eval$acc_probs)
    acc_vals <- c(acc_vals, acc)
    args[[args_tune]] = args[[args_tune]]*(1 + (acc - opt_acc))
    args[["x0"]] = method_eval$samples[,length(method_eval$acc_probs)]
    args[["x0"]] = args[["x0"]] + rnorm(n = length(args[["x0"]]), sd = 0.2)
  }
  return(list(acc = acc_vals, args = args))
}
calc_model_ppi <- function(thetas){
  df <- data.table::data.table(t(abs(thetas)) > 1e-10)
  return(as.matrix(df[,.(COUNT = .N), by = names(df)]))
}
generate.logistic.data <- function(beta, n.obs, Sig) {
  p <- length(beta)
  dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
  vals <- dataX %*% as.vector(beta)
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}

library('rjpdmp', lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")
library('nimble', lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")

### Simulation
print(sprintf("N_index is %d.\n",N_index))
table_of_sim <- as.matrix(expand.grid(c(100,200,400,800), c(100, 200, 400), c(1,2,3,4,5)))
colnames(table_of_sim) <-c('n','p','sim')
(current_sim <- table_of_sim[N_index,])
n <- as.numeric(current_sim[1])
p <- as.numeric(current_sim[2])
sim <- as.numeric(current_sim[3])
from_ind <- 1

## Prior information
rj_val <- 0.6
prior_on <- 10/p
prior_sigma2 <- 10

## Compute for use in zigzag or other.
# Simulation settings
maxT <- 120
maxI <- 10^6
nSim <- 100

name <- paste0("MCMC/replication_logit_n%d_p%d_sim_%d.RData")
name <- sprintf(name,n,p,sim)
datname <- paste0("Data/data_Logit_n%d_p%d_sim_%d.RData")
datname <- sprintf(datname,n,p,sim)
load(datname)

## Nimble setup
methods <- c('Gibbs', 'RJ', 'ZZ', 'HMC', 'BPS_N', 'HMC_2')
results <- array(0, dim = c(p, 2, length(methods), nSim))
nIters <- times <- matrix(0, nrow = length(methods), ncol = nSim)
nmodels <- matrix(0, nrow = length(methods), ncol = nSim)
code <- nimbleCode({
  for(j in 1:p) {
    theta[j] ~ dnorm(0,sd = sqrt(prior_sigma2))               # Conditional Regression coefficient
  }
  for(i in 1:n) {
    mu[i] <- inprod(theta[], X[i,])
    lambda[i] <- 1/(1+exp(-mu[i]))
    y[i] ~ dbin(lambda[i],1)          # Likelihood
  }
})
dimensions = list(theta = p,lambda = n,X = c(n,p),mu = n)
consts <- list(n = n, p=p,prior_sigma2= prior_sigma2)
dat <- list(X = data$dataX, y = data$dataY)
print("NIM")
library('methods')
rNoIndicatorModel <- nimbleModel(code, constants = consts,
                                 data = dat, dimensions = dimensions)

noIndicatorModelConf <- configureMCMC(rNoIndicatorModel)
configureRJ(conf = noIndicatorModelConf,      ## model configuration
            targetNodes = "theta", ## coefficients for selection
            priorProb = prior_on,                   ## prior probability of inclusion
            control = list(adaptInterval = 200,#mean = 0, scale = 1,
                           adaptive = TRUE))
rNoIndicatorMCMC <- buildMCMC(noIndicatorModelConf)
cNoIndicatorModel <- compileNimble(rNoIndicatorModel)
cNoIndicatorMCMC <- compileNimble(rNoIndicatorMCMC, project = rNoIndicatorModel)

log_logit_prior <- function(X, y, sigma2, x){
  return(log_logit(X, y, x) - sum(dnorm(x,sd = sqrt(sigma2), log=T)))
}
inds_fit <- which(abs(beta_true) > 0)
stoch_time <- 0.1
args <- list(maxTime = 1, nmax = 1e4,x0 = beta_true[inds_fit],
             X = data$dataX[,inds_fit, drop=F],
             y = data$dataY,prior_sigma2 = prior_sigma2, epsilon = 0.1,
             stoch_time = stoch_time, burn = 1, ppi = 1)
r2 <-tune_it(n_eval = 20, args = args, args_tune = "epsilon", opt_acc = 0.6, method = "hmc_logit") #list()
print(r2$acc)
epsilon <- r2$args$epsilon #0.58
print(epsilon)
accmat <- matrix(0, nrow = 3, ncol = nSim)

###
for( sim in from_ind:nSim){
  print(sim)
  cvglmn <- rep(0,p)
  beta0 <- rep(0,p)
  gamma_glm <- rep(0,p)

  ## HMC
  ptm.start <- proc.time()[3]
  set.seed(sim);hmc_fit <- hmc_logit_rj2(nmax = maxI, X = data$dataX, y = data$dataY,thin = 1,
                                         prior_sigma2 = prior_sigma2, epsilon = epsilon, x0 = beta0,
                                         stoch_time = stoch_time, ppi = prior_on, rw_sig2 = 0.5, maxTime = maxT)
  ptm.stop <- proc.time()[3]
  n_eval <- length(hmc_fit$acc_probs)
  accmat[1,sim] <- mean(hmc_fit$acc_probs)
  nIters[4,sim] <- n_eval
  ## calc PPI
  results[,1,4,sim] <- rowMeans(abs(hmc_fit$samples[,-c(1:ceiling(0.1*n_eval))]) >1e-10 )
  ## calc mean
  results[,2,4,sim] <- rowMeans(hmc_fit$samples[,-c(1:ceiling(0.1*n_eval))])
  ## Save Timings
  times[4,sim] <- ptm.stop - ptm.start
  # Save nmodels
  nmodels[4,sim] <- nrow(calc_model_ppi(abs(hmc_fit$samples) > 1e-10))
  rm(hmc_fit)
  gc()
  save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)

  # ## HMC Alternative
  # ptm.start <- proc.time()[3]
  # set.seed(sim);hmc_fit <- hmc_logit_rj_better(nmax = maxI, X = data$dataX, y = data$dataY,thin = 1,
  #                                              prior_sigma2 = prior_sigma2, epsilon = epsilon, x0 = beta0,
  #                                              stoch_time = stoch_time, ppi = prior_on, rw_sig2 = 0.5, maxTime = maxT)
  # ptm.stop <- proc.time()[3]
  # n_eval <- length(hmc_fit$acc_probs)
  # accmat[2,sim] <- mean(hmc_fit$acc_probs)
  # nIters[6,sim] <- n_eval
  # ## calc PPI
  # results[,1,6,sim] <- rowMeans(abs(hmc_fit$samples[,-c(1:ceiling(0.1*n_eval))]) >1e-10 )
  # ## calc mean
  # results[,2,6,sim] <- rowMeans(hmc_fit$samples[,-c(1:ceiling(0.1*n_eval))])
  # ## Save Timings
  # times[6,sim] <- ptm.stop - ptm.start
  # # Save nmodels
  # nmodels[6,sim] <- nrow(calc_model_ppi(abs(hmc_fit$samples) > 1e-10))
  # rm(hmc_fit)
  # gc()
  # save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)

  ptm.start <- proc.time()[3]
  set.seed(sim);system.time(zigzag_fit <- zigzag_logit(nmax = maxI, dataX = data$dataX, datay = data$dataY,
                                                       prior_sigma2 = prior_sigma2,theta0 = gamma_glm,
                                                       x0 = beta0, rj_val = rj_val,ppi = prior_on,
                                                       maxTime = maxT))
  ptm.stop <- proc.time()[3]
  n_evalz <- length(zigzag_fit$positions[1,])
  nIters[3,sim] <- n_evalz
  times[3,sim] <- ptm.stop - ptm.start
  results[,1,3,sim] <- model_probabilities(zigzag_fit$times,zigzag_fit$positions, marginals = 1:p, burnin = ceiling(0.1*n_evalz))$marginal_prob
  results[,2,3,sim] <- marginal_mean(zigzag_fit$times,zigzag_fit$positions,zigzag_fit$theta, marginals = 1:p, burnin = ceiling(0.1*n_evalz))
  nmodels[3,sim] <- nrow(calc_model_ppi(zigzag_fit$theta))
  rm(zigzag_fit)
  gc()
  save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)


  ## BPS N
  ptm.start <- proc.time()[3]
  set.seed(sim);system.time(ref_bps <- bps_n_logit(nmax = maxI, dataX = data$dataX, datay = data$dataY,
                                                   prior_sigma2 = prior_sigma2,theta0 = gamma_glm*rnorm(p),ref = 0.1,
                                                   x0 = beta0,rj_val = rj_val, ppi = prior_on, maxTime = maxT))
  ptm.stop <- proc.time()[3]
  n_eval <- length(ref_bps$positions[1,])
  nIters[5,sim] <- n_eval
  times[5,sim] <- ptm.stop - ptm.start
  results[,1,5,sim] <- model_probabilities(ref_bps$times,ref_bps$positions, marginals = 1:p, burnin = ceiling(0.1*n_eval))$marginal_prob
  results[,2,5,sim] <- marginal_mean(ref_bps$times,ref_bps$positions,ref_bps$theta, marginals = 1:p, burnin = ceiling(0.1*n_eval))
  nmodels[5,sim] <- nrow(calc_model_ppi(ref_bps$theta))
  rm(ref_bps)
  gc()
  save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)

  ## Gibbs
  set.seed(sim);
  ptm.start <- proc.time()[3]
  set.seed(sim);system.time(gibbs_fit <- gibbs_logit(datay = data$dataY, dataX = data$dataX, beta = beta0, gamma = gamma_glm,
                                                     ppi = prior_on, prior_sigma2 = prior_sigma2,
                                                     nsamples = maxI, maxTime = maxT))
  ptm.stop <- proc.time()[3]
  n_eval <- length(gibbs_fit$gamma[1,])
  nIters[1,sim] <- nIters[2,sim] <- n_eval
  ## calc PPI
  results[,1,1,sim] <- rowMeans(gibbs_fit$gamma[,-c(1:ceiling(0.1*n_eval))])
  ## calc mean
  results[,2,1,sim] <- rowMeans(gibbs_fit$beta[,-c(1:ceiling(0.1*n_eval))]*gibbs_fit$gamma[,-c(1:ceiling(0.1*n_eval))])
  ## Save Timings
  times[1,sim] <- ptm.stop - ptm.start
  # Save nmodels
  nmodels[1,sim] <- nrow(calc_model_ppi(gibbs_fit$gamma))
  rm(gibbs_fit)
  gc()
  save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)

  ptm.start <- proc.time()[3]
  samples <- runMCMC(cNoIndicatorMCMC, n_eval, inits = list(theta = beta0))
  ptm.stop <- proc.time()[3]
  times[2,sim] <- (ptm.stop - ptm.start)
  accmat[3,sim] <- mean(coda::rejectionRate(coda::as.mcmc(samples))[which(abs(beta_true) > 0)])
  results[,1,2,sim] <- colMeans((abs(samples)>1e-10)[-c(1:ceiling(0.1*n_eval)),])
  results[,2,2,sim] <- colMeans(samples[-c(1:ceiling(0.1*n_eval)),])
  nmodels[2,sim] <- nrow(calc_model_ppi(t(samples)))
  rm(samples)
  gc()
  save(results, times, nmodels, p, n, methods, nIters, accmat, file = name)
}