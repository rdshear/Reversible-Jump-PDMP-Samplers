library('rjpdmp', lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")
n.dimension <- 15

con_mean <- function(times, positions, thetas, theta_c, burnin = 1){
  eps = 1e-12
  nonZeroindicies <- which(abs(theta_c) > 1e-10)
  nZ <- length(nonZeroindicies)
  maxIter <- length(times)
  cond_mean <- rep(0, nZ)
  for( Zi in 1:nZ){
    zi <- nonZeroindicies[Zi]
    beta_mean <- 0
    total_time <- 0
    for(i in (burnin):(maxIter-1)){
      tauv <- (times[i+1] - times[i])
      if(all(abs(theta_c-abs(thetas[,i])) <1e-10)){
        total_time <- total_time + tauv
        beta_mean = beta_mean + (tauv*positions[zi,i] + thetas[zi,i]*tauv^2/2)
        cond_mean[Zi] = beta_mean/total_time
      }
    }
  }
  return(cond_mean)
}

n.experiments <- 100
exponents <- 7:14
n.epochs <- 1e6
maxT <- 15

datname <- "Data/data_subsample_Logit.RData"
name <- "MCMC/rep_logit_scale_n.RData"
i<-j<-1
nmethods <- 6
times_all <- array(NA, dim = c(nmethods, length(exponents), n.experiments))
coef_all <- array(NA, dim = c(nmethods, 2, length(exponents), n.experiments))
Nit_all <- array(NA, dim = c(nmethods, length(exponents), n.experiments))
coef_ref <- matrix(0, nrow = length(exponents), 2)

prior_on <- 1/n.dimension
rj_val <- 0.6
prior_sigma2 <- 10
s <- 2
beta <- c(rep(1, s), rep(0, n.dimension-s))
gamma <- as.numeric(abs(beta) >= 1e-10)
ind_anal <- 1

load(datname)

library(Rcpp, lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")
cppFunction('double log_logit(arma::mat X, arma::vec y, arma::vec x){
  arma::vec z = X*x;
  return(  sum( log(1+exp(z)) - y%z)) ;
}', depends = c("RcppArmadillo"))

## Return the negative log posterior for logistic regression -- used in finding the CV
log_logitistic <- function(X, y, sigma2, x){
  return(log_logit(X, y, x) - sum(dnorm(x,sd = sqrt(sigma2), log=T)))
}

methods <- c("ZZ","Gibbs","","","SS","CV")
rownames(coef_all) <- methods
print(paste("max T:",maxT))

for (j in 1:length(exponents)) {
  n.observations <- 2^exponents[j]
  cat(n.observations,'observations.\n')
  coef_zigzag <- matrix(nrow = n.experiments, ncol = 1)
  coef_gibbs <- matrix(nrow = n.experiments, ncol = 1)
  coef_ss <- matrix(nrow = n.experiments, ncol = 1)
  coef_cv <- matrix(nrow = n.experiments, ncol = 1)
  ess_zigzag <- matrix(nrow = n.experiments, ncol = 1)
  ess_gibbs <- matrix(nrow = n.experiments, ncol = 1)
  ess_ss <- matrix(nrow = n.experiments, ncol = 1)
  ess_cv <- matrix(nrow = n.experiments, ncol = 1)
  essps_zigzag <- matrix(nrow = n.experiments, ncol = 1)
  essps_gibbs <- matrix(nrow = n.experiments, ncol = 1)
  essps_ss <- matrix(nrow = n.experiments, ncol = 1)
  essps_cv <- matrix(nrow = n.experiments, ncol = 1)
  
  # new data set
  logisticData <- list(dataX = logisticDataf$dataX[1:n.observations,],
                       dataY = logisticDataf$dataY[1:n.observations])
  
  cvref <- rep(0, n.dimension)
  cvref[which(gamma>0)] <- optim(par = beta[which(gamma>0)], fn = log_logitistic,
                                 X = logisticData$dataX[,which(gamma>0), drop=F],
                                 y = logisticData$dataY, sigma2=prior_sigma2, method = "BFGS")$par
  print(cvref)
  beta0 <- cvref
  
  print(j)
  for (i in 1:n.experiments) {
    
    # canonical zig zag
    cat(1)
    ptm.start <- proc.time()[3]
    set.seed(i);result.ZZ <- zigzag_logit(nmax = n.epochs, 
                                          dataX = logisticData$dataX,
                                          datay = logisticData$dataY,
                                          prior_sigma2 = prior_sigma2,
                                          theta0 = gamma,
                                          x0 = beta0, 
                                          rj_val = rj_val,
                                          ppi = prior_on, 
                                          maxTime = maxT)
    ptm.stop <- proc.time()[3]
    times_all[1,j,i] <- (ptm.stop - ptm.start)
    Nit_all[1,j,i] <- length(result.ZZ$times)
    
    coef_all[1,,j,i] <- con_mean(result.ZZ$times, 
                                 result.ZZ$positions,
                                 result.ZZ$theta,
                                 theta_c = gamma)
    rm(result.ZZ)
    gc()
    
    # # Compare gibbs
    ptm.start <- proc.time()[3]
    set.seed(i);gibbs_fit <- gibbs_logit(logisticData$dataX,
                                         logisticData$dataY,
                                         beta = beta0, 
                                         gamma = gamma,
                                         ppi = prior_on, 
                                         prior_sigma2 = prior_sigma2,
                                         nsamples = n.epochs, maxTime = maxT)
    ptm.stop <- proc.time()[3]
    times_all[2,j,i] <- (ptm.stop - ptm.start)
    Nit_all[2,j,i] <- length(gibbs_fit$beta[1,])
    True_modg <- apply(gibbs_fit$gamma, 2, function(a) all(abs(a - gamma) < 1e-10))
    
    coef_all[2,,j,i] <- rowMeans(gibbs_fit$beta[1:2,which(True_modg)])
    
    cat(5)
    # Subsample basic
    ptm.start <- proc.time()[3]
    set.seed(i);result.SS <- zigzag_logit_ss(nmax = n.epochs, 
                                             dataX = logisticData$dataX, 
                                             datay = logisticData$dataY,
                                             cvref = 0*cvref,
                                             prior_sigma2 = prior_sigma2,
                                             theta0 = gamma,
                                             x0 = beta0, rj_val = rj_val,
                                             ppi = prior_on, 
                                             maxTime = maxT)
    ptm.stop <- proc.time()[3]
    times_all[5,j,i] <- (ptm.stop - ptm.start)
    Nit_all[5,j,i] <- length(result.SS$times)
    coef_all[5,,j,i] <- con_mean(result.SS$times, 
                                 result.SS$positions,
                                 result.SS$theta, 
                                 theta_c = gamma)
    rm(result.SS)
    gc()
    
    cat(6)
    # Subsample CV
    ptm.start <- proc.time()[3]
    set.seed(i);result.CV <- zigzag_logit_ss(nmax = n.epochs, 
                                             dataX = logisticData$dataX, 
                                             datay = logisticData$dataY,
                                             cvref = cvref,
                                             prior_sigma2 = prior_sigma2,
                                             theta0 = gamma,
                                             x0 = beta0, rj_val = rj_val,
                                             ppi = prior_on, 
                                             maxTime = maxT)
    ptm.stop <- proc.time()[3]
    times_all[6,j,i] <- (ptm.stop - ptm.start)
    Nit_all[6,j,i] <- length(result.CV$times)
    coef_all[6,,j,i] <- con_mean(result.CV$times, 
                                 result.CV$positions,
                                 result.CV$theta, 
                                 theta_c = gamma)
    coef_cv[i] <- coef_all[6,ind_anal,j,i]
    
    rm(result.CV)
    gc()
    save(coef_all, Nit_all, n.dimension, beta, times_all, exponents, coef_ref, file = name)
  }
}
