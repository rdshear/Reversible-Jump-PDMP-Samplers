library('rjpdmp', lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")
n.dimension <- 15
exponents <- 7:14
n.epochs <- 1e6
maxT <- 60*60*6

datname <- "Data/data_subsample_Logit.RData"
name <- "MCMC/gs_rep_logit_scale_n.RData"

i<-j<-1
coef_ref <- matrix(0, nrow = length(exponents), 2)
times <- Nit <- rep(0, length(exponents))

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


for (j in 1:length(exponents)) {
  n.observations <- 2^exponents[j]
  cat(n.observations,'observations.\n')

  # new data set
  logisticData <- list(dataX = logisticDataf$dataX[1:n.observations,],
                       dataY = logisticDataf$dataY[1:n.observations])
  
  cvref <- rep(0, n.dimension)
  cvref[which(gamma>0)] <- optim(par = beta[which(gamma>0)], fn = log_logitistic,
                                 X = logisticData$dataX[,which(gamma>0), drop=F],
                                 y = logisticData$dataY, sigma2=prior_sigma2, method = "BFGS")$par
  print(cvref)
  beta0 <- cvref
  
  
  ptm.start <- proc.time()[3]
  set.seed(i);gibbs_fit <- gibbs_logit(logisticData$dataX,
                                       logisticData$dataY,
                                       beta = beta0, 
                                       gamma = gamma,
                                       ppi = prior_on, 
                                       prior_sigma2 = prior_sigma2,
                                       nsamples = n.epochs, maxTime = maxT)
  ptm.stop <- proc.time()[3]
  times[j] <- (ptm.stop - ptm.start)
  Nit[j] <- length(gibbs_fit$beta[1,])
  True_modg <- apply(gibbs_fit$gamma, 2, function(a) all(abs(a - gamma) < 1e-10))
  
  coef_ref[j,] <- rowMeans(gibbs_fit$beta[1:2,which(True_modg)])
  
  save(coef_ref, Nit, times, file = name)
}

