## Read argument
library(utils, methods)
library(rstan)
library(rjpdmp, lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.6/")
load("Data/dataRR.Rdata")
dataX <- data$dataX; dataY <- data$dataY
dataY_hold <- data$dataY2;  dataX_hold <- data$dataX2
p0 <- 4; p <- ncol(dataX)
beta1 <- c(rep(2,p0), rep(0,p-p0))
Mode_Mod <- paste(c(rep(1,p0), rep(0, p-p0)), collapse = '')
ppi_val <- p0/p

nDatasets <- 50
nSamples <- 2^c(4:11)

mean_estimates_ZZ <- array(0, dim = c(p, length(nSamples), nDatasets))
mean_estimates_BPS_N <- array(0, dim = c(p, length(nSamples), nDatasets))
mean_estimates_BPS_S <- array(0, dim = c(p, length(nSamples), nDatasets))
mean_estimates_STAN <- array(0, dim = c(p, length(nSamples), nDatasets))

median_estimates_ZZ <- array(0, dim = c(p, length(nSamples), nDatasets))
median_estimates_BPS_N <- array(0, dim = c(p, length(nSamples), nDatasets))
median_estimates_BPS_S <- array(0, dim = c(p, length(nSamples), nDatasets))
median_estimates_STAN <- array(0, dim = c(p, length(nSamples), nDatasets))

pred_ZZ <- array(0, dim = c(length(nSamples), nDatasets))
pred_BPS_N <- array(0, dim = c(length(nSamples), nDatasets))
pred_BPS_S <- array(0, dim = c(length(nSamples), nDatasets))
pred_STAN <- array(0, dim = c(length(nSamples), nDatasets))

burn_time <- comp_time <- matrix(0,length(nSamples), nDatasets)
Niters<- Mode_found <- array(0, dim = c(4,length(nSamples), nDatasets))

## Test the number of iterations to use for burnin...
stan_data <- list(n=nrow(dataX), p=p, y=dataY, x=dataX,
                  nu_local = 1,
                  tau = (p0/(p - p0))/sqrt(nrow(dataX)),
                  c = 10)

mean_pred <- function(y, X, beta){
  MSE <- rep(0,nrow(beta))
  for( i in 1:nrow(beta)){
    MSE[i] <- mean((y - X%*%beta[i,])^2)
  }
  return(mean(MSE))
}

ns <- ds <- 1
set.seed(ds)
theta0 <- rep(0,p)
x0 <- rep(0,p)

res_trunc <- function(res, percent_samples){
  nEvents <- length(res$times)
  rnE <- floor(nEvents*percent_samples)
  res$times <- res$times[1:rnE]
  print(rnE)
  res$positions <- res$positions[,1:rnE]
  res$theta <- res$theta[,1:rnE]
  return(res)
}

ref_N <- .6
ref_S <- .6
rj_burn <- .6 
rj_val <- 0.6

fname <- paste0('output/RobustRegression/rep_rr.Rdata')

smod <- stan_model(file = 'code/Generate Output/RR_HorseShoe.stan')
seed_cnt <- 0
for( ds in 1:nDatasets ){
  ns_0 <- 1
  
  for( ns in ns_0:length(nSamples)){
    seed_cnt = seed_cnt + 1
    set.seed(seed_cnt)
    
    fit <- sampling(smod, data = stan_data, iter = 1000 + nSamples[ns], warmup = 1000, chains = 1, cores = 1)
    print(" time: ")
    print(get_elapsed_time(fit))
    sfit <- extract(fit, 'beta')$beta
    
    time_burn <- get_elapsed_time(fit)[1]
    time_to_run <- get_elapsed_time(fit)[2]
    comp_time[ns,ds] <- time_to_run
    burn_time[ns,ds] <- time_burn
    
    ## Burn session
    print(" bps_n burn: ")
    set.seed(seed_cnt)
    resn_burn <- bps_n_rr(maxTime = time_burn, dataX = dataX, datay = dataY,
                          prior_sigma2 = 10^2, x0 = x0, theta0 = theta0, ref = ref_N,
                          rj_val = rj_burn, ppi = ppi_val, nmax = 1e6)
    
    x0_N <- resn_burn$position[,length(resn_burn$times)]
    theta0_N <- resn_burn$theta[,length(resn_burn$times)]
    rm(resn_burn)
    
    print(" bps_s burn: ")
    set.seed(seed_cnt)
    ress_burn <- bps_s_rr(maxTime = time_burn, dataX = dataX, datay = dataY,
                          prior_sigma2 = 10^2, x0 = x0, theta0 = theta0, ref = ref_S,
                          rj_val = rj_burn, ppi = ppi_val, nmax = 1e6)
    x0_S <- ress_burn$position[,length(ress_burn$times)]
    theta0_S <- ress_burn$theta[,length(ress_burn$times)]
    rm(ress_burn)
    
    print(" zz burn: ")
    set.seed(seed_cnt)
    resz_burn <- zigzag_rr(maxTime = time_burn, dataX = dataX, datay = dataY,
                           prior_sigma2 = 10^2, x0 = x0, theta0 = theta0,
                           rj_val = rj_burn, ppi = ppi_val, nmax = 1e6)
    x0_Z <- resz_burn$position[,length(resz_burn$times)]
    theta0_Z <- resz_burn$theta[,length(resz_burn$times)]
    rm(resz_burn)
    
    ## Run Methods
    set.seed(seed_cnt)
    system.time(resN <- bps_n_rr(maxTime = time_to_run, dataX = dataX, datay = dataY,
                                 prior_sigma2 = 10^2, x0 = x0_N, theta0 = theta0_N, ref = ref_N,
                                 rj_val = rj_val, ppi = ppi_val, nmax = 1e6, burn = -1))
    set.seed(seed_cnt)
    system.time(resS <- bps_s_rr(maxTime = time_to_run, dataX = dataX, datay = dataY,
                                 prior_sigma2 = 10^2, x0 = x0_S, theta0 = theta0_S, ref = ref_S,
                                 rj_val = rj_val, ppi = ppi_val, nmax = 1e6, burn = -1))
    set.seed(seed_cnt)
    system.time(resZ <- zigzag_rr(maxTime = time_to_run, dataX = dataX, datay = dataY,
                                  prior_sigma2 = 10^2, x0 = x0_Z, theta0 = theta0_Z,
                                  rj_val = rj_val, ppi = ppi_val, nmax = 1e6, burn = -1))
    
    
    Mode_found[1,ns,ds] <- mean(apply(sfit, 1, function(s) all(as.numeric(abs(s) > 0.1) == c(rep(1,p0),rep(0,p-p0)) )))
    
    mean_estimates_STAN[,ns,ds] <- colMeans(sfit)
    median_estimates_STAN[,ns,ds] <- apply(sfit,2,median)
    pred_STAN[ns,ds] <- mean_pred(dataY_hold, dataX_hold, sfit)
    #rm(sfit)
    #gc()
    
    ### BPS_N
    Niters[2,ns,ds] <- length(resN$times)
    post_mean <- marginal_mean(resN$times, resN$positions, resN$theta, marginals = 1:p)
    Mode_found[2,ns,ds] <- model_probabilities(resN$times,resN$theta,models = c(1,1,1,1,rep(0,p-4)))$prob_mod
    sres <- gen_sample(resN$positions,resN$times, theta = resN$theta, nsample = nSamples[ns])
    rm(resN)
    gc()
    mean_estimates_BPS_N[,ns,ds] <- post_mean
    median_estimates_BPS_N[,ns,ds] <- apply(sres$x,1,median)
    pred_BPS_N[ns,ds] <- mean_pred(dataY_hold, dataX_hold, t(sres$x))
    
    ### BPS_S
    Niters[3,ns,ds] <- length(resS$times)
    sres <- gen_sample(resS$positions,resS$times, theta = resS$theta, nsample = nSamples[ns])
    post_mean <- marginal_mean(resS$times, resS$positions, resS$theta, marginals = 1:p)
    Mode_found[3,ns,ds] <- model_probabilities(resS$times,resS$theta,models = c(1,1,1,1,rep(0,p-4)))$prob_mod
    rm(resS)
    gc()
    
    mean_estimates_BPS_S[,ns,ds] <- post_mean
    median_estimates_BPS_S[,ns,ds] <- apply(sres$x,1,median)
    pred_BPS_S[ns,ds] <- mean_pred(dataY_hold, dataX_hold, t(sres$x))
    
    ### Zig-Zag
    sres <- gen_sample(resZ$positions,resZ$times, theta = resZ$theta, nsample = nSamples[ns])
    post_mean <- marginal_mean(resZ$times, resZ$positions, resZ$theta, marginals = 1:p)
    Mode_found[4,ns,ds] <- model_probabilities(resZ$times,resZ$theta,models = c(1,1,1,1,rep(0,p-4)))$prob_mod
    Niters[4,ns,ds] <- length(resZ$times)
    rm(resZ)
    gc()
    
    mean_estimates_ZZ[,ns,ds] <- post_mean
    median_estimates_ZZ[,ns,ds] <- apply(sres$x,1,median)
    pred_ZZ[ns,ds] <- mean_pred(dataY_hold, dataX_hold, t(sres$x))
    rm(sres)
    gc()
    
    save(file = fname,
         mean_estimates_ZZ,
         mean_estimates_BPS_N,
         mean_estimates_BPS_S,
         mean_estimates_STAN,
         median_estimates_ZZ,
         median_estimates_BPS_N,
         median_estimates_BPS_S,
         median_estimates_STAN,
         pred_ZZ,
         pred_BPS_N,
         pred_BPS_S,
         pred_STAN,
         comp_time,
         burn_time,
         Mode_found,
         Niters)
  }
}
