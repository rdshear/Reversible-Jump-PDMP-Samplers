## Read argument
args<-commandArgs(TRUE)
print("begin")
print(args)
print("end")
N_index <- eval( parse(text=args[1]) )
library(utils, methods)
library(Rcpp, lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")

## Aux Functions
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

### Simulation
print(sprintf("N_index is %d.\n",N_index))
table_of_sim <- as.matrix(expand.grid(c(100,200,400,800), c(100,200,400), c(1,2,3,4,5)))
colnames(table_of_sim) <-c('n','p','sim')
(current_sim <- table_of_sim[N_index,])
n <- as.numeric(current_sim[1])
p <- as.numeric(current_sim[2])
sim <- as.numeric(current_sim[3])
from_ind <- 1

library('rjpdmp', lib.loc = "/home/suttonm/R/x86_64-pc-linux-gnu-library/3.4/")

## Prior information
rj_val <- 0.6
prior_on <- 10/p
prior_sigma2 <- 10

## Compute for use in zigzag or other.
# Simulation settings
maxT <- 60*60*48
maxI <- 5*10^6

name <- paste0("ref/logit_ref_n%d_p%d_sim_%d.RData")
name <- sprintf(name,n,p,sim)
datname <- paste0("Data/data_Logit_n%d_p%d_sim_%d.RData")
datname <- sprintf(datname,n,p,sim)
load(datname)

gibbs_b <- gibbs_logit(datay = data$dataY, dataX = data$dataX, beta = rep(0,p), gamma = rep(0,p),
                                                   ppi = prior_on, prior_sigma2 = prior_sigma2,
                                                   nsamples = 100, maxTime = maxT)

ptm.start <- proc.time()[3]
set.seed(sim);system.time(gibbs_fit <- gibbs_logit_light(datay = data$dataY, dataX = data$dataX, 
						beta = gibbs_b$beta[,100], gamma = gibbs_b$gamma[,100],
                                                   ppi = prior_on, prior_sigma2 = prior_sigma2,
                                                   nsamples = maxI, maxTime = maxT))
ptm.stop <- proc.time()[3]
n_eval <- gibbs_fit$niter
## calc PPI
ref_p <- gibbs_fit$marg_gamma
## calc mean
ref_m <- gibbs_fit$marg_beta
## Save Timings
time_gibbs <- ptm.stop - ptm.start
save(ref_p,ref_m, time_gibbs, n_eval,  file = name)
