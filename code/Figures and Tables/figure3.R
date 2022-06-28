generate.rr.data <- function(beta, n, Sig, noise) {
  p <- length(beta)
  dataX <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sig)
  dataY <- dataX %*% as.vector(beta)+rnorm(n, sd = sqrt(noise))
  return(list(dataX = dataX, dataY = dataY))
}

p <- 4
n<- 120
beta <- c(0.5,0.5, rep(0,p-2))
set.seed(1)
data <- generate.rr.data(beta,n,diag(1,p), 2)
dataX <- data$dataX; dataY <- data$dataY

library(rstan)
library(bayesplot)
m = 0.5
x0<- c(-1,2, rep(0.5,p-2))
library(rjpdmp)
set.seed(1)
theta0 <- rnorm(p)*(abs(x0)>1e-10)
p0<-1
ppi_val <- p0/ncol(data$dataX)

set.seed(1)
stan_data <- list(n=n, p=p, y=as.numeric(dataY), x=dataX,
                  nu_local = 1,
                  tau = (p0/(ncol(data$dataX) -p0))/sqrt(n),
                  c = 10)
fit <- stan(file = 'code/Generate Output/RR_HorseShoe.stan', data = stan_data, warmup = 10^3,
            iter = 2*10^3, chains = 1, cores = 1, thin = 1,init = list(list(beta = x0)),
            control=list(adapt_delta=0.99, max_treedepth=15))
sfit<- extract(fit, 'beta')$beta
colnames(sfit) <- paste('beta',1:p)
time_b <- rstan::get_elapsed_time(fit)[1]
time <- rstan::get_elapsed_time(fit)[2]

set.seed(1)
resn_burn <- bps_n_rr(maxTime = time_b, dataX = dataX, datay = dataY,
                      prior_sigma2 = 10^2, x0 = x0, theta0 = theta0, ref = 0.6,
                      rj_val = 0.6, ppi = ppi_val, nmax = 10^5)
resn <- bps_n_rr(maxTime = time, dataX = dataX, datay = dataY,
                 prior_sigma2 = 10^2, x0 = resn_burn$positions[,length(resn_burn$times)],
                 theta0 = resn_burn$theta[,length(resn_burn$times)], ref = 0.6,
                 rj_val = 0.6, ppi = ppi_val, nmax = 10^5)

resz_b <- zigzag_rr(maxTime = time_b, dataX = dataX, datay = dataY,
                    prior_sigma2 = 10^2, x0 = x0, theta0 = rep(1,ncol(dataX)),
                    rj_val = 0.6, ppi = ppi_val, nmax = 10^5)
resz <- zigzag_rr(maxTime = time, dataX = dataX, datay = dataY,
                  prior_sigma2 = 10^2, x0 = resz_b$positions[,length(resz_b$times)],
                  theta0 = resz_b$theta[,length(resz_b$times)],
                  rj_val = 0.6, ppi = ppi_val, nmax = 10^5)

library(latex2exp)
par(mfrow = c(2,3))
np <- 10^3
inds <- -seq(1, 10^4)
xlim <- c(-0.05,0.85)
ylim <- c(0.14,1.15)

indsZ <- seq(1, floor(length(resz$times)*0.3))
indsN <- seq(1, floor(length(resn$times)*0.3))
plot(resn$positions[1,indsN], resn$positions[2,indsN], main='RJ-BPS (N)',xlim = xlim, ylim = ylim,
     xlab=TeX('$\\theta_1$'), ylab=TeX('$\\theta_2$'), type = 'l')
sresn <- gen_sample(resn$positions[,indsN], resn$times[indsN], resn$theta[,indsN], nsample = np)
points(sresn$x[1,], sresn$x[2,], col = 'red', pch = 20)
plot(resz$positions[1,indsZ], resz$positions[2,indsZ], main='RJ-ZigZag',xlim = xlim, ylim = ylim,
     xlab=TeX('$\\theta_1$'), ylab=TeX('$\\theta_2$'), type = 'l')
sresz <- gen_sample(resz$positions[,indsZ], resz$times[indsZ], resz$theta[,indsZ], nsample = np)
points(sresz$x[1,], sresz$x[2,], col = 'red', pch = 20)
plot(sfit[1:np,1], sfit[1:np,2], main='HorseShoe',xlim = xlim, ylim = ylim,
     xlab=TeX('$\\theta_1$'), ylab=TeX('$\\theta_2$'), pch = 20, col = 'red')
##
pa1 <- 2
pa2 <- 3
xlim <- c(-0.15, 0.61)
plot(resn$positions[pa2,indsN], resn$positions[pa1,indsN], xlim = xlim, ylim = ylim,
     main='RJ-BPS (N)',xlab=TeX('$\\theta_3$'), ylab=TeX('$\\theta_2$'), type = 'l')
points(sresn$x[pa2,], sresn$x[pa1,], col = 'red', pch = 20)
plot(resz$positions[pa2,indsZ], resz$positions[pa1,indsZ], main='RJ-ZigZag', xlim = xlim, ylim = ylim,
     xlab=TeX('$\\theta_3$'), ylab=TeX('$\\theta_2$'), type = 'l')
points(sresz$x[3,], sresz$x[2,], col = 'red', pch = 20)
plot(sfit[1:np,pa2], sfit[1:np,pa1], main='HorseShoe', xlim = xlim, ylim = ylim,
     xlab=TeX('$\\theta_3$'), ylab=TeX('$\\theta_2$'), pch = 20, col = 'red')





