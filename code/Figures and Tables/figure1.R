library(rjpdmp)
library(latex2exp)

## Sample the Continuous Spike-Slab 
zigzag_cts_ss <- function(nmax = 100, Tmax = Inf, d = 0.5, Cval = 1, tauval = 1, alpha = 0.5){
  
  Uprior <- function(b){
    expB <- exp(log(d/((1-d)*Cval))+b^2*(Cval^2-1)/(2*tauval^2*Cval^2))
    res1<- log(1 + expB/Cval^2)-log(1+expB)
    if(expB > 10^15){
      return(b/(Cval^2*tauval^2))
    } else {
      return(b/tauval^2*exp(res1))
    }
  }
  
  ## Init
  theta <- 1
  beta <- rnorm(1)
  t <- 0 
  
  ## store
  tt <- rep(0, nmax)
  bbeta <- matrix(0, nrow = 1, nmax)
  bbeta[,1] <- beta
  maxGrad <- (1/tauval^2 + (1-d)/(d*tauval^2*Cval^3))/(1+(1-d)/(Cval*d))
  
  ## simulate time:
  a <- Uprior(beta)
  b <- maxGrad*1.01
  tau <- rjpdmp::linear_inv_t(a, b, runif(1))
  
  n <- 2
  while( n  <= nmax & t < Tmax){
    beta <- beta + tau*theta
    t <- t + tau
    upper <- max(0, a+b*tau) 
    val <- max(0, theta*(Uprior(beta))) 
    if( val > upper){
      print(paste("upper ",upper," val ",val))
    }
    if( runif(1) < val/upper){
      theta = -theta
      tt[n] <- t
      bbeta[,n] <- beta
      n <- n+1
    }
    a <- theta*Uprior(beta) # small recalculate
    tau <- rjpdmp::linear_inv_t(a, b, runif(1))
  }
  print(n)
  tt <- tt[1:(n-1)]
  bbeta <- bbeta[,1:(n-1), drop = F]
  res <- list( tt = tt, xx = bbeta)
  return( res )
}

d = 0.5
Ctau <- sqrt(16);
Cval <- 5;  tauval <- Ctau/Cval
par(mfrow = c(1,3))
set.seed(2);reszz <- zigzag_cts_ss(Tmax = 110, d = d, Cval = Cval, tauval = tauval)
sres1 <- rjpdmp::gen_sample(positions = reszz$xx,times = reszz$tt, nsample = 800)
plot(sres1$t, sres1$x[1,], type = 'l', xlab = 't', ylab = TeX('$\\theta(t)$'), main = paste('Continuous spike and slab (c = 1/5)',collapse = '',sep=''), ylim = c(-12,12), xlim = c(0,100))
Cval <- 30;  tauval <- Ctau/Cval
set.seed(1);reszz2 <- zigzag_cts_ss(Tmax = 107,nmax = 10^6, d = d, Cval = Cval, tauval = tauval)
sres2 <- rjpdmp::gen_sample(positions = reszz2$xx,times = reszz2$tt, nsample = 800)
plot(sres2$t, sres2$x[1,], type = 'l', xlab = 't', ylab = TeX('$\\theta(t)$'), main = paste('Continuous spike and slab (c = 1/30)',collapse = '',sep=''), ylim = c(-12,12), xlim = c(0,100))

## Spike-slab part:
time_to_zerozz <- function(x, theta){
  if(abs(x) < 1e-10) {
    return(Inf)
  } else if( x*theta > 0){
    # if traveling away from 0
    return(Inf)
  } else {
    return( -x/theta )
  }
}
get_time <- function(x, theta, a, b, alpha, d, mu, sd){
  if( abs(theta) > 1e-10){
    hit_time <- time_to_zerozz(x,theta)
    t_norm <-  rjpdmp::linear_inv_t(a, b, runif(1))
    return(min(hit_time, t_norm))
  } else{
    u = -log(runif(1))
    frac <- dnorm(0, mean = mu, sd = sd)
    lambda = alpha*frac*d/(1-d)
    return(u/lambda)
  }
}
zigzag_ss <- function(nmax = 10^5, Tmax = Inf, sigma, x0 = 0, mu = 0, d = 0.5, alpha = 0.5){
  
  ## Init
  theta <- -1
  beta <- x0
  t <- 0 
  
  ## store
  tt <- rep(0, nmax)
  bbeta <- matrix(0, nrow = 1, nmax)
  bbeta[,1] <- beta
  
  ## simulate time:
  a <- theta*(beta - mu)/sigma^2
  b <- theta^2/sigma^2
  tau <- get_time(beta, theta, a, b, alpha, d, mu, sigma)
  
  n <- 2
  while( n  <= nmax & t < Tmax){
    beta <- beta + tau*theta
    t <- t + tau
    
    if( abs(beta) < 1e-10 ){
      beta = 0
      # if currently have this index removed
      if( abs(theta) < 1e-10 ){
        theta = (runif(1) < 0.5)*2-1 # +-1 with equal prob
        tt[n] <- t
        bbeta[,n] <- beta
        n <- n+1
      } else {
        if(runif(1) < alpha){
          theta = 0
          tt[n] <- t
          bbeta[,n] <- beta
          n <- n+1
        }
      }
    } else {
      theta = -theta
      tt[n] <- t
      bbeta[,n] <- beta
      n <- n+1
    }
    a <- theta*(beta - mu)/sigma^2
    b <- theta^2/sigma^2
    tau <- get_time(beta, theta, a, b, alpha, d, mu, sigma)
  }
  tt <- tt[1:(n-1)]
  bbeta <- bbeta[,1:(n-1), drop=F]
  res <- list( tt = tt, xx = bbeta)
  return( res )
}

set.seed(1);reszz3 <- zigzag_ss(nmax = 10^5,Tmax = 110,sigma = Ctau, d = 0.5, alpha = 0.5)
sres3 <- rjpdmp::gen_sample(positions = reszz3$xx,times = reszz3$tt, nsample = 800)
plot(sres3$t, sres3$x[1,], type = 'l', xlab = 't', ylab = TeX('$\\theta(t)$'), main = paste('Dirac spike and slab',collapse = '',sep=''), ylim = c(-12,12), xlim = c(0,100))




