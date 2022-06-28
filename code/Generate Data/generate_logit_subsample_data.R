
generate.logistic.data <- function(beta, n.obs, Sig = diag(1,length(beta)), alpha = 1) {
  p <- length(beta)
  dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
  dataX <- dataX * (matrix(runif(n.obs*p), n.obs, p) < alpha)
  vals <- dataX %*% as.vector(beta)
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}
n.dimension <- 15
s <- 2
beta <- c(rep(1, s), rep(0, n.dimension-s))

set.seed(1)
logisticDataf <- generate.logistic.data(beta, 2^14)
save(logisticDataf, file = "Data/data_subsample_Logit.RData")