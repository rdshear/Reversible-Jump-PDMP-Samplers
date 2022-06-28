
generate.rr.data <- function(beta1, n1, n2, Sig) {
  p <- length(beta1)
  dataX <- MASS::mvrnorm(n=n1,mu=rep(0,p),Sigma=Sig)
  dataY <- dataX %*% as.vector(beta1)+ rt(n1, df = 1)
  
  dataX2 <- MASS::mvrnorm(n=n2,mu=rep(0,p),Sigma=Sig)
  dataY2 <- dataX2 %*% as.vector(beta1)+ rt(n2, df = 1)
  
  return(list(dataX = dataX, dataY = c(dataY),
              dataX2 = dataX2, dataY2 = c(dataY2)))
}
p <- 200
n1<- 100
n2 <- 100
p0 <- 4
beta1 <- c(rep(2,p0), rep(0,p-p0))
Sig <- solve(as.matrix(ar.matrix::Q.AR1(p,1,0.5)))
set.seed(1)
data <- generate.rr.data(beta1,n1,n2,Sig)
