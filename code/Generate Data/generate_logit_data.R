## Aux Functions
generate.logistic.data <- function(beta, n.obs, Sig) {
  p <- length(beta)
  dataX <- MASS::mvrnorm(n=n.obs,mu=rep(0,p),Sigma=Sig)
  vals <- dataX %*% as.vector(beta)
  generateY <- function(p) { rbinom(1, 1, p)}
  dataY <- sapply(1/(1 + exp(-vals)), generateY)
  return(list(dataX = dataX, dataY = dataY))
}
### Simulation
table_of_sim <- as.matrix(expand.grid(c(100,200,400,800), c(100, 200, 400), c(5)))
colnames(table_of_sim) <-c('n','p','sim')

for( N_index in 1:nrow(table_of_sim)){
  (current_sim <- table_of_sim[N_index,])
  n <- as.numeric(current_sim[1])
  p <- as.numeric(current_sim[2])
  sim <- as.numeric(current_sim[3])
  
  if(sim == 1){
    ## Strongly correlated 2 variables
    s <- 1
    beta_true <- c(rep(1,s),rep(0, p-s))
    Sig <- diag(1,p,p)
    Sig[1,2] <- Sig[2,1] <- 0.9
  }
  if(sim == 2){
    ## Structured correlation
    s <- 6
    beta_true <- c(c(3,3,-2,3,3,-3),rep(0, p-s))
    Sig <- sapply(1:p, function(i) sapply(1:p, function(j) exp(-abs(i-j))))
  }
  if(sim == 3){
    ## No correlation
    s <- 6
    beta_true <- c(c(3,3,-2,3,3,-3),rep(0, p-s))
    Sig <- diag(1,p,p)
  }
  if(sim == 4){
    ## Block correlation
    rho <- 0.9
    Sig <- diag(1,p,p)
    for (i in 1:6){Sig[i,i+floor(p/2)]<-rho; Sig[i+floor(p/2),i]<-rho}
    require(MASS)
    s <- 6
    beta_true <- c(c(3,3,-2,3,3,-3),rep(0, p-s))
  }
  if(sim == 5){
    ## Block correlation
    rho <- 0.99
    Sig <- diag(1,p,p)
    for (i in 1:6){Sig[i,i+floor(p/2)]<-rho; Sig[i+floor(p/2),i]<-rho}
    require(MASS)
    s <- 6
    beta_true <- c(c(3,3,-2,3,3,-3),rep(0, p-s))
  }
  
  set.seed(1)
  data <- generate.logistic.data(beta_true, n, Sig)
  name <- paste0("Data/data_Logit_n%d_p%d_sim_%d.RData")
  name <- sprintf(name,n,p,sim)
  save(data,beta_true, file = name, version = 2)
}
