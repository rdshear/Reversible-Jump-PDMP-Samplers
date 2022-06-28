#### Run for table:
get_coefthing <- function(coef_all, times = NULL, ref_m = NULL, ref_p = NULL){
  if(is.null(times)){
    times <- rep(1, dim(results)[1])
  }
  p <- length(coef_all[,1,1])
  res <- matrix(0,p,2)
  
  for( i in 1:p ){
    # ppi
    rs <- if(is.null(ref_p)) mean(coef_all[i,1,]) else ref_p[i]
    res[i,1] <- mean((coef_all[i,1,] - rs)^2*times)
    # coef
    rs <- if(is.null(ref_m)) mean(coef_all[i,2,]) else ref_m[i]
    res[i,2] <- mean((coef_all[i,2,] - rs)^2*times)
  }
  return(res)
}
#calc <- 1 # 1 = ppi, 2 = mean
n_p_vals <- expand.grid(c(100,200,400,800), c(100, 200,400))
methods <- c("Gibbs", "RJ", "ZZ", "HMC", "BPS_N", "HMC_2")

for(sim in 1:5){
  tab <- matrix(0, ncol = 12, nrow(n_p_vals))
  nameCol <- rep(methods, each = 2)
  colnames(tab) <- nameCol
  rownames(tab) <- apply(n_p_vals,1,function(x) paste(x, collapse = ', '))
  for(np in 1:nrow(n_p_vals)){
    n <- n_p_vals[np,1]
    p <- n_p_vals[np,2]
    name <- paste0("output/Logistic/replication_logit_n",n,"_p",p,"_sim_",sim,".RData")
    load(name)
    inds <- which(abs(results[1,1,1,])>0)
    print(length(inds))
    print(paste("n:",n,"p:",p))
    namer <- paste0("output/gold_standards/logit_ref_n",n,"_p",p,"_sim_",sim,".RData")
    load(namer)
    for( m in 1:6){
      if(length(inds) < 1){
        print(sprintf("num on n= %d, p=%d, sim =%d, with inds = %d",n,p,sim,length(inds)))
      } else {
        # print(sprintf("num on n= %d, p=%d, sim =%d, with inds = %d",n,p,sim,length(inds)))
        coef_mse <- get_coefthing(results[,,m,inds], times[m,inds],
                                  ref_m = ref_m, ref_p = ref_p)
        tab[np,which(methods[m] == colnames(tab))[1]] <- median(coef_mse[,1])
        tab[np,which(methods[m] == colnames(tab))[2]] <- median(coef_mse[,2])
      }
    }
  }
  nameCol[2*c(1:6)-1] <- paste('ppi', methods)
  nameCol[2*c(1:6)] <- paste('mean', methods)
  colnames(tab) <- nameCol
  stab <- tab

  notRJmean <- c(6,10,2,8,12)
  notRJppi <- c(6,10,2,8,12)-1

  stab[,notRJppi] <- round(apply(tab[,notRJppi], 2, function(s) tab[,3]/s),2)
  stab[,notRJmean] <- round(apply(tab[,notRJmean], 2, function(s) tab[,4]/s),2)
  print(sprintf("Table on sim =%d",sim))
  print(knitr::kable(round(stab[,c(5,6,9,10,1,2,7,8)],2)))
  print(xtable::xtable(round(stab[,c(5,6,9,10,1,2,7,8)],2)))
}

