load("output/Logistic/rep_logit_scale_n.RData")
load("output/gold_standards/gs_rep_logit_scale_n.RData")

exponents <- exponents[-8]
stat_means <- matrix(0, nrow = length(exponents), ncol = 6)
overall_means <- matrix(0, nrow = length(exponents), ncol = 6)

for(method in c(1,2,5,6)){
  for( j in 1:length(exponents)){
    overall_eff <- stat_eff <- rep(0,2)
    for(i in 1:2){
      stat_eff[i] <- mean((coef_all[method,i,j,] - coef_ref[j,i])^2*Nit_all[method, j,])
      overall_eff[i] <- mean((coef_all[method,i,j,] - coef_ref[j,i])^2*times_all[method, j,])
    }
    stat_means[j,method] <- median(stat_eff)
    overall_means[j,method] <- median(overall_eff)
  }
}
par(mfrow = c(1,2), mar=c(4,4,3,2))
matplot(x = exponents,log2(stat_means[,2]/stat_means[,c(1,2,5,6)]), 
        xlab = '', ylab ='',ylim = c(-16,2.5), 
        type = 'o', pch=1, lty = 1,
        col=c('black','green','blue','magenta'))

legend('bottomleft', legend = c("ZZ","Gibbs","SS","CV"), lwd=1, col=1:4)
title(xlab="log2(number of observations)", ylab="log2(RSE)", main = 'Relative statistical efficiency')

matplot(x = exponents,log2(overall_means[,2]/overall_means[,c(1,2,5,6)]), 
        xlab = '', ylab ='',ylim = c(-5,8),
        type = 'o', pch=1, lty = 1,
        col=c('black','green','blue','magenta'))

title(xlab="log2(number of observations)", ylab="log2(RE)", main = 'Relative efficiency')
legend('topleft',c("ZZ","Gibbs","SS","CV"),
       col = c('black','green','blue','magenta'), lwd = 1)

