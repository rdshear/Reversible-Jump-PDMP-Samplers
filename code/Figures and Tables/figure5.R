library(ggplot2)
library(grid) 
library(gridExtra)
nit <- 50
indx <- 1
nSamples <- 2^c(4:11)
NS <- 1
load("output/RobustRegression/rep_rr.Rdata")

## Pred plot
NS <- 1
df <- data.frame(Method = c(rep('stan',nit),rep('zz',nit),rep('bps_n',nit), rep('bps_s',nit)),
                 Prediction = c(pred_STAN[NS,1:nit],pred_ZZ[NS,1:nit],
                                pred_BPS_N[NS,1:nit], pred_BPS_S[NS,1:nit]),
                 Iterations = nSamples[NS])
for( NS in 2:8){
  df <- rbind(df, data.frame(Method = c(rep('stan',nit),rep('zz',nit),rep('bps_n',nit), rep('bps_s',nit)),
                             Prediction = c(pred_STAN[NS,1:nit],pred_ZZ[NS,1:nit],
                                            pred_BPS_N[NS,1:nit], pred_BPS_S[NS,1:nit]),
                             Iterations = nSamples[NS]))
}
library(Rmisc)
df_summary <- summarySE(df, measurevar="Prediction", groupvars=c("Method","Iterations"))
df_summary$Iterations <- as.factor(df_summary$Iterations)
df_summary$Method <- factor(df_summary$Method, levels =  c('stan', 'bps_n', 'bps_s', 'zz'))

pred <- ggplot(df_summary, aes(x = Iterations, y = Prediction, fill = Method, group = Method)) +
  geom_line(aes(y = Prediction, color = Method))+
  geom_ribbon(aes(ymin = Prediction-se, ymax = Prediction+se), alpha = 0.1)

pred
