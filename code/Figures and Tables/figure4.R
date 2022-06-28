library(ggplot2)
library(grid) 
library(gridExtra)
nit <- 50
indx <- 1
nSamples <- 2^c(4:11)
NS <- 1
load("output/RobustRegression/rep_rr.Rdata")

## Process Data
df_estimate_coord <- data.frame(Method = c(rep('stan',nit),rep('zz',nit),rep('bps_n',nit), rep('bps_s',nit)),
                 Estimate = c(mean_estimates_STAN[indx,NS,1:nit],mean_estimates_ZZ[indx,NS,1:nit],
                              mean_estimates_BPS_N[indx,NS,1:nit], mean_estimates_BPS_S[indx,NS,1:nit]),
                 Iterations = nSamples[NS])
df_mode_found <- data.frame(Method = c(rep('stan',1),rep('zz',1),rep('bps_n',1), rep('bps_s',1)),
                            MissedMode = c(sum(Mode_found[1,NS,1:nit]< 0.1),
                                           sum(Mode_found[4,NS,1:nit]< 0.1),
                                           sum(Mode_found[2,NS,1:nit]< 0.1),
                                           sum(Mode_found[3,NS,1:nit]< 0.1)),
                            Iterations = nSamples[NS])
for( NS in 2:8){
  df_estimate_coord <- rbind(df_estimate_coord, data.frame(Method = c(rep('stan',nit),rep('zz',nit),rep('bps_n',nit), rep('bps_s',nit)),
                             Estimate = c(mean_estimates_STAN[indx,NS,1:nit],mean_estimates_ZZ[indx,NS,1:nit],
                                          mean_estimates_BPS_N[indx,NS,1:nit], mean_estimates_BPS_S[indx,NS,1:nit]),
                             Iterations = nSamples[NS]))
  ## Order in the mode_found data (Stan, BPS_N, BPS_S, ZigZag)
  df_mode_found <- rbind(df_mode_found, data.frame(Method = c(rep('stan',1),rep('zz',1),rep('bps_n',1), rep('bps_s',1)),
                                                   MissedMode = c(sum(Mode_found[1,NS,1:nit]< 0.1),
                                                                  sum(Mode_found[4,NS,1:nit]< 0.1),
                                                                  sum(Mode_found[2,NS,1:nit]< 0.1),
                                                                  sum(Mode_found[3,NS,1:nit]< 0.1)),
                                                   Iterations = nSamples[NS]))
}
df_estimate_coord$Iterations <- as.factor(df_estimate_coord$Iterations)
df_estimate_coord$Method <- factor(df_estimate_coord$Method, levels =  c('stan', 'bps_n', 'bps_s', 'zz'))
df_mode_found$Iterations <- as.factor(df_mode_found$Iterations)
df_mode_found$Method <- factor(df_mode_found$Method, levels =  c('stan', 'bps_n', 'bps_s', 'zz'))
df_mode_found$`Didn't find the Mode` <- df_mode_found$MissedMode

# settings for zoom plot
zoomtheme <- theme(legend.position="none", axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),axis.title.y=element_blank(),
                   panel.grid.major = element_line(colour = "light grey"),#element_blank(),
                   panel.grid.minor = element_line(colour = "light grey"),
                   panel.background = element_rect(color='red', fill="white"),
                   plot.margin = unit(c(0,0,1,1),"mm"))

p.zoom <- ggplot(df_estimate_coord, aes(x = Iterations, y = Estimate, colour=Method)) +
  geom_boxplot() + zoomtheme
p.full <- ggplot(df_estimate_coord, aes(x = Iterations, y = Estimate,  colour = Method)) +
  geom_boxplot()+
  coord_cartesian(ylim = c(2.25,2.43))
g <- ggplotGrob(p.zoom)
box <- p.full +
  annotation_custom(grob = g, xmin = 5, xmax = Inf,
                    ymin = 2.25, ymax = 2.3)

ggbar <- ggplot(df_mode_found, aes(x = Iterations, y = `Didn't find the Mode`, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge())

boxplots_RR <- ggpubr::ggarrange(box,ggbar, nrow = 2, heights = c(2, 0.7))
boxplots_RR
