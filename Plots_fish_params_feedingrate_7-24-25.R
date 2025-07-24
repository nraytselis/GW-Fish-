library(dplyr)
library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)
library(cowplot)
library(wesanderson)

joint_fitB <- readRDS(file = "Joint_GW_dp50kfB.RDA")
joint_fitC <- readRDS(file = "Joint_GW_dp50kfC.RDA")

fish_fit = c(joint_fitB, joint_fitC)

parameter_chains = data.frame(fish_fit[[16]]$samples)

parameter_chains_feeding = select(parameter_chains, c('f','f_J','f_N'))

parameter_chains_feeding$f_Juv = parameter_chains_feeding$f * parameter_chains_feeding$f_J
parameter_chains_feeding$f_Naup = parameter_chains_feeding$f * parameter_chains_feeding$f_N

parameter_chains_feeding = select(parameter_chains_feeding, c('f','f_Juv','f_Naup'))

longdf = parameter_chains_feeding %>% pivot_longer(cols=1:3,names_to = "Stage",values_to = "Value")

longdfsummary = longdf %>% group_by(Stage) %>% summarise(mean = mean(Value), sd_F = sd(Value),
                                                                           n_F = n(),
                                                                           se_F = sd_F / sqrt(n_F))
ggplot(longdfsummary, aes(x = Stage, y = mean)) +
  geom_point() +
  theme_minimal() +
  geom_errorbar(aes(ymin=mean-se_F, ymax=mean+se_F, alpha =0.2)) + theme_classic() +theme(legend.position="none") +
  theme(text = element_text(size = 15)) + labs(x = "Copepod Stage",y="Feeding Rate")
