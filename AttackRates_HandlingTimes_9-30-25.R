library(dplyr)
library(tidyverse)
library(ggpubr)
#bring in mcmc chains
setwd("~/Desktop/Rscripts/Data/Fish")
chainsA <- readRDS("Joint_GW_full_A2.RDA")
chainsB <- readRDS("Joint_GW_full_B2.RDA")
chainsC <- readRDS("Joint_GW_full_C2.RDA")

fishchains = c(chainsA,chainsB,chainsC)

#get best fit
get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump,
       max(chain.list[[which.max(chain.scores)]]$log.p))
  
}

samps = get_best_fit(fishchains)
parameters = samps[[1]]
variances = samps[[2]]

f = numeric()
f_J = numeric()
f_N = numeric()
im = numeric() #maximum ingestion rate
imJ = numeric()
imN = numeric()

for(i in 1:length(fishchains)){
  
  f = c(f,fishchains[[i]]$samples[50001:250000,"f"]) 
  f_J = c(f_J,fishchains[[i]]$samples[50001:250000,"f"]*fishchains[[i]]$samples[50001:250000,"f_J"]) 
  f_N = c(f_N,fishchains[[i]]$samples[50001:250000,"f"]*fishchains[[i]]$samples[50001:250000,"f_N"]) 
  im = c(im,1/fishchains[[i]]$samples[50001:250000,"h"]/40) 
  imJ = c(imJ,1/(fishchains[[i]]$samples[50001:250000,"h"]*fishchains[[i]]$samples[50001:250000,"h_J"]/40))
  imN = c(imN,1/(fishchains[[i]]$samples[50001:250000,"h"]*fishchains[[i]]$samples[50001:250000,"h_N"]/40))
}

attackratesdf = data.frame(f,f_J,f_N)

meanattack = colMeans(attackratesdf)

maxingestionrate = data.frame(im,imJ,imN)

meanmaxingestion = colMeans(maxingestionrate)


SEM = function(x){
  sd(x)/ sqrt(length(x)) 
} 

sdattack = apply(X=attackratesdf,MARGIN=2,FUN=SEM) 
sdmaxingestionrate = apply(X=maxingestionrate,MARGIN=2,FUN=SEM) 


#ggplot(data=maxingestionrate, aes(x=1,y=im)) + geom_boxplot() + scale_y_log10()

pal = c("#77AADD","#44BB99","#BBCC33") 

maxingestionratelong = maxingestionrate %>% pivot_longer(cols = c(1:3),names_to = "Stage",values_to = "Values" )

plotingestion = ggplot(data=maxingestionratelong, aes(x=Stage,y=Values,fill=Stage)) + geom_boxplot() + scale_y_log10() + theme_classic() + scale_fill_manual(values = pal) +
  theme(axis.text = element_text(size = 15))


attackratesdflong = attackratesdf %>% pivot_longer(cols = c(1:3),names_to = "Stage",values_to = "Values" )

plotattack = ggplot(data=attackratesdflong, aes(x=Stage,y=Values,fill=Stage)) + geom_boxplot() + scale_y_log10() + theme_classic() + scale_fill_manual(values = pal) +
  theme(axis.text = element_text(size = 15))


ggarrange(plotingestion,plotattack,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          legend = "none")



#get best fit
get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump,
       max(chain.list[[which.max(chain.scores)]]$log.p))
  
}


samps = get_best_fit(fishchains)
parameters = samps[[1]]
variances = samps[[2]]

n_samples <- nrow(chainsA[[1]]$samples)

paramdf <- data.frame(
  ID = names(parameters),
  estimate = parameters,
  sd = sqrt(diag(variances)),
  se = sqrt(diag(variances)) / sqrt(n_samples)  
)

paramdfattackrates = paramdf %>% filter(ID %in% c("f", "f_J","f_N"))

paramdfattackrates$ID = gsub("f_N", "N", paramdfattackrates$ID)
paramdfattackrates$ID = gsub("f_J", "J", paramdfattackrates$ID)
paramdfattackrates$ID = gsub("f", "A", paramdfattackrates$ID)


pal = c("#77AADD","#44BB99","#BBCC33") 

attack = ggplot(data=paramdfattackrates,aes(x=as.factor(ID),y = estimate, group = as.factor(ID), color = as.factor(ID))) + geom_point(size = 5) +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.4) + theme_classic() + scale_y_log10()+
  labs(x="Copepod Stage", y ="log(Attack Rate)") + scale_color_manual(values = pal) 


handlingtime = paramdf %>% filter(ID %in% c("h", "h_J","h_N"))

handlingtime$ID = gsub("h_N", "N", handlingtime$ID)
handlingtime$ID = gsub("h_J", "J", handlingtime$ID)
handlingtime$ID = gsub("h", "A", handlingtime$ID)

handle = ggplot(data=handlingtime,aes(x=as.factor(ID),y = estimate,group = as.factor(ID), color = as.factor(ID))) + geom_point(size = 5) +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.4) + theme_classic() + scale_y_log10()+
  labs(x="Copepod Stage", y ="log(Handling Rate)") + scale_color_manual(values = pal) 


ggarrange(attack,handle,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          legend = "none")
