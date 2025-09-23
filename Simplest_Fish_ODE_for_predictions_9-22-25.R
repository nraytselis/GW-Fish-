library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)


FishODE =function(t, y, parameters) { 
  
  N=y[1]; J=y[2]; A=y[3]
  
  with(as.list(parameters),{  
    
    Preds = ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
    Pred_A = f*(Preds/VOL)/(1 + f*h*(A+f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL) 
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    dNdt = b_M*A/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A)) - (m_N_c+d_N_c)*N - cann*(VOL/15)*A*N - Pred_N*N 
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J - Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A - Pred_A*A
    
    result = c(dNdt,dJdt,dAdt) 
    
    return(list(result)) 
  } 
  )  
}  

setwd("~/Desktop/Rscripts/Data")

fullA = readRDS("Rebound_parameters_full_disp250kff_5.RDA")
fullB = readRDS("Rebound_parameters_full_disp250kff_3.RDA")
fullC = readRDS("Rebound_parameters_full_disp250kff_05.RDA")
fullD = readRDS("Rebound_parameters_full_disp250kff_01.RDA")
full <- c(fullA, fullB, fullC, fullD)

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
samps = get_best_fit(full)
pars = samps[[1]]
Initial_conditions = c(N = 8.704566e+03, J=1.200548e+04 , A=6.386437e+02)
timespan = 365 

#A) h and a (f) both are same for all stages 
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=1,h_N=1,i_P=0)
parameters = c(pars,other_parameters)
sim1 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#B) h is same for all stages, a (f) is largest for for adults (visual preds) and smallest for nauplii
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=1,h_N=1,i_P=0)
parameters = c(pars,other_parameters)
sim2 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#C) h is largest for for adults and smallest for nauplii, a is same for all stages
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim3 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#D) h is largest for for adults and a (f) is largest for for adults (visual preds) and smallest for nauplii
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim4 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)

#Control (no fish) 
other_parameters = c(Preds = 0, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.001,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim5 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)


colnames(sim1) = c("time","N1","J1","A1")
colnames(sim2) = c("time","N2","J2","A2")
colnames(sim3) = c("time","N3","J3","A3")
colnames(sim4) = c("time","N4","J4","A4")
colnames(sim5) = c("time","N5","J5","A5")
allscenarios = bind_cols(sim1,sim2,sim3,sim4,sim5)
names(allscenarios)[1] <- "Time"

#Plot scenarios A and D
Naup = allscenarios %>% select(c(1,2,14,18)) #last one is the control
Juv = allscenarios %>% select(c(1,3,15,19)) #last one is the control
Ad = allscenarios %>% select(c(1,4,16,20)) #last one is the control
NaupLong = Naup %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )
JuvLong = Juv %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )
AdLong = Ad %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )

#relabel scenarios to be more informative 
NaupLong$scenario <- gsub("N1", "equalacrossstages", NaupLong$scenario)
NaupLong$scenario <- gsub("N4", "bothlargestforadults", NaupLong$scenario)
NaupLong$scenario <- gsub("N5", "control", NaupLong$scenario)
JuvLong$scenario <- gsub("J1", "equalacrossstages", JuvLong$scenario)
JuvLong$scenario <- gsub("J4", "bothlargestforadults", JuvLong$scenario)
JuvLong$scenario <- gsub("J5", "control", JuvLong$scenario)
AdLong$scenario <- gsub("A1", "equalacrossstages", AdLong$scenario)
AdLong$scenario <- gsub("A4", "bothlargestforadults", AdLong$scenario)
AdLong$scenario <- gsub("A5", "control", AdLong$scenario)

palette = c(equalacrossstages = "#CC5500", bothlargestforadults = "gold2",control = "olivedrab")

a = ggplot(data=NaupLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + ylim(0,1500) + theme_classic() +
  labs(x="Time (days)",y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")
b = ggplot(data=JuvLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) +ylim(0,6000) + theme_classic() +
  labs(x="Time (days)",y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")
c = ggplot(data=AdLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + ylim(0,200)+ theme_classic()+
  labs(x="Time (days)",y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")

library(ggpubr)
ggarrange(a,b,c,
          nrow = 1, ncol = 3,
          labels = c("N","J","A"),
          common.legend = TRUE,
          legend = "right")

#predator gradient
preds_gradient <- seq(0, 0.3, by = 0.01)

#scenarios
scenario_params <- list(
  equalacrossstages = c(VOL = 1, f = 1, f_J = 1, f_N = 1,
        h = 0.006, h_J = 1, h_N = 1, i_P = 0),
  
  bothlargestforadults = c(VOL = 1, f = 1, f_J = 0.1, f_N = 0.01,
        h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0),
  
  control = c(VOL = 1, f = 1, f_J = 0.1, f_N = 0.01,
              h = 0.001, h_J = 0.1, h_N = 0.01, i_P = 0)
)


sim_results <- list()

for (s in names(scenario_params)) {
  for (p in preds_gradient) {
    
    # update parameters
    parameters <- c(pars,
                    Preds = p,
                    scenario_params[[s]])
    
    # run simulation
    sim <- ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = parameters,
      method = "lsoda",
      func = FishODE
    )
    
    # store results
    sim_df <- as.data.frame(sim) %>%
      mutate(Preds = p,
             Scenario = s)
    
    sim_results[[paste(s, p, sep = "_")]] <- sim_df
  }
}

# bind into one dataframe
all_sims <- bind_rows(sim_results, .id = "run_id")

library(dplyr)

#take the last time point from each simulation
endpoints <- all_sims %>%
  group_by(Scenario, Preds) %>%
  slice_tail(n = 1) %>%   # last time step
  ungroup()

PredsN = ggplot(endpoints, aes(x = Preds, y = N, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")

PredsJ = ggplot(endpoints, aes(x = Preds, y = J, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")

PredsA = ggplot(endpoints, aes(x = Preds, y = A, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Copepod Density, L' ^ -1)) + scale_color_manual(values = palette, name = "scenario")

ggarrange(a,b,c,PredsN,PredsJ,PredsA,
          nrow = 2, ncol = 3,
          labels = c("N","J","A","N","J","A"),
          common.legend = TRUE,
          legend = "right")

# Stage <- c("A","J","N","A","J","N")
# Scenario <-c("A","A","A","B","B","B")
# HValues <- c("0.006","0.006","0.006","0.006","0.0006","0.00006")
# handling <- data.frame(Stage=Stage,Scenario=Scenario,HValue=HValues) 
# handlingtime= ggplot(data=handling,aes(x=Stage,y=HValue,group=Scenario,color=Scenario)) + geom_point(size=3) + theme_classic() + labs(y="Handling Time")
# 
# Stage <- c("A","J","N","A","J","N")
# Scenario <-c("A","A","A","B","B","B")
# AValues <- c("40","40","40","40","4","0.4")
# attach <- data.frame(Stage=Stage,Scenario=Scenario,AValues=AValues) 
# attackrates = ggplot(data=handling,aes(x=Stage,y=AValues,group=Scenario,color=Scenario)) + geom_point(size=3) + theme_classic() + labs(y="Attack Rate")
# 
# ggarrange(handlingtime,attackrates,
#           nrow = 1, ncol = 2,
#           labels = c("handlingtime","attackrates"),
#           common.legend = TRUE,
#           legend = "right")
