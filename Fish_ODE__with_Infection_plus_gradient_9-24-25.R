# Loading packages
library(Matrix)
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(deSolve)
library(coda)

#need to fix disconnect between Preds and L3F 

Pond_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]; L3F = y[5+latent_stages] #Preds = y[5+latent_stages]; 
    VOL = 1
    
    Preds = Preds #ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
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
    
    #presence of adults in initial conditions automatically seeds the exposed class (Es) through lambda*A
    # development of all stages
    latent_progression = latent_rate*Es
    # lost to next stage   #death      #gained from last stage
    dEsdt = -latent_progression - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)]) - Pred_A*Es
    
    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I - Pred_A*I
    
    #dPredsdt = 0 #convEff*(Pred_N*N + Pred_J*J + Pred_A*A + Pred_A*sum(Es)) - d_F*Preds 
    
    dL3Fdt <- ifelse(Preds>0,Pred_A*I - d_W*L3F - d_F*L3F,0)  
    
    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt,dL3Fdt) #dPredsdt
    
    
    return(list(result))
  }
  )
}

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

parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032
parameters["d_W"] = 0.05
parameters["d_F"] = 0.05
#parameters["convEff"] = 0.001 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)
parameters["Preds"] = 0.125
#parameters["Period"] = 100000

parameters = unlist(parameters)
Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values

###
Initial_conditions = c(N = 8.704566e+03, J=1.200548e+04 , A=6.386437e+02, Exposed_values, I = 0, L3F = 0)
timespan = 365

PondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                         method="lsoda", func=Pond_ODE)) 



#reformat sim data frame
PondSim[,"Es"] = rowSums(PondSim) - PondSim[,"N"] - PondSim[,"J"]- PondSim[,"A"] - PondSim[,"I"] - PondSim[,"time"] - PondSim[,"L3F"]
PondSim[,"Exposed"] = apply(X=PondSim[,which(str_detect(colnames(PondSim),"E"))],MARGIN=1,FUN = sum) 

ggplot(PondSim, aes(x = time)) +
  geom_line(aes(y = L3F), color = "red") +
  theme_minimal()

#reformat for plotting
PondSim = PondSim %>% select(time,N,J,A,Exposed,I,L3F) #Preds
PondSim = PondSim %>% pivot_longer(cols = c(N,J,A,Exposed,I,L3F)) #Preds

#plot sim 
p1 = ggplot(PondSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + 
  geom_line() + ylab("density per L") + theme_minimal() + scale_y_log10() +xlim(0,100)
p1



####Pond sim over a fish density gradient####
timespan = 365
preds_gradient <- seq(0, 1, by = 0.01)
Initial_conditions = c(N = 8.704566e+03, J=1.200548e+04 , A=6.386437e+02, Exposed_values, I = 0, L3F = 0)

# create a list to store results
results_list <- list()

for (p in preds_gradient) {
  params <- c(parameters, Preds = p)
  
  sim <- data.frame(
    ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = params,
      hmax = 1,
      method = "lsoda",
      func = Pond_ODE
    )
  )
  
  # add Preds value as a column
  sim$Preds <- p
  
  results_list[[as.character(p)]] <- sim
}

PondSimB <- do.call(rbind, results_list)

PondSimB[,"Es"] = rowSums(PondSimB) - PondSimB[,"N"] - PondSimB[,"J"]- PondSimB[,"A"] - PondSimB[,"I"] - PondSimB[,"time"] - PondSimB[,"L3F"]
PondSimB[,"Exposed"] = apply(X=PondSimB[,which(str_detect(colnames(PondSimB),"E"))],MARGIN=1,FUN = sum) 

CleanedPondSimB = PondSimB[,c("time","Preds","N","J","A","I","L3F","Exposed")] 

ggplot(CleanedPondSimB, aes(x = time, y = L3F, group = as.factor(Preds), color = as.factor(Preds))) + geom_line()




timespan = 365*5
fish_vec = seq(from=0, to=0.075, by=0.01)
Initial_conditions = c(N = 500, J = 200, A = 25, Exposed_values, I = 0, L3F = 0)
L3F_results <- numeric(length(fish_vec))  
N_results <- numeric(length(fish_vec))  
J_results <- numeric(length(fish_vec))  
A_results <- numeric(length(fish_vec))
I_results <- numeric(length(fish_vec))


for(i in 1:length(fish_vec)){
  #introduction_times = (1:60)*30
  parameters["Preds"] = fish_vec[i]
  PondSimPreds = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, method="lsoda", func=Pond_ODE))
  
  #introduction_times[i] = introduction_time
  L3F_results[i] <- round(PondSimPreds$L3F[nrow(PondSimPreds)], 5) 
  N_results[i] <- round(PondSimPreds$N[nrow(PondSimPreds)], 5) 
  J_results[i] <- round(PondSimPreds$J[nrow(PondSimPreds)], 5) 
  A_results[i] <- round(PondSimPreds$A[nrow(PondSimPreds)], 5) 
  I_results[i] <- round(PondSimPreds$I[nrow(PondSimPreds)], 5) 
 }


plot(fish_vec, L3F_results, type = "l", xlab = "Fish Density", ylab = "Final L3F", main = "L3F vs Fish Density")

fish_vec = data.frame(fish_vec)
L3s = data.frame(L3F_results)
Nsim = data.frame(N_results)
Jsim = data.frame(J_results)
Asim = data.frame(A_results)
Isim = data.frame((I_results))

data = cbind(fish_vec,L3s,Nsim,Jsim,Asim,Isim)

data = data %>% mutate(totalcopes = N_results + J_results + A_results)


L3F = ggplot(data=data,aes(x=fish_vec,y=L3F_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Total Parasites in Fish")

Infected = ggplot(data=data,aes(x=fish_vec,y=I_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Infected Copepods")

ggarrange(L3F,Infected, 
          nrow = 1, ncol = 2,
          labels = c("A","B"),
          common.legend = TRUE,
          legend = "right")


datainfected = data[,c(1,2,6)]

datainfectedlong = datainfected %>% pivot_longer(cols=c(2,3),names_to = "InfectionRoute", values_to = "ParasiteDensity")

datainfectedlong$InfectionRoute = gsub("L3F_results", "L3F", datainfectedlong$InfectionRoute)
datainfectedlong$InfectionRoute = gsub("X.I_results.", "InfectedCope", datainfectedlong$InfectionRoute)

colors = c(
  L3F = "#77AADD",
  InfectedCope = "#EE8866") 

SelectivePred = ggplot(data=datainfectedlong, aes(x = fish_vec,y=ParasiteDensity,group=InfectionRoute,color=InfectionRoute)) + geom_line(linewidth=2) +
  theme_classic() + labs(x=expression('Fish Density, L' ^ -1),y=expression('Density, L' ^ -1)) 

#no predator specificity

timespan = 365*5
fish_vec = seq(from=0, to=0.075, by=0.01)
Initial_conditions = c(N = 500, J = 200, A = 25, Exposed_values, I = 0, L3F = 0)
L3F_results2 <- numeric(length(fish_vec))  
N_results2 <- numeric(length(fish_vec))  
J_results2 <- numeric(length(fish_vec))  
A_results2 <- numeric(length(fish_vec))
I_results2 <- numeric(length(fish_vec))

parameters = parameters
parameters2 = parameters

parameters2["f_J"] = 1
parameters2["f_N"] = 1
parameters2["h_J"] = 1
parameters2["h_N"] = 1


for(i in 1:length(fish_vec)){
  #introduction_times = (1:60)*30
  parameters2["Preds"] = fish_vec[i]
  PondSimPreds = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters2, method="lsoda", func=Pond_ODE))
  
  #introduction_times[i] = introduction_time
  L3F_results2[i] <- round(PondSimPreds$L3F[nrow(PondSimPreds)], 5) 
  N_results2[i] <- round(PondSimPreds$N[nrow(PondSimPreds)], 5) 
  J_results2[i] <- round(PondSimPreds$J[nrow(PondSimPreds)], 5) 
  A_results2[i] <- round(PondSimPreds$A[nrow(PondSimPreds)], 5) 
  I_results2[i] <- round(PondSimPreds$I[nrow(PondSimPreds)], 5) 
}


fish_vec = data.frame(fish_vec)
L3s = data.frame(L3F_results2)
Nsim = data.frame(N_results2)
Jsim = data.frame(J_results2)
Asim = data.frame(A_results2)
Isim = data.frame((I_results2))

data2 = cbind(fish_vec,L3s,Nsim,Jsim,Asim,Isim)

data2 = data2 %>% mutate(totalcopes = N_results + J_results + A_results)


L3F = ggplot(data=data2,aes(x=fish_vec,y=L3F_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Total Parasites in Fish")

Infected = ggplot(data=data2,aes(x=fish_vec,y=I_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Infected Copepods")

ggarrange(L3F,Infected, 
          nrow = 1, ncol = 2,
          labels = c("A","B"),
          common.legend = TRUE,
          legend = "right")


datainfected2 = data2[,c(1,2,6)]

datainfectedlong2 = datainfected2 %>% pivot_longer(cols=c(2,3),names_to = "InfectionRoute", values_to = "ParasiteDensity")

datainfectedlong2$InfectionRoute = gsub("L3F_results2", "L3F", datainfectedlong2$InfectionRoute)
datainfectedlong2$InfectionRoute = gsub("X.I_results2.", "InfectedCope", datainfectedlong2$InfectionRoute)

colors = c(
  L3F = "#77AADD",
  InfectedCope = "#EE8866") 

NonSelectivePred = ggplot(data=datainfectedlong2, aes(x = fish_vec,y=ParasiteDensity,group=InfectionRoute,color=InfectionRoute)) + geom_line(linewidth=2) +
  theme_classic() + labs(x=expression('Fish Density, L' ^ -1),y=expression('Density, L' ^ -1)) 



#no selectivity in feeding/attack rate
timespan = 365*5
fish_vec = seq(from=0, to=0.075, by=0.01)
Initial_conditions = c(N = 500, J = 200, A = 25, Exposed_values, I = 0, L3F = 0)
L3F_results3 <- numeric(length(fish_vec))  
N_results3 <- numeric(length(fish_vec))  
J_results3 <- numeric(length(fish_vec))  
A_results3 <- numeric(length(fish_vec))
I_results3 <- numeric(length(fish_vec))


parameters3 = parameters

parameters3["f_J"] = 1
parameters3["f_N"] = 1


for(i in 1:length(fish_vec)){
  #introduction_times = (1:60)*30
  parameters3["Preds"] = fish_vec[i]
  PondSimPreds = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters3, method="lsoda", func=Pond_ODE))
  
  #introduction_times[i] = introduction_time
  L3F_results3[i] <- round(PondSimPreds$L3F[nrow(PondSimPreds)], 5) 
  N_results3[i] <- round(PondSimPreds$N[nrow(PondSimPreds)], 5) 
  J_results3[i] <- round(PondSimPreds$J[nrow(PondSimPreds)], 5) 
  A_results3[i] <- round(PondSimPreds$A[nrow(PondSimPreds)], 5) 
  I_results3[i] <- round(PondSimPreds$I[nrow(PondSimPreds)], 5) 
}


fish_vec = data.frame(fish_vec)
L3s = data.frame(L3F_results3)
Nsim = data.frame(N_results3)
Jsim = data.frame(J_results3)
Asim = data.frame(A_results3)
Isim = data.frame((I_results3))

data3 = cbind(fish_vec,L3s,Nsim,Jsim,Asim,Isim)

data3 = data3 %>% mutate(totalcopes = N_results + J_results + A_results)


L3F = ggplot(data=data3,aes(x=fish_vec,y=L3F_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Total Parasites in Fish")

Infected = ggplot(data=data3,aes(x=fish_vec,y=I_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Infected Copepods")

ggarrange(L3F,Infected, 
          nrow = 1, ncol = 2,
          labels = c("A","B"),
          common.legend = TRUE,
          legend = "right")


datainfected3 = data3[,c(1,2,6)]

datainfectedlong3 = datainfected3 %>% pivot_longer(cols=c(2,3),names_to = "InfectionRoute", values_to = "ParasiteDensity")

datainfectedlong3$InfectionRoute = gsub("L3F_results", "L3F", datainfectedlong3$InfectionRoute)
datainfectedlong3$InfectionRoute = gsub("X.I_results.", "InfectedCope", datainfectedlong3$InfectionRoute)

colors = c(
  L3F = "#77AADD",
  InfectedCope = "#EE8866") 

NonSelectivePredAttack = ggplot(data=datainfectedlong3, aes(x = fish_vec,y=ParasiteDensity,group=InfectionRoute,color=InfectionRoute)) + geom_line(linewidth=2) +
  theme_classic()


#no selectivity in handling time
timespan = 365*5
fish_vec = seq(from=0, to=0.075, by=0.01)
Initial_conditions = c(N = 500, J = 200, A = 25, Exposed_values, I = 0, L3F = 0)
L3F_results4 <- numeric(length(fish_vec))  
N_results4 <- numeric(length(fish_vec))  
J_results4 <- numeric(length(fish_vec))  
A_results4 <- numeric(length(fish_vec))
I_results4 <- numeric(length(fish_vec))


parameters4 = parameters

parameters4["h_J"] = 1
parameters4["h_N"] = 1


for(i in 1:length(fish_vec)){
  #introduction_times = (1:60)*30
  parameters4["Preds"] = fish_vec[i]
  PondSimPreds = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters4, method="lsoda", func=Pond_ODE))
  
  #introduction_times[i] = introduction_time
  L3F_results4[i] <- round(PondSimPreds$L3F[nrow(PondSimPreds)], 5) 
  N_results4[i] <- round(PondSimPreds$N[nrow(PondSimPreds)], 5) 
  J_results4[i] <- round(PondSimPreds$J[nrow(PondSimPreds)], 5) 
  A_results4[i] <- round(PondSimPreds$A[nrow(PondSimPreds)], 5) 
  I_results4[i] <- round(PondSimPreds$I[nrow(PondSimPreds)], 5) 
}


fish_vec = data.frame(fish_vec)
L3s = data.frame(L3F_results4)
Nsim = data.frame(N_results4)
Jsim = data.frame(J_results4)
Asim = data.frame(A_results4)
Isim = data.frame(I_results4)

data4 = cbind(fish_vec,L3s,Nsim,Jsim,Asim,Isim)

data4 = data3 %>% mutate(totalcopes = N_results + J_results + A_results)


L3F = ggplot(data=data4,aes(x=fish_vec,y=L3F_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Total Parasites in Fish")

Infected = ggplot(data=data4,aes(x=fish_vec,y=I_results)) + geom_line(size=2, color = "#77AADD") + theme_classic() +
  labs(x=expression('Fish Density, L' ^ -1),y="Infected Copepods")

ggarrange(L3F,Infected, 
          nrow = 1, ncol = 2,
          labels = c("A","B"),
          common.legend = TRUE,
          legend = "right")


datainfected4 = data4[,c(1,2,6)]

datainfectedlong4 = datainfected3 %>% pivot_longer(cols=c(2,3),names_to = "InfectionRoute", values_to = "ParasiteDensity")

datainfectedlong4$InfectionRoute = gsub("L3F_results", "L3F", datainfectedlong4$InfectionRoute)
datainfectedlong4$InfectionRoute = gsub("X.I_results.", "InfectedCope", datainfectedlong4$InfectionRoute)

colors = c(
  L3F = "#77AADD",
  InfectedCope = "#EE8866") 

NonSelectivePredHandling = ggplot(data=datainfectedlong4, aes(x = fish_vec,y=ParasiteDensity,group=InfectionRoute,color=InfectionRoute)) + geom_line(linewidth=2) +
  theme_classic() 



ggarrange(SelectivePred,NonSelectivePred,NonSelectivePredAttack,NonSelectivePredHandling,
          nrow = 2, ncol = 2,
          labels = c("Selective","NonSelective","SelectiveAttack","SelectiveHandling"),
          common.legend = TRUE,
          legend = "right")



ggarrange(SelectivePred,NonSelectivePred,
          nrow = 2, ncol = 1,
          labels = c("Selective","Non-Selective"),
          common.legend = TRUE,
          legend = "right",
          label.x = 0.1,
          align = "v")



