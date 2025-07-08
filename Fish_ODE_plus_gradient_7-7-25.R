# Loading packages
library(Matrix)
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(deSolve)


Pond_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]; Preds = y[5+latent_stages]; L3F = y[6+latent_stages]
    VOL = 1
    
    Pred_A = f*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*Preds/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es))) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es))) #density dependence in maturation
    
    dNdt = b_M*(A + sum(Es))/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A + sum(Es))) - (m_N_c+d_N_c)*N - cann*(A + I + sum(Es))*N - Pred_N*N
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J -Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A - lambda*A - Pred_A*A
    
    # development of all stages
    latent_progression = latent_rate*Es
    # lost to next stage   #death      #gained from last stage
    dEsdt = -latent_progression - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)]) - Pred_A*Es
    
    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I - Pred_A*I
    
    dPredsdt = 0 #convEff*(Pred_N*N + Pred_J*J + Pred_A*A + Pred_A*sum(Es)) - d_F*Preds 
    
    dL3Fdt <- Pred_A*I - d_W*L3F - d_F*L3F 
    
    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt,dPredsdt, dL3Fdt)
    
    
    return(list(result))
  }
  )
}

#bring in mcmc chains
setwd("~/Desktop/Rscripts/Data")
joint_fit <- readRDS(file = "Joint_GW_50000p.RDA")

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

joint_fit = get_best_fit(joint_fit)

parameters = joint_fit[[1]]

parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 1.3
parameters["d_W"] = 0.05
parameters["d_F"] = 0.05
parameters["convEff"] = 0.001 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)

parameters = unlist(parameters)

Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values


Initial_conditions = c(N = 500, J = 200, A = 25, Exposed_values, I = 0, Preds = 0.05,L3F = 0)
timespan = 365*50

PondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                         method="lsoda", func=Pond_ODE)) 



#reformat sim data frame
PondSim[,"Es"] = rowSums(PondSim) - PondSim[,"N"] - PondSim[,"J"]- PondSim[,"A"] - PondSim[,"I"] - PondSim[,"time"]
PondSim[,"Exposed"] = apply(X=PondSim[,which(str_detect(colnames(PondSim),"E"))],MARGIN=1,FUN = sum) 

#check to make sure that events worked, NEED TO CHANGE EVENT SO THAT CANNOT GO NEGATIVE 
ggplot(PondSim, aes(x = time)) +
  geom_line(aes(y = Preds), color = "blue") +
  geom_line(aes(y = L3F), color = "red") +
  labs(y = "Fish (blue = total, red = infected)", title = "Predator Fish Dynamics Over Time") +
  theme_minimal()

#reformat for plotting
PondSim = PondSim %>% select(time,N,J,A,Exposed,I,Preds, L3F)
PondSim = PondSim %>% pivot_longer(cols = c(N,J,A,Exposed,I,Preds, L3F)) 

#plot sim 
p1 = ggplot(PondSim, aes(x=time, y = value + 0.01 , group = name, color = name)) + 
  geom_line() + ylab("density per L") + theme_minimal() + scale_y_log10()
p1



###Pond sim over a fish density gradient 

timespan = 365*10
fish_vec = seq(from=0, to=1, by=0.01)
L3F_results <- numeric(length(fish_vec))  


for(i in 1:length(fish_vec)) {
  Initial_conditions <- c(
    N = 500,
    J = 200,
    A = 25,
    Exposed_values,
    I = 0,
    Preds = fish_vec[i],
    L3F = 0
  )
  
  PondSimPreds <- data.frame(ode(
    y = Initial_conditions,
    times = 1:timespan,
    parms = parameters,
    method = "lsoda",
    func = Pond_ODE
  ))
  
  # Assuming you want the last time point value of L3F
  L3F_results[i] <- round(PondSimPreds$L3F[nrow(PondSimPreds)], 5) 
}

plot(fish_vec, L3F_results, type = "l", xlab = "Fish Density", ylab = "Final L3F", main = "L3F vs Fish Density")

fish_vec = data.frame(fish_vec)
L3s = data.frame(L3F_results)

data = cbind(fish_vec,L3s)

ggplot(data=data, aes(x = fish_vec, y = L3F_results)) + geom_line(linewidth=1) + theme_classic() + xlab("Fish Density (per L)") + ylab("L3 Parasites in Fish (per Fish)") +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15))







