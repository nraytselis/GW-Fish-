#Packages
library(panelPomp)
library(adaptMCMC)

print("pt 1")
##### Imports current best fit parameters and jump covariances #####
# Gets the best fit parameters and the covariance matrix
get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump)
  
}

setwd("~/Desktop/Rscripts/Data/Fish")
chainsA <- readRDS("Joint_GW_full_A2.RDA")
chainsB <- readRDS("Joint_GW_full_B2.RDA")
chainsC <- readRDS("Joint_GW_full_C2.RDA")

fishchains = c(chainsA,chainsB,chainsC)

pars = get_best_fit(fishchains)[[1]]
variances = get_best_fit(fishchains)[[2]] 

###############################################################

##### Brings in data for fish experiment #####
# Brings in data, note this has been processed in other scripts to avoid using tidyverse on cluster
Fish_Time_Series = read.csv("Fish_experiment_data.csv")
###############################################################

##### Pomp construction Fish experiment #####
# Initial conditions function for pomp
GW_F_rinit <- Csnippet("
  N = rpois(N0_F);
  J = rpois(J0_F);
  A = rpois(A0_F);
")

# r process model
GW_F_step_process <- Csnippet("
  double VOL = 40;
  
  /* These three equations specify predation rates on A,J,N respectiviely. f is predation rate on adults. f_j and f_n, are relative predation rates on juveniles and nauplii.*/
  /* Similarly, h is handeling time. h_J and h_N are relative handeling times. Lastly, i_P representats strength of interference.*/
  double Pred_A = f*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  double Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  double Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J + f_N*h_N*N)/VOL  + i_P*fmax(Preds/VOL - 1/VOL, 0)/VOL);
  
  double birth_rate = b_M*exp(-comp_b / VOL * (c_N * N + c_J * J + A));

  double d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A));
  double d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + cann*A/VOL;
  
  double m_N_c = m_N*exp(-comp_m / VOL * (c_N * N + c_J * J + A));
  double m_J_c = m_J*exp(-comp_m / VOL * (c_N * N + c_J * J + A));


  /* Specify adult rates */
  const double A_rates[2] = {d_A_c, Pred_A};
  double A_trans[2];
  reulermultinom(2, A, &A_rates, dt, &A_trans);
  
  /* Specify juvenile rates */
  const double J_rates[3] = {m_J_c, d_J_c, Pred_J};
  double J_trans[3];
  reulermultinom(3, J, &J_rates, dt, &J_trans);
  
  /* Specify nauplii rates */
  const double N_rates[3] = {m_N_c, d_N_c, Pred_N};
  double N_trans[3];
  reulermultinom(3, N, &N_rates, dt, &N_trans);
  
  /* births */
  int births = rpois(birth_rate * A / 2 * dt);
  
  /* New totals */
  int A_out = A - A_trans[0] - A_trans[1] + J_trans[0];
  int J_out = J - J_trans[0] - J_trans[1] - J_trans[2] + N_trans[0];
  int N_out = N - N_trans[0] - N_trans[1] - N_trans[2] + births;
  N = N_out >= 0 ? N_out : 0;
  J = J_out >= 0 ? J_out : 0;
  A = A_out >= 0 ? A_out : 0;
")

# r measure model - for some reason, Pomp will not compile when we make this a C-snippet (all others work just fine)
GW_F_rmeasure = function(N,J,A,VOL=40,SAMPLEVOL=0.5, counted_volume, aliquot_volume=10, k_F, ...){
  c(
    N1 = rnbinom(n=1, mu=N*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F), 
    JOA1 = rnbinom(n=1, mu=(J + 0.5*(1 + 1/3)*A)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    AF1 = rnbinom(n=1, mu=0.5*2/3*A*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    N2 = rnbinom(n=1, mu=N*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F), 
    JOA2 = rnbinom(n=1, mu=(J + 0.5*(1 + 1/3)*A)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F),
    AF2 = rnbinom(n=1, mu=0.5*2/3*A*counted_volume/aliquot_volume*SAMPLEVOL/VOL, size=k_F))
}

# d measure model
GW_F_dmeasure = Csnippet("
  double VOL = 40;
  double SAMPLEVOL = 0.5;
  double aliquot_volume = 10;
  lik = dnbinom_mu(N1, k_F,fmax(N,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(JOA1, k_F,  fmax(J+ 0.5*(1 + 1.0/3.0)*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(AF1, k_F, fmax(0.5*2.0/3.0*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(N2, k_F, fmax(N,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(JOA2, k_F, fmax(J+ 0.5*(1 + 1.0/3.0)*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log) +
        dnbinom_mu(AF2, k_F, fmax(0.5*2.0/3.0*A,1)*counted_volume/aliquot_volume*SAMPLEVOL/VOL, give_log);"
)
###############################################################

##### Panel pomp construction for Fish experiment #####

# Template pomp to eventually construct the panelpomp
template_F = pomp(
  data =data.frame( #the template structure needs a dummy data set with the right column names
    Day = c(28, 42, 49, 56, 63, 70, 77, 84, 91, 98), #days are based on experimental sampling days
    Tank = NA, 
    total_observedN = NA, total_JA = NA), times="Day", t0=28,
  rprocess = discrete_time(GW_F_step_process, delta.t=1),
  rmeasure = GW_F_rmeasure,
  dmeasure = GW_F_dmeasure,
  rinit = GW_F_rinit,
  obsnames = c("N1", "JOA1", "AF1", "N2", "JOA2", "AF2"),
  statenames = c("N","J","A"),
  covarnames = c("Preds", "counted_volume"),
  paramnames = c("N0_F", "J0_F", "A0_F", "f", "f_J", "f_N", "h", "h_J", "h_N",  "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", 
                 "k_F", "i_P"),
  params = pars[c("N0_F", "J0_F", "A0_F", "f", "f_J", "f_N", "h", "h_J", "h_N",  "b_M", "comp_b", "comp_d", "comp_m", "cann", "c_N", "c_J", "d_A", "m_J", "d_J", "m_N", "d_N", 
                  "k_F", "i_P")]
)

# Checking how many tanks there are
U <- length(unique(Fish_Time_Series$Tank))
# Empty list for the pomp templates
poList <- setNames(vector(mode = "list", length = U),
                   nm = paste0("unit", 1:U))

# This loop populates each element of the list with the model template, covariates, and the data (Fish_Time_Series)
for (u in seq_len(U)) {
  # For each tank, we need to grab the covariates, add the starting condition covariates, and get the covariate timing right
  #                for fish introduction
  covariates <- subset(Fish_Time_Series, Tank == u, select=c("Day", "Live_Fish", "counted_volume"))
  Pred_treat <- subset(Fish_Time_Series, Tank == u, select=c("Preds"))
  covariates <- rbind(covariates, c(50, max(Pred_treat, 0), 0)) # Get covariates for day of fish addition
  colnames(covariates)[1:2] <- c("FishDay", "Preds")  # time for covariates cannot have the same name as time for observed data
  # orders covariate table by time, which pomp() wants
  covariates = covariates[order(covariates$FishDay),]
  cov_table = covariate_table(covariates, times="FishDay")
  poList[[u]] <- pomp(template_F, covar = cov_table)
  data_u <- pomp(
    data = subset(Fish_Time_Series, Tank == u, select = c("Day", "Tank", "N1", "JOA1", "AF1", "N2", "JOA2", "AF2")),
    times="Day",
    t0=timezero(template_F)
  )
  poList[[u]]@data <- data_u@data
}

# Creates the panel Pomp object
GW_F_panel <- panelPomp(object = poList, shared = coef(template_F))
#plot(simulate(GW_F_panel, nsim = 1))
###############################################################


###testing model fit to best parameters
params = list()
for(l in 1:24){
  for(i in (5:25)*10) 
    params = c(params, list(fishchains[[l]]$samples[i,]))
}

post_sim = function(x){simulate(GW_F_panel, shared=x, nsim=1)} 

sim = lapply(X=params,FUN=post_sim) 

#sim = simulate(GW_panel, nsim = 1000)

# Create an empty list to hold data frames from each simulation
df_list <- vector("list", length(sim))

# Loop over each simulation (1-100)
for (s in seq_along(sim)) {
  sim_obj <- sim[[s]]@unit_objects
  
  # Loop over units in the simulation (1-60)
  unit_dfs <- lapply(seq_along(sim_obj), function(u) {
    unit_obj <- sim_obj[[u]]
    mat <- unit_obj@states  # 3 x 18 matrix (latent states - true abundance)
    
    # Create a data frame with 18 rows (one per time step aka week)
    data.frame(
      time = unit_obj@times,
      sim = s,
      unit = u,
      N = mat["N", ],
      J = mat["J", ],
      A = mat["A", ]
    )
  }) 
  
  # Combine all 60 units for this simulation into one data frame
  df_list[[s]] <- do.call(rbind, unit_dfs)
}


# Combine all simulations into one big data frame
full_df <- do.call(rbind, df_list)

#make new columns for AF1, JOA1,N1
full_df$N1 <- full_df$N
full_df$JOA1 <- round(full_df$J + 2/3*(full_df$A))
full_df$AF1 <- round(1/3*(full_df$A)) 

#add column for total copepods
full_df$total = full_df$AF1 + full_df$JOA1 + full_df$N1

#bring in fish treatment info
Fishtreatments = read_csv("FishTreatments.csv") 

#merge dataframes with treatment info and sim
merged_df <- left_join(full_df, Fishtreatments, by = "unit")

#summary stats 
merged_df_summary = merged_df %>% group_by(FishDensity,time) %>% summarise(mean_A = mean(AF1), 
                                                                                          mean_J = mean(JOA1), 
                                                                                          mean_N = mean(N1),
                                                                                          mean_total = mean(total),
                                                                                          sd_A = sd(AF1),
                                                                                          sd_J = sd(JOA1),
                                                                                          sd_N = sd(N1),
                                                                                          sd_total = sd(total),
                                                                                          n_A = n(),
                                                                                          n_J = n(),
                                                                                          n_N = n(),
                                                                                          n_total = n(),
                                                                                          se_A = sd_A / sqrt(n_A),
                                                                                          se_J = sd_J / sqrt(n_J),
                                                                                          se_N = sd_N / sqrt(n_N),
                                                                                          se_total = sd_total / sqrt(n_total),
                                                                                          lower_ci_A = quantile(AF1, prob=0.025),
                                                                                          upper_ci_A = quantile(AF1, prob=0.975),
                                                                                          lower_ci_J = quantile(JOA1, prob=0.025),
                                                                                          upper_ci_J = quantile(JOA1, prob=0.975),
                                                                                          lower_ci_N =quantile(N1, prob=0.025),
                                                                                          upper_ci_N = quantile(N1, prob=0.975),
                                                                                          lower_ci_total =quantile(total, prob=0.025)) 
  
                                                                           
                                                                           
                                                                           
####visualizations                                                                                                                                                                   upper_ci_total = quantile(total, prob=0.975)) 
#facet wrapped

A = ggplot(merged_df_summary, aes(x = time, y = mean_A)) +
  geom_line() +
  facet_wrap(~ FishDensity) +
  theme_minimal() +
  labs(title = "Mean A over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_A, ymax=upper_ci_A, alpha =0.2),fill = "grey") + theme_classic()

J = ggplot(merged_df_summary, aes(x = time, y = mean_J)) +
  geom_line() +
  facet_wrap(~ FishDensity) +
  theme_minimal() +
  labs(title = "Mean J over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_J, ymax=upper_ci_J, alpha =0.2),fill = "grey") + theme_classic()

N = ggplot(merged_df_summary, aes(x = time, y = mean_N)) +
  geom_line() +
  facet_wrap(~ FishDensity) +
  theme_minimal() +
  labs(title = "Mean N over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_N, ymax=upper_ci_N, alpha =0.2),fill = "grey") + theme_classic()

total = ggplot(merged_df_summary, aes(x = time, y = mean_total)) +
  geom_line() +
  facet_wrap(~ FishDensity) +
  theme_minimal() +
  labs(title = "Mean Total Copepods over Time by Treatment Group") +
  geom_ribbon(aes(ymin=lower_ci_total, ymax=upper_ci_total, alpha =0.2),fill = "grey") + theme_classic()

library(ggpubr)
ggarrange(total,A,J,N,
          nrow = 2,
          ncol =2,
          labels = c("A","B","C","D"), 
          legend = "none"
)


#all on same plot
A = ggplot(merged_df_summary, aes(x = time, y = mean_A, group = as.factor(FishDensity), color = as.factor(FishDensity))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Mean A over Time by Treatment Group") 

J = ggplot(merged_df_summary, aes(x = time, y = mean_J,group = as.factor(FishDensity), color = as.factor(FishDensity))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Mean J over Time by Treatment Group") 

N = ggplot(merged_df_summary, aes(x = time, y = mean_N, group = as.factor(FishDensity), color = as.factor(FishDensity))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Mean N over Time by Treatment Group") 
total = ggplot(merged_df_summary, aes(x = time, y = mean_total,group = as.factor(FishDensity), color = as.factor(FishDensity))) +
  geom_line() + 
  theme_minimal() +
  labs(title = "Mean Total Copepods over Time by Treatment Group") 
library(ggpubr)
ggarrange(total,A,J,N,
          nrow = 2,
          ncol =2,
          labels = c("A","B","C","D"), 
          legend = "none"
)


