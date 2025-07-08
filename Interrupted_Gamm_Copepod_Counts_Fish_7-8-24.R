library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(wesanderson)
library(scales)
library(dplyr)
library(readr)
library(ggplot2)
library(mgcv)
library(deSolve)
library(tidyverse)
library(stats)
library(itsadug)
library(car)
library(broom)
library(viridis)
library(RColorBrewer)
library(pomp)
library(stringr)


setwd("~/Desktop/Rscripts/Data")

CopepodCounts <- read_csv("CopepodCounts_Experiment1_NR.csv",locale=locale(encoding="latin1"))
CopepodCounts[which(CopepodCounts[,"Fish_Treatment"] == 11),"Fish_Treatment"] = 12

#Does fish treatment have an effect on copepod densities and size structure?

CopepodCounts$Fish_Treatment <- factor(CopepodCounts$Fish_Treatment)
CopepodCounts$Tank <- as.factor(CopepodCounts$Tank)


# #add columns for time since intervention and intervention

#assign values to Time Since Intervention
CopepodCounts <- CopepodCounts %>%
  mutate(Time_Since_Intervention = case_when(
    Week == 2 ~ 0,
    Week == 3 ~ 0,
    Week == 4 ~ 0,
    Week == 5 ~ 1,
    Week == 6 ~ 2,
    Week == 7 ~ 3,
    Week == 8 ~ 4,
    Week == 9 ~ 5,
    Week == 10 ~ 6,
    Week == 11 ~ 7,
    TRUE ~ NA_integer_  # Ensure that other dates get NA in Week
  ))

CopepodCounts <- CopepodCounts %>%
  mutate(Time_Since_Intervention = case_when(
    Fish_Treatment == 0 ~ 0,
    Fish_Treatment != 0 ~ Time_Since_Intervention,
    TRUE ~ NA_integer_  # Ensure that other dates get NA in Week
  ))

#assign values to Intervention
#same as Fish_Treatment after week 4
CopepodCounts <- CopepodCounts %>%
  mutate(Intervention = case_when(
    Week <= 4 ~ 0,
    Week > 4 & Fish_Treatment == 0 ~ 1,
    Week > 4 & Fish_Treatment == 1 ~ 2,  
    Week > 4 & Fish_Treatment == 2 ~ 3, 
    Week > 4 & Fish_Treatment == 4 ~ 4, 
    Week > 4 & Fish_Treatment == 7 ~ 5, 
    Week > 4 & Fish_Treatment == 12 ~ 6, 
    TRUE ~ NA_real_  # Default case
  ))


# Convert Intervention to a factor
CopepodCounts <- CopepodCounts %>%
 mutate(Intervention = factor(Intervention, levels = 0:6))

#make sure data columns in correct format
CopepodCounts$Tank = as.factor(CopepodCounts$Tank)

# Remove or impute NA/NaN/Inf values if found
CopepodCounts <- CopepodCounts %>%
  filter(!is.na(Week) & !is.na(Intervention) & !is.na(`Copepods/L`) & !is.na(`Tank`))

#remove column with all NAs
CopepodCounts <- subset(CopepodCounts, select = -12)

#add additional columns for smooths for each treatment

#fish treatment 1 
CopepodCounts <- CopepodCounts %>%
  mutate(F1 = case_when(
    Week <= 4 ~ 0,
    Week == 5 & Fish_Treatment == 1 ~ 1,
    Week == 6 & Fish_Treatment == 1 ~ 2,  
    Week == 7 & Fish_Treatment == 1 ~ 3, 
    Week == 8 & Fish_Treatment == 1 ~ 4, 
    Week == 9 & Fish_Treatment == 1 ~ 5, 
    Week == 10 & Fish_Treatment == 1 ~ 6,
    Week == 11 & Fish_Treatment == 1 ~ 7,
    Fish_Treatment != 1 ~ 0,
    TRUE ~ NA_real_  # Default case
  ))

#fish treatment 2
CopepodCounts <- CopepodCounts %>%
  mutate(F2 = case_when(
    Week <= 4 ~ 0,
    Week == 5 & Fish_Treatment == 2 ~ 1,
    Week == 6 & Fish_Treatment == 2 ~ 2,  
    Week == 7 & Fish_Treatment == 2 ~ 3, 
    Week == 8 & Fish_Treatment == 2 ~ 4, 
    Week == 9 & Fish_Treatment == 2 ~ 5, 
    Week == 10 & Fish_Treatment == 2 ~ 6,
    Week == 11 & Fish_Treatment == 2 ~ 7,
    Fish_Treatment != 2 ~ 0,
    TRUE ~ NA_real_  # Default case
  ))

#fish treatment 4
CopepodCounts <- CopepodCounts %>%
  mutate(F4 = case_when(
    Week <= 4 ~ 0,
    Week == 5 & Fish_Treatment == 4 ~ 1,
    Week == 6 & Fish_Treatment == 4 ~ 2,  
    Week == 7 & Fish_Treatment == 4 ~ 3, 
    Week == 8 & Fish_Treatment == 4 ~ 4, 
    Week == 9 & Fish_Treatment == 4 ~ 5, 
    Week == 10 & Fish_Treatment == 4 ~ 6,
    Week == 11 & Fish_Treatment == 4 ~ 7,
    Fish_Treatment != 4 ~ 0,
    TRUE ~ NA_real_  # Default case
  ))

#fish treatment 7
CopepodCounts <- CopepodCounts %>%
  mutate(F7 = case_when(
    Week <= 4 ~ 0,
    Week == 5 & Fish_Treatment == 7 ~ 1,
    Week == 6 & Fish_Treatment == 7 ~ 2,  
    Week == 7 & Fish_Treatment == 7 ~ 3, 
    Week == 8 & Fish_Treatment == 7 ~ 4, 
    Week == 9 & Fish_Treatment == 7 ~ 5, 
    Week == 10 & Fish_Treatment == 7 ~ 6,
    Week == 11 & Fish_Treatment == 7 ~ 7,
    Fish_Treatment != 7 ~ 0,
    TRUE ~ NA_real_  # Default case
  ))


#fish treatment 12
CopepodCounts <- CopepodCounts %>%
  mutate(F12 = case_when(
    Week <= 4 ~ 0,
    Week == 5 & Fish_Treatment == 12 ~ 1,
    Week == 6 & Fish_Treatment == 12 ~ 2,  
    Week == 7 & Fish_Treatment == 12 ~ 3, 
    Week == 8 & Fish_Treatment == 12 ~ 4, 
    Week == 9 & Fish_Treatment == 12 ~ 5, 
    Week == 10 & Fish_Treatment == 12 ~ 6,
    Week == 11 & Fish_Treatment == 12 ~ 7,
    Fish_Treatment != 12 ~ 0,
    TRUE ~ NA_real_  # Default case
  ))

#Add column with sum of all copepod classes
CopepodCounts$TotalCopes <- CopepodCounts$`Adult with Eggs` + CopepodCounts$`Copepodite No eggs` + CopepodCounts$Nauplii 

#change volume column name
colnames(CopepodCounts)[8] = "vol"

colnames(CopepodCounts)[10] = "CopesPerL"

#interrupted GAMM
TotalCopepodsInterrupted <- gam(
  round(CopesPerL) ~ s(Week) + s(Tank, bs="re") + s(F1, k = 5) + s(F2, k = 5) + s(F4, k = 5) + s(F7, k = 5) + s(F12, k = 5),
  family = nb(),
  data = CopepodCounts
)

summary(TotalCopepodsInterrupted)

#Make a data frame for treatments I want to predict for
#make note of 1 tank for each treatment (tank 1,2,3,4,5,8)
#remove sample B because we assume that all samples have the same treatment mean
df_for_predictions = subset(CopepodCounts, Tank %in% c(1,2,3,4,5,8) & str_detect(Sample, "A")) 

### predict function to get the model fit 
pred <- data.frame(predict(TotalCopepodsInterrupted, newdata=df_for_predictions, type = 'link', se.fit = TRUE, exclude="s(Tank)")) 
pred$ll = exp(pred$fit - pred$se.fit)
pred$ul = exp(pred$fit + pred$se.fit)
pred$fit = exp(pred$fit)

manage_cope_density = cbind(df_for_predictions, pred)

#summarize original data to get treatment means to plot
CopepodCountsSummary = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(CopesPerL), SE = sd(CopesPerL) / sqrt(n()))

CopesDensityplot <- ggplot(data=manage_cope_density, aes(x=Week, y=fit, group=Fish_Treatment, colour=Fish_Treatment)) +
  theme_classic() + 
  geom_line(data = manage_cope_density, aes(x = Week, y = fit, group = Fish_Treatment, colour = Fish_Treatment, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density, aes(ymin = ll, ymax = ul, fill = Fish_Treatment, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummary, aes(x=Week, y= mean, colour=Fish_Treatment), inherit.aes=F, show.legend = FALSE) +
  labs(y = "Copepods Per L") +
  geom_errorbar(data=CopepodCountsSummary, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=Fish_Treatment), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Treatment"), alpha = "none")


print(CopesDensityplot)

#recode fish treatment as fish density 
manage_cope_density$FishPerL <- as.numeric(as.character(manage_cope_density$Fish_Treatment)) / 40
manage_cope_density$FishPerL <- as.factor(manage_cope_density$FishPerL)

CopepodCountsSummary$FishPerL <- as.numeric(as.character(CopepodCountsSummary$Fish_Treatment)) / 40
CopepodCountsSummary$FishPerL <- as.factor(CopepodCountsSummary$FishPerL)


ggplot(data=manage_cope_density, aes(x=Week, y=fit, group=Fish_Treatment, colour=Fish_Treatment)) +
  theme_classic() + 
  geom_line(data = manage_cope_density, aes(x = Week, y = fit, group = Fish_Treatment, colour = Fish_Treatment, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density, aes(ymin = ll, ymax = ul, fill = Fish_Treatment, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummary, aes(x=Week, y= mean, colour=Fish_Treatment), inherit.aes=F, show.legend = FALSE) +
  labs(y = "Copepods Per L") +
  geom_errorbar(data=CopepodCountsSummary, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=Fish_Treatment), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Density"), alpha = "none") +
  scale_fill_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +   scale_color_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15))  

#fish density 
ggplot(data=manage_cope_density, aes(x=Week, y=fit, group=FishPerL, colour=FishPerL)) +
  theme_classic() + 
  geom_line(data = manage_cope_density, aes(x = Week, y = fit, group = FishPerL, colour = FishPerL, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density, aes(ymin = ll, ymax = ul, fill = FishPerL, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummary, aes(x=Week, y= mean, colour=FishPerL), inherit.aes=F, show.legend = FALSE) +
  labs(y = "Copepods Per L") +
  geom_errorbar(data=CopepodCountsSummary, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=FishPerL), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Density"), alpha = "none") +
  scale_fill_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +   scale_color_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15))  


###add counts per liter for other stages
#divide by counted volume and multiply by 20
colnames(CopepodCounts)[5] <- "AF"
colnames(CopepodCounts)[6] <- "JOA"

CopepodCounts$NperL = 20*(CopepodCounts$Nauplii/CopepodCounts$vol) 
CopepodCounts$JOAperL = 20*(CopepodCounts$JOA/CopepodCounts$vol) 
CopepodCounts$AFperL = 20*(CopepodCounts$AF/CopepodCounts$vol) 

NaupliiCountsSummary = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(NperL), SE = sd(NperL) / sqrt(n()))
JOACountsSummary = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(JOAperL), SE = sd(JOAperL) / sqrt(n()))
AFCountsSummary = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(AFperL), SE = sd(AFperL) / sqrt(n()))


###Run for Nauplii 
#interrupted GAMM Naup
NInterrupted <- gam(
  round(NperL) ~ s(Week) + s(Tank, bs="re") + s(F1, k = 5) + s(F2, k = 5) + s(F4, k = 5) + s(F7, k = 5) + s(F12, k = 5),
  family = nb(),
  data = CopepodCounts
)

df_for_predictions_N = subset(CopepodCounts, Tank %in% c(1,2,3,4,5,8) & str_detect(Sample, "A")) 

### predict function to get the model fit 
pred <- data.frame(predict(NInterrupted, newdata=df_for_predictions_N, type = 'link', se.fit = TRUE, exclude="s(Tank)")) 
pred$ll = exp(pred$fit - pred$se.fit)
pred$ul = exp(pred$fit + pred$se.fit)
pred$fit = exp(pred$fit)

manage_cope_density_N = cbind(df_for_predictions_N, pred)

#summarize original data to get treatment means to plot
CopepodCountsSummaryNaup = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(NperL), SE = sd(NperL) / sqrt(n()))

ggplot(data=manage_cope_density_N, aes(x=Week, y=fit, group=Fish_Treatment, colour=Fish_Treatment)) +
  theme_classic() + 
  geom_line(data = manage_cope_density_N, aes(x = Week, y = fit, group = Fish_Treatment, colour = Fish_Treatment, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density_N, aes(ymin = ll, ymax = ul, fill = Fish_Treatment, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummaryNaup, aes(x=Week, y= mean, colour=Fish_Treatment), inherit.aes=F, show.legend = FALSE) +
  labs(y = "Nauplii Per L") +
  geom_errorbar(data=CopepodCountsSummaryNaup, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=Fish_Treatment), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Density"), alpha = "none") +
  scale_fill_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +   scale_color_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15))  

###Run for J0A
#interrupted GAMM Naup
JOAInterrupted <- gam(
  round(JOAperL) ~ s(Week) + s(Tank, bs="re") + s(F1, k = 5) + s(F2, k = 5) + s(F4, k = 5) + s(F7, k = 5) + s(F12, k = 5),
  family = nb(),
  data = CopepodCounts
)

df_for_predictions_JOA = subset(CopepodCounts, Tank %in% c(1,2,3,4,5,8) & str_detect(Sample, "A")) 

### predict function to get the model fit 
pred <- data.frame(predict(JOAInterrupted, newdata=df_for_predictions_JOA, type = 'link', se.fit = TRUE, exclude="s(Tank)")) 
pred$ll = exp(pred$fit - pred$se.fit)
pred$ul = exp(pred$fit + pred$se.fit)
pred$fit = exp(pred$fit)

manage_cope_density_JOA = cbind(df_for_predictions_JOA, pred)

#summarize original data to get treatment means to plot
CopepodCountsSummaryJOA = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(JOAperL), SE = sd(JOAperL) / sqrt(n()))

ggplot(data=manage_cope_density_JOA, aes(x=Week, y=fit, group=Fish_Treatment, colour=Fish_Treatment)) +
  theme_classic() + 
  geom_line(data = manage_cope_density_JOA, aes(x = Week, y = fit, group = Fish_Treatment, colour = Fish_Treatment, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density_JOA, aes(ymin = ll, ymax = ul, fill = Fish_Treatment, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummaryJOA, aes(x=Week, y= mean, colour=Fish_Treatment), inherit.aes=F, show.legend = FALSE) +
  labs(y = "JOA Per L") +
  geom_errorbar(data=CopepodCountsSummaryJOA, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=Fish_Treatment), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Density"), alpha = "none") +
  scale_fill_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +   scale_color_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15))  

###Run for AF
#interrupted GAMM Naup
AFInterrupted <- gam(
  round(AFperL) ~ s(Week) + s(Tank, bs="re") + s(F1, k = 5) + s(F2, k = 5) + s(F4, k = 5) + s(F7, k = 5) + s(F12, k = 5),
  family = nb(),
  data = CopepodCounts
)

df_for_predictions_AF = subset(CopepodCounts, Tank %in% c(1,2,3,4,5,8) & str_detect(Sample, "A")) 

### predict function to get the model fit 
pred <- data.frame(predict(AFInterrupted, newdata=df_for_predictions_AF, type = 'link', se.fit = TRUE, exclude="s(Tank)")) 
pred$ll = exp(pred$fit - pred$se.fit)
pred$ul = exp(pred$fit + pred$se.fit)
pred$fit = exp(pred$fit)

manage_cope_density_AF = cbind(df_for_predictions_AF, pred)

#summarize original data to get treatment means to plot
CopepodCountsSummaryAF = CopepodCounts %>% group_by(Week,Fish_Treatment) %>% summarise(mean = mean(AFperL), SE = sd(AFperL) / sqrt(n()))

ggplot(data=manage_cope_density_AF, aes(x=Week, y=fit, group=Fish_Treatment, colour=Fish_Treatment)) +
  theme_classic() + 
  geom_line(data = manage_cope_density_AF, aes(x = Week, y = fit, group = Fish_Treatment, colour = Fish_Treatment, show.legend = FALSE)) +
  geom_ribbon(data=manage_cope_density_AF, aes(ymin = ll, ymax = ul, fill = Fish_Treatment, colour = NA, alpha = 0.2, show.legend = TRUE)) + 
  geom_point(data=CopepodCountsSummaryAF, aes(x=Week, y= mean, colour=Fish_Treatment), inherit.aes=F, show.legend = FALSE) +
  labs(y = "AF Per L") +
  geom_errorbar(data=CopepodCountsSummaryAF, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=Fish_Treatment), width=0.2, inherit.aes=FALSE, show.legend = FALSE) +
  guides(colour = "none", fill = guide_legend(title = "Fish Density"), alpha = "none") +
  scale_fill_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +   scale_color_manual( values = c("#648FFF","#785EF0","#DC267F","#FE6100","#FFB000","#DDCC77")) +
  theme(axis.text = element_text(size = 13)) + theme(axis.title = element_text(size = 15)) 
