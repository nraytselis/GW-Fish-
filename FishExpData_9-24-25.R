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
library(ggpubr)


setwd("~/Desktop/Rscripts/Data")

CopepodCounts <- read_csv("CopepodCounts_Experiment1_NR.csv",locale=locale(encoding="latin1"))

#remove rows with missing values for tank 18
CopepodCounts <- CopepodCounts[CopepodCounts$Initials %in% c("NR", "MR","JM","TH","AS"), ]

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

CopepodCounts$FishPerL <- as.numeric(as.character(CopepodCounts$Fish_Treatment)) / 40
CopepodCounts$FishPerL <- as.factor(CopepodCounts$FishPerL)

CopepodCountsSummary = CopepodCounts %>% group_by(Week,FishPerL) %>% summarise(mean = mean(CopesPerL), SE = sd(CopesPerL) / sqrt(n()))

colnames(CopepodCounts)[5] <- "AF"
colnames(CopepodCounts)[6] <- "JOA"
CopepodCounts$NperL = 20*(CopepodCounts$Nauplii/CopepodCounts$vol) 
CopepodCounts$JOAperL = 20*(CopepodCounts$JOA/CopepodCounts$vol) 
CopepodCounts$AFperL = 20*(CopepodCounts$AF/CopepodCounts$vol) 

CopepodCountsSummaryNaup = CopepodCounts %>% group_by(Week,FishPerL) %>% summarise(mean = mean(NperL), SE = sd(NperL) / sqrt(n()))
CopepodCountsSummaryJOA = CopepodCounts %>% group_by(Week,FishPerL) %>% summarise(mean = mean(JOAperL), SE = sd(JOAperL) / sqrt(n()))
CopepodCountsSummaryAF = CopepodCounts %>% group_by(Week,FishPerL) %>% summarise(mean = mean(AFperL), SE = sd(AFperL) / sqrt(n()))

combined_palette <- c("#77AADD","#EE8866","#EEDD88","#44BB99","#BBCC33","#DDDDDD") 

total = ggplot(data=CopepodCountsSummary, aes(x=Week, y= mean, colour=FishPerL), inherit.aes=F, show.legend = FALSE) + geom_line(linewidth=1.5) + geom_point(size = 3) + 
  geom_errorbar(data=CopepodCountsSummary, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=FishPerL), width=0.2, inherit.aes=FALSE, show.legend = TRUE) +
  theme_classic() + labs(y = expression('Total Copepod Density, L' ^ -1)) +
  guides(color = guide_legend(title = expression('Fish Density, L' ^ -1),  ))

N = ggplot(data=CopepodCountsSummaryNaup, aes(x=Week, y= mean, colour=FishPerL), inherit.aes=F, show.legend = FALSE) +
  labs(y = expression('N Density, L' ^ -1)) + geom_line(linewidth=1.5) + geom_point(size = 3) + 
  geom_errorbar(data=CopepodCountsSummaryNaup, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=FishPerL), width=0.2, inherit.aes=FALSE, show.legend = TRUE) +
  theme_classic() + guides(color = guide_legend(title = expression('Fish Density, L' ^ -1),  ))

J = ggplot(data=CopepodCountsSummaryJOA, aes(x=Week, y= mean, colour=FishPerL), inherit.aes=F, show.legend = FALSE) +
  labs(y = expression('JOA Density, L' ^ -1)) + geom_line(linewidth=1.5) + geom_point(size = 3) + 
  geom_errorbar(data=CopepodCountsSummaryJOA, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=FishPerL), width=0.2, inherit.aes=FALSE, show.legend = TRUE) +
  theme_classic() + guides(color = guide_legend(title = expression('Fish Density, L' ^ -1),  ))

A = ggplot(data=CopepodCountsSummaryAF, aes(x=Week, y= mean, colour=FishPerL), inherit.aes=F, show.legend = FALSE) +
  labs(y = expression('AF Copepod Density, L' ^ -1)) + geom_line(linewidth=1.5) + geom_point(size = 3) + 
  geom_errorbar(data=CopepodCountsSummaryAF, aes(x=Week, ymin=mean-SE, ymax=mean+SE, colour=FishPerL), width=0.2, inherit.aes=FALSE, show.legend = TRUE) +
  theme_classic() + guides(color = guide_legend(title = expression('Fish Density, L' ^ -1),  ))

ggarrange(total,N,J,A,
          nrow = 2, ncol = 2,
          labels = c("A","B","C","D"),
          common.legend = TRUE,
          legend = "right")


