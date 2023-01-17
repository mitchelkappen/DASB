######## New Analysis DASB Paper ######## 
# Written by Mitchel Kappen #
# Date: 22-12-2021 #
######## Package getting ##########
rm(list = ls()) # Clear environment
cat("\014") # Clear console
dev.off() # Clear plot window

library(yarrr)
library(lme4)
library(emmeans)
library(pander)

library(reshape)
library(lme4)
library(lmerTest)
library(pander)

library(ggpubr)
library(car)

library(arrow)
library(tibble)
library(effects)


########## Settings and data getting ###########
nAGQ = 1 # Set to 1 for eventual analysis

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

data <- read.csv(paste0("dataDASB.csv")) # Load data

names(data)[names(data) == "Dier.nr"] <- "Diernr" # Rename this column because working with dots is annoying

######## DATA preparation #########
# Factorize the right variables
data$Diernr <- as.factor(data$Diernr)
data$Dier <- as.factor(data$Dier)
data$Treatment <- as.factor(data$Treatment)
data$Time <- as.factor(data$Time)
data$time <- as.factor(data$time)

## Overview different groups ##
dogs5active <- unique(data$Dier[data$Treatment == "5 sessions active"])
dogs20active <- unique(data$Dier[data$Treatment == "20 sessions active"])
data <- data[data$Treatment!= "20 sessions sham",] # Taking out sham data - not using that here

# Define a function to calculate overlaps
iselement <- function(x, A) {
  if (x %in% A) {
    return(TRUE)
  }
  return(FALSE)
}

# 5 vs 20 active
count = 0
for (i in 1:length(dogs20active)) {
  if (iselement(dogs20active[i], dogs5active) == TRUE){
    count = count + 1
  }
}
print(paste0("Out of all ", length(dogs20active), " dogs in group 20 sessions active, ", count, " were also in group 5 sessions active"))

########### DATA analysis ########
# Time effect, x group
data$Time2[data$Time == "Baseline"] = "1) baseline"
data$Time2[data$Time == "24h "] = "2) 24h"
data$Time2[data$Time == "1m"] = "3) 1m"
data$Time2[data$Time == "3m"] = "4) 3m "

# Define formula
formula <- subgen.bis ~ Time2 * Treatment + (1|Diernr)

# Model fitting
d0.1 <- lmer(formula,data=data)

d0.2 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 10000000)),nAGQ = nAGQ)

d0.3 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

modelNames = c(d0.1, d0.2, d0.3)

# Model Selection
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III') # Run Anova, double square brackets because of list properties

# Main effect: Time
emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ Time2, adjust ="fdr", type = "response")
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
plot(effect("Time2", chosenModel[[1]])) #just to check

# Interaction effect
emmeans0.2 <- emmeans(chosenModel[[1]], pairwise ~ Treatment | Time2, adjust ="fdr", type = "response")
emm0.2 <- summary(emmeans0.2)$emmeans
emmeans0.2$contrasts
plot(effect("Time2:Treatment ", chosenModel[[1]])) #just to check

## Visualisation
tiff("../figures/figure2.tiff", units="in", width=10, height=5, res=300) # Save to tiff
pd <- position_dodge(0.01) # Plotting setting: move them .05 to the left and right
# Main effect: Time
ggplot(emm0.1, aes(x=Time2, y=emmean)) +
  geom_point(size = 1) + 
  geom_line(aes(group = 1),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  # ggtitle("Subgenual ACC at different time points over groups")+
  labs(y = "sgACC", x = "")+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title = element_text(size = 20))  

dev.off() # Close jpeg file and save it

# Interaction effect
pd <- position_dodge(0.03) # Plotting setting: move them .05 to the left and right
tiff("../figures/figure3.tiff", units="in", width=10, height=5, res=300) # Save to tiff
ggplot(emm0.2, aes(x=Time2, y=emmean, color=Treatment)) +
  geom_point(size = 1) + 
  geom_line(aes(group = Treatment),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  # ggtitle("Treatment")+
  labs(y = "sgACC", x = "")+
  theme(axis.text.x = element_text(size = 16))+ # X Axis ticks
  theme(axis.text.y = element_text(size = 10))+ # Y axis ticks
  theme(axis.title = element_text(size = 16))+ # Axis titles
  theme(legend.text = element_text(size = 16))+ # Legend text
  theme(legend.title = element_text(size = 14)) # Legend title
dev.off() # Close tiff file and save it
