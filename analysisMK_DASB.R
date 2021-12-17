# New Analysis DASB Paper

######### TODO ############
# Yangfeng:
#   
#   Make a table of how balanced the dogs are in every group
# 
# To test:
#   
#   Compare between groups: sham vs active at 20 sessions
# active 20 vs 5 sessions on group levels
# 
# compare over time effect, x group
# 
# Pons Left thalamus Presubgenual cortex Subgenual cortex
# 
# page 125
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
plotPrefix = "Plots/"

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

# After loading Robrecht .RData
# write.csv(data_rTMS_anes, paste0("MKData/dataDASB.csv"))

data <- read.csv(paste0("MKData/dataDASB.csv"))

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
# dogs5sham <- unique(data$Dier[data$Treatment == "5 sessions sham"])
dogs20active <- unique(data$Dier[data$Treatment == "20 sessions active"])
dogs20sham <- unique(data$Dier[data$Treatment == "20 sessions sham"])

# Define a function to calculate overlaps
iselement <- function(x, A) {
  if (x %in% A) {
    return(TRUE)
  }
  return(FALSE)
}

# 20 vs 20
count = 0
for (i in 1:length(dogs20active)) {
  if (iselement(dogs20active[i], dogs20sham) == TRUE){
    count = count + 1
  }
}
print(paste0("Out of all ", length(dogs20active), " dogs in group 20 sessions active, ", count, " were also in group 20 sessions sham"))

# 5 vs 20 active
count = 0
for (i in 1:length(dogs20active)) {
  if (iselement(dogs20active[i], dogs5active) == TRUE){
    count = count + 1
  }
}
print(paste0("Out of all ", length(dogs20active), " dogs in group 20 sessions active, ", count, " were also in group 5 sessions active"))
print("This is a very big overlap so we will not be able to compare these.. ")

########### Chapter 1: full analysis ########
## First do all analysis?? ##
dataBackup <- data # Here the full dataframe is stored
## 20 active vs 20 sham ##
# data <- data[data$Treatment!= "5 sessions active",]
data <- data[data$Treatment!= "20 sessions sham",]

# time effect, x group
data$Time2[data$Time == "Baseline"] = "1) baseline"
data$Time2[data$Time == "24h "] = "2) 24h"
data$Time2[data$Time == "1m"] = "3) 1m"
data$Time2[data$Time == "3m"] = "4) 3m "

# testing = 'Main effect' # Else set to Interaction effect (or Main effect)
formula <- subgen.bis ~ Time2 * Treatment + (1|Diernr)
# formula <- presubgen.bis ~ Time2 * Treatment + (1|Diernr)
# formula <- pons ~ Time2 * Treatment + (1|Diernr)
# formula <- hippocampus.L ~ Time2 * Treatment + (1|Diernr)


# Model
d0.1 <- lmer(formula,data=data)
d0.2 <- glmer(formula,data=data, family = gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.3 <- glmer(formula,data=data, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.5 <- glmer(formula,data=data, family = Gamma(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.6 <- glmer(formula,data=data, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
  
d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.8 <- glmer(formula,data=data, family = inverse.gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.9 <- glmer(formula,data=data, family = inverse.gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
 
d0.1 <- lmer(formula,data=data)
d0.2 <- glmer(formula,data=data, family = gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))
d0.3 <- glmer(formula,data=data, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))

d0.4 <- glmer(formula,data=data, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))
d0.5 <- glmer(formula,data=data, family = Gamma(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))
d0.6 <- glmer(formula,data=data, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))

d0.7 <- glmer(formula,data=data, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))
d0.8 <- glmer(formula,data=data, family = inverse.gaussian(link = "inverse"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))
d0.9 <- glmer(formula,data=data, family = inverse.gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)))


 
modelNames = c(d0.1)
# modelNames = c(d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9)

# Model Selection
# tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3), AIC(d0.4), AIC(d0.5), AIC(d0.6), AIC(d0.7), AIC(d0.8), AIC(d0.9))
tabel <- cbind(AIC(d0.1))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
# chosenModel = d0.1

Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ Treatment | Time2, adjust ="fdr", type = "response")
emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ Treatment | Time2, type = "response")

# emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ Time2 | Treatment, adjust ="fdr", type = "response")
# emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ Time2 * Treatment, type = "response")

emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts
plot(effect("Time2:Treatment ", chosenModel[[1]])) #just to check

pd <- position_dodge(0.01) # move them .05 to the left and right

## Visualisation
ggplot(emm0.1, aes(x=Time2, y=emmean, color=Treatment)) +
  geom_point(size = 1) + 
  geom_line(aes(group = Treatment),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  # geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Treatment")+
  labs(y = "Subgenial")
  # annotate(geom="text", x=2, y=1, label="*", color="#000000")+ #24h
  # annotate(geom="text", x=3, y=.84, label="**", color="#000000") #1m

# Reponse???
ggplot(emm0.1, aes(x=Time2, y=response, color=Treatment)) +
  geom_point(size = 1) + 
  geom_line(aes(group = Treatment),size = 1)+
  geom_errorbar(width=.125, aes(ymin=response-SE, ymax=response+SE), position=pd)+
  # geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Treatment")+
  labs(y = "Subgenial")


########### Chapter 2: Baseline correction makes more sense #######
require(data.table)
data <- dataBackup
data <- data[data$Treatment!= "20 sessions sham",]

data <- setDT(data) # Turn into data.frame
# data[, delta.subgen_bis:=c(NA, subgen.bis[-.N]), by=Diernr]
data[ , delta.subgen_bis := subgen.bis - shift(subgen.bis), by = Diernr] # Create delta scores
data <- data[is.na(data$delta.subgen_bis) == FALSE, ] # Get rid of created NA's (becaus baselines dont have delta's)
data$deltaTimes[data$Time == "24h "] = "1) Delta 24h - baseline"
data$deltaTimes[data$Time == "1m"] = "2) Delta 1m - 24h"
data$deltaTimes[data$Time == "3m"] = "3) Delta 3m - 1m"

# Define formula
formula <- delta.subgen_bis ~ Treatment * deltaTimes + (1|Diernr)

data$deltaTimes <- as.ordered(data$deltaTimes)

# Model generation
d1.1 <- lmer(formula,data=data)

# Stat eval
emmeans1.1 <- emmeans(d1.1, pairwise ~ Treatment | deltaTimes, adjust ="fdr", type = "response") # Compute a variable containing all emmeans/contrasts
# emmeans1.1 <- emmeans(d1.1, pairwise ~ deltaTimes| Treatment, adjust ="fdr", type = "response") # Compute a variable containing all emmeans/contrasts

emm1.1 <- summary(emmeans1.1)$emmeans
emmeans1.1$contrasts

## Visualisation
ggplot(emm1.1, aes(x=deltaTimes, y=emmean, color=Treatment)) +
  geom_point(size = 1) + 
  geom_line(aes(group = Treatment),size = 1)+
  geom_errorbar(width=.125, aes(ymin=emmean-SE, ymax=emmean+SE), position=pd)+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw(base_size = 8)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  ggtitle("Treatment")+
  labs(y = "Delta subgenial")
  # annotate(geom="text", x=xplotPosition, y=-37.5, label="**", color="#000000")+ #IBI3
  # annotate(geom="text", x=xplotPosition + 1, y=-34.5, label="**", color="#000000")+ #IBI4
  # annotate(geom="text", x=xplotPosition + 2, y=-26, label="***", color="#000000")+ #IBI5
  # annotate(geom="text", x=xplotPosition + 3, y=-17, label="***", color="#000000")+ #IBI6
  # annotate(geom="text", x=xplotPosition + 4, y=-12.5, label="**", color="#000000") #IBI7




