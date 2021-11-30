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

nAGQ = 0 # Set to 1 for eventual analysis
plotPrefix = "Plots/"

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location

# After loading Robrecht .RData
# write.csv(data_rTMS_anes, paste0("MKData/dataDASB.csv"))

data <- read.csv(paste0("MKData/dataDASB.csv"))

names(data)[names(data) == "Dier.nr"] <- "Diernr" # Rename this column because working with dots is annoying

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


## First do all analysis?? ##
dataBackup <- data # Here the full dataframe is stored
## 20 active vs 20 sham ##
# data <- data[data$Treatment!= "5 sessions active",]
data <- data[data$Treatment!= "20 sessions sham",]

# time effect, x group

# testing = 'Main effect' # Else set to Interaction effect (or Main effect)
testing = 'Interaction effect'

if (testing == 'Main effect'){
  formulas = c('pons ~ Treatment', 'Thalamus.L ~ Treatment', 'presubgen.bis ~ Treatment', 'subgen.bis ~ Treatment')
} else if (testing == 'Interaction effect') {
  formulas = c('pons ~ time * Treatment', 'Thalamus.L ~ time * Treatment', 'presubgen.bis ~ time * Treatment', 'subgen.bis ~ time * Treatment')
}

plotTitles = c('Pons', 'Thalamus.L', 'presubgen.bis', 'subgen.bis')

for(i in 4) {
  formula <- paste0(formulas[i], ' + (1|Diernr)')
  # formula <- paste0(formulas[i], ' + (1|Dier)')
  
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
  
  modelNames = c(d0.1,d0.2,d0.3,d0.4,d0.5,d0.6,d0.7,d0.8,d0.9)
  
  # Model Selection
  tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.3), AIC(d0.4), AIC(d0.5), AIC(d0.6), AIC(d0.7), AIC(d0.8), AIC(d0.9))
  chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC
  
  Anova(chosenModel[[1]]) # Run Anova, double square brackets because of list properties
  print("Stats baseline vs control vs stress:")
  
  if (testing == 'Main effect'){
    print(emmeans(chosenModel[[1]], pairwise ~ Treatment , adjust ="fdr", type="response")) # This is the right one for self-reports
  } else if (testing == 'Interaction effect') {
    print(emmeans(chosenModel[[1]], pairwise ~ time * Treatment , adjust ="fdr", type="response")) # This is the right one for self-reports
  }
  
  # Plotting
  dpi=600    #pixels per square inch
  # jpeg(paste0(plotPrefix, "Figure", "_", plotTitles[i], ".jpeg"), width=8*dpi, height=4*dpi, res=dpi)
  par(mfcol = c(1, 1))
  plotAROUSAL <- pirateplot(
    formula = formulas[i],
    data = data,
    theme = 1,
    pal = "info",
    main = plotTitles[i],
    bean.f.o = .6, # Bean fill
    point.o = .3,  # Points
    inf.f.o = .7,  # Inference fill
    inf.b.o = .8,  # Inference border
    avg.line.o = 1,  # Average line
    # bar.f.o = .5, # Bar
    inf.f.col = "white",  # Inf fill col
    inf.b.col = "black",  # Inf border col
    avg.line.col = "black",  # avg line col
    bar.f.col = gray(.8),  # bar filling color
    point.pch = 21,
    point.bg = "white",
    point.col = "black",
    point.cex = .7,
    
    xlab = "",
  )
  # abline(lm(formulas[i], data=data), lwd=4, lty=2, col = "red")
  
  # Secondary points and axis:
  # mtext("fileNum",1,line=1,at=0.2,col="red")
  
  # mtext("Experiment phase",1,line=3,at=0.2,col="blue")
  # axis(side=1, at=c(1:9), line = 3, labels=c('prerest','control1','control2','control3','midrest','stress1','stress2','stress3','postrest' ))
  # axis(side=1, at=c(1:6), line = 3, labels=c('control1','control2','control3','stress1','stress2','stress3' ))
  
  # dev.off()
}
