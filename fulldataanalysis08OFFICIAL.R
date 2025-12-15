
library(here)
library(lme4)
library(ggplot2)
library(DHARMa)
library(gridExtra)
library(ggeffects)
library(tidyverse)
library(nlme)
library(predictmeans)
library(sjPlot)
library(glmmTMB)
library(car)
library(emmeans)
library(MASS)
library(sjlabelled)
library(sjmisc)
library(ggrepel)
library(maps)
library(xtable)
library(vegan)
library(carData)
library(car)
library(performance)
library(MuMIn)
library(glmm.hp)
library(eulerr)
library(ggeffects)
library(vegan)
library(reshape2)
library(pheatmap)
library(dichromat)

dfmalesClimate <- readRDS(file = "dfmalesClimate.rds")
dffemalesClimate <- readRDS(file = "dffemalesClimate.rds")
dfmalesphenotypingClimate <- readRDS(file = "dfmalesphenotypingClimate.rds")
dffemalesphenotypingClimate <- readRDS(file = "dffemalesphenotypingClimate.rds")
sumdfClimate <- readRDS(file = "sumdfClimate.rds")
dfSClimate <- readRDS(file = "dfSClimate.rds")
dfClimate <- readRDS(file = "dfClimate.rds")
dfmales <- dfmalesClimate  %>% drop_na(biomass) 
dffemales <- dffemalesClimate  %>% drop_na(biomass) 
dfmales <- dfmalesClimate  %>% drop_na(Height..cm.) 
dffemales <- dffemalesClimate  %>% drop_na(Height..cm.) 
dfSClimate <- dfSClimate  %>% drop_na(weight1seed) 
dffemales$biomassSC <- sqrt(dffemales$biomass)
dfmales$biomassSC <- sqrt(dfmales$biomass)

######### MAIN TEXT
### FULL MODELS SEP SEX AIR

glmmIsLeakyM <- glmmTMB(is.leaky  ~ Air * Hormone.Tag  + 
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ Air * Hormone.Tag + 
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmLDM <- glmmTMB(.Female.flowers.or.fruits ~ Air * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dfmalesClimate, 
                   family = nbinom2, ziformula = ~ Air)
glmmLDF <- glmmTMB(.Male.flowers  ~ Air * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dffemalesClimate, 
                   family = nbinom2, ziformula = ~ Air)
glmmREM <- glmmTMB(.Male.flowers  ~  Air * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  Air * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)
glmmBiomassF <- glmmTMB(biomassSC  ~ Air * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
glmmBiomassM <- glmmTMB(biomassSC  ~ Air * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightM <- glmmTMB(Height..cm.  ~ Air * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightF <- glmmTMB(Height..cm.  ~ Air * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
sweight <- glmmTMB(weight1seed  ~ Air, 
                   data = dfSClimate)

sarea <- glmmTMB(Area  ~ 
                   Air +
                   (1|ID)
                 , 
                 data = dfClimate)

glmmBiomassF <- glmmTMB(biomassSC  ~ Air * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
glmmBiomassF <- glmmTMB(biomassSC  ~ Air * Hormone.Tag + (1|Population), data = dffemales)
#glmmBiomassF <- glmmTMB(biomassSC  ~ Air * Hormone.Tag + (1|Table), data = dffemales)

summary(glmmIsLeakyM)
summary(glmmIsLeakyF)
summary(glmmLDM)
summary(glmmLDF)
summary(glmmREM)
summary(glmmREF)
summary(glmmBiomassM)
summary(glmmBiomassF)
summary(glmmHeightM)
summary(glmmHeightF)
summary(sweight)
summary(sarea)

darmaglmmRS=simulateResiduals(fittedModel = glmmIsLeakyM,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmIsLeakyF,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmLDM,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmLDF,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmREM,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmREF,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmBiomassM,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmBiomassF,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmHeightM,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = glmmHeightF,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = sweight,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)

darmaglmmRS=simulateResiduals(fittedModel = sarea,n=5000)
plot(darmaglmmRS) 
testDispersion(darmaglmmRS)
testZeroInflation(darmaglmmRS)


### FULL MODELS SEP SEX AIR PCs

glmmIsLeakyM <- glmmTMB(is.leaky  ~ Air  + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ Air  + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmREM <- glmmTMB(.Male.flowers  ~  Air + PC1 + PC2 +
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  Air + PC1 + PC2 +
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)

summary(glmmIsLeakyM)
summary(glmmIsLeakyF)
summary(glmmREM)
summary(glmmREF)

check_collinearity(glmmIsLeakyM)
check_collinearity(glmmIsLeakyF)
check_collinearity(glmmREM)
check_collinearity(glmmREF)


varpartitionwohormoneAIR <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Geo Dist.&PC1" = round(max(0,ab), 3),
    "Geo Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Geo Dist.&PC1&PC2" = round(max(0,abc), 3)
  )
  # Generate and plot the diagram
  fit <- euler(venn_data, shape = "ellipse")
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("salmon1", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, # list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionlmAIR <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Geo Dist.&PC1" = round(max(0,ab), 3),
    "Geo Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Geo Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- euler(venn_data, shape = "ellipse")
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("#E31A1C", "#1F78B4", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionwohormoneAIRnum <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  # Generate and plot the diagram
  x <- c(residua,ua,ub,uc,bc, 0)
  return(x)
} 
varpartitionlmAIRnum <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  x <- c(residua,ua,ub ,uc,bc, 0)
  return(x)
} 

varpartitionwohormoneAIR <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Geo Dist.&PC1" = round(max(0,ab), 3),
    "Geo Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Geo Dist.&PC1&PC2" = round(max(0,abc), 3)
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("indianred1", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, # list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionlmAIR <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Air-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Geo Dist.&PC1" = round(max(0,ab), 3),
    "Geo Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Geo Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("indianred1", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 

varpartitionwohormoneAIR(glmmIsLeakyF)
varpartitionwohormoneAIR(glmmIsLeakyM)
varpartitionwohormoneAIR(glmmREF)
varpartitionwohormoneAIR(glmmREM)

LFA <- varpartitionwohormoneAIRnum(glmmIsLeakyF)
LMA <- varpartitionwohormoneAIRnum(glmmIsLeakyM)
RFA <- varpartitionwohormoneAIRnum(glmmREF)
RMA <- varpartitionwohormoneAIRnum(glmmREM)




a <- 2000
b <- 20
c <- 10
newdat <- data.frame(x = x, y =  1/(1+exp(-(-4.07464+2.35702*x))))
x = seq(min(dfSClimate$Air), max(dfSClimate$Air), 
        by = ((max(dfSClimate$Air) - min(dfSClimate$Air))/
                a))
newdatF <- data.frame(x = x, y =  plogis((-3.4625869+0.0008605*x)))
newdatM <- data.frame(x = x, y =  plogis((-5.9861965+0.0011948*x)))

ggplot() +
  geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                     y = leakymalescountwhormone / (sumMwH )), color = "#223F47", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                     y = leakyfemalescountwhormone /(sumFwH)), color = "#A56B1F", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                     y = leakymalescountnohormone / (sumMnoH)), color = "#5A99A5", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                     y = leakyfemalescountnohormone /(sumFnoH)), color = "#BFA144", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data=newdatM, aes(x = x, y = y), size = 0.5, color = "#376270", shape = 15) +
  geom_point(data=newdatF, aes(x = x, y = y), size = 0.5, color = "goldenrod1", shape = 15) +
  xlab(label = 'Geodesic distance to core (km)') +
  ylab(label = 'Incidence of leakiness') + theme_bw() 

  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())

newdatM <- data.frame(x = x, y =  exp(-1.792e+00 + 1.333e-04 *x))
newdatF <- data.frame(x = x, y =  exp(-4.190e+00 + 3.475e-04*x))

ggplot() +
  geom_point(data= dffemalesphenotypingClimate[dffemalesphenotypingClimate$Hormone.Tag =="x",], 
             aes(x = Air, y = .Female.flowers.or.fruits/ biomass), color = "#A56B1F", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data= dffemalesphenotypingClimate[dffemalesphenotypingClimate$Hormone.Tag =="-",], 
             aes(x = Air, y = .Female.flowers.or.fruits/ biomass), color = "#BFA144", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data=newdatF, aes(x = x, y = y), size = 0.5, color = "goldenrod1", shape = 15) +
  xlab(label = 'Geodesic distance to core (km)') +
  ylab(label = 'Female reproductive effort') + theme_bw() 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())

ggplot() +
  geom_point(data= dfmalesphenotypingClimate[dfmalesphenotypingClimate$Hormone.Tag =="x",], 
             aes(x = Air,y = .Male.flowers / biomass), color = "#223F47", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data= dfmalesphenotypingClimate[dfmalesphenotypingClimate$Hormone.Tag =="-",], 
             aes(x = Air, y = .Male.flowers/ biomass), color = "#5A99A5", size = 2, alpha = 0.8, width = 0.5, height = 0) +
  geom_point(data=newdatM, aes(x = x, y = y), size = 0.5, color = "#376270", shape = 15) +
  xlab(label = 'Geodesic distance to core (km)') +
  ylab(label = 'Male reproductive effort') + theme_bw() 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())
  

  ggplot() +
    geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                       y = leakymalescountwhormone / (sumMwH )), color = "#223F47", size = 0.5, alpha = 1) +
    geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                       y = leakyfemalescountwhormone /(sumFwH)), color = "#A56B1F", size = 0.5, alpha = 1) +
    geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                       y = leakymalescountnohormone / (sumMnoH)), color = "#5A99A5", size = 0.5, alpha = 1) +
    geom_point(data= sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                                       y = leakyfemalescountnohormone /(sumFnoH)), color = "#BFA144", size = 0.5, alpha = 1) +
    geom_function(fun = \(x) plogis((-3.4625869+0.0008605*x)), color = "darkgoldenrod4")+
    geom_function(fun = \(x) plogis((-5.9861965+0.0011948*x)), color = "#376270")+
    xlab(label = 'Geodesic distance to core (km)') +
    ylab(label = 'Incidence of leakiness') + theme_bw() 
  
  ggplot() +
    geom_point(data= dffemalesphenotypingClimate[dffemalesphenotypingClimate$Hormone.Tag =="x",], 
               aes(x = Air, y = .Female.flowers.or.fruits/ biomass), color = "#A56B1F", size = 0.5, alpha = 1) +
    geom_point(data= dffemalesphenotypingClimate[dffemalesphenotypingClimate$Hormone.Tag =="-",], 
               aes(x = Air, y = .Female.flowers.or.fruits/ biomass), color = "#BFA144", size = 0.5, alpha = 1) +
    geom_function(fun = \(x) exp(-4.190e+00 + 3.475e-04*x), color = "darkgoldenrod4") +
    xlab(label = 'Geodesic distance to core (km)') +
    ylab(label = 'Female reproductive effort') + theme_bw() 
  
  ggplot() +
    geom_point(data= dfmalesphenotypingClimate[dfmalesphenotypingClimate$Hormone.Tag =="x",], 
               aes(x = Air,y = .Male.flowers / biomass), color = "#223F47", size = 0.5, alpha = 1) +
    geom_point(data= dfmalesphenotypingClimate[dfmalesphenotypingClimate$Hormone.Tag =="-",], 
               aes(x = Air, y = .Male.flowers/ biomass), color = "#5A99A5", size = 0.5, alpha = 1 )+
    geom_function(fun = \(x) exp(-1.792e+00 + 1.333e-04 *x), color = "#376270") +
    xlab(label = 'Geodesic distance to core (km)') +
    ylab(label = 'Male reproductive effort') + theme_bw() 


###CORRELATION 
ggplot(data = sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular)) +
  geom_point(aes(y = PC1), color = "black", size = 2, shape = 16) +
  geom_point(aes(y = PC2), color = "black", size = 2, shape = 1)  +
  xlab(label = 'Geodesic distance to core (km)')+
         ylab(label = 'Climate principal component scores') + theme_bw()

forcor <- data.frame(cbind(sumdfClimate$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                           sumdfClimate$PC1,
                           sumdfClimate$PC2))
colnames(forcor) <- c("GeoD", "PC1", "PC2")
forcor2 <- data.frame(round(cor(forcor),3)) 
pheatmap(forcor2, main = "", display_numbers = T, 
         cluster_rows = F, cluster_cols = F, fontsize_number = 12 ,
         fontface_number = 1 ,number_color = "#000000", 
         color = colorRampPalette(c("cornflowerblue", "white", "salmon2"))(50),
         show_rownames = T, na_col = "white")

varfull <- data.frame(LFA,LMA, RFA, RMA )
colnames(varfull) <- c("Leakiness probablity F","Leakiness probablity M", "Reproductive effort F" ,"Reproductive effort M"  
)
varfull[6,] <- 1- colSums(varfull)
rownames(varfull) <- c("residuals","ua","ub", "uc","bc", "Others_joint")
varfull[varfull < 0] <- 0
varfull$ID <- c("residuals","Geo Dist.", "PC1","PC2", "PC joint", "Others joint")
varfull2 <- varfull[-1,]
varfull2freq <- apply(varfull2[,-5],2,function(x){x/sum(x)})
varfull2freq <- data.frame(varfull2freq)
varfull2freq$ID <- c("Geo Dist." , "PC1","PC2", "PC joint", "Others joint")
varfullstack <- melt(varfull2freq, id.vars=c('ID'), variable.name='Response_Trait', value.name='Explained_Marginal_Variance')
ggplot(varfullstack, aes(fill=factor(ID, levels=c("Others joint" ,
                                                  "PC joint",
                                                  "PC2",
                                                  "PC1",
                                                  "Geo Dist.")), y=Explained_Marginal_Variance, x=Response_Trait)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#999999", "#33A1C9", "#A6CEE3", "cadetblue3", "indianred1")) +
  labs(fill = "Predictors") + theme_bw() 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())



######### APPENDIX A2

### FULL MODELS SEP SEX WALK

glmmIsLeakyM <- glmmTMB(is.leaky  ~ Walk * Hormone.Tag  + 
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ Walk * Hormone.Tag + 
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmLDM <- glmmTMB(.Female.flowers.or.fruits ~ Walk * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dfmalesClimate, 
                   family = nbinom2, ziformula = ~ Walk)
glmmLDF <- glmmTMB(.Male.flowers  ~ Walk * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dffemalesClimate, 
                   family = nbinom2, ziformula = ~ Walk)
glmmREM <- glmmTMB(.Male.flowers  ~  Walk * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  Walk * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)
glmmBiomassF <- glmmTMB(biomassSC  ~ Walk * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
glmmBiomassM <- glmmTMB(biomassSC  ~ Walk * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightM <- glmmTMB(Height..cm.  ~ Walk * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightF <- glmmTMB(Height..cm.  ~ Walk * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
sweight <- glmmTMB(weight1seed  ~ Walk, 
                   data = dfSClimate)
sarea <- glmmTMB(Area  ~ 
                   Walk +
                   (1|ID)
                 , 
                 data = dfClimate)

ggplot() + geom_point(data= dfSClimate, 
                      aes(x = Walk , y = weight1seed * 1000), color = "#BFA144", size = 0.7, alpha = 1) +
  geom_function(fun = \(x) 1000*(1.998e-03 - 1.086e-07*x), color = "darkgoldenrod4")+
  xlab(label = 'Land distance from core (km)') +
  ylab(label = 'Seed weight (mg)') + theme_bw()

ggplot() +
  geom_point(data= sumdfClimate, aes(x = WalktoHaifavisGGmap, 
                                     y = leakymalescountwhormone / (sumMwH )), color = "#223F47", size = 0.7, alpha = 1) +
  geom_point(data= sumdfClimate, aes(x = WalktoHaifavisGGmap, 
                                     y = leakyfemalescountwhormone /(sumFwH)), color = "#A56B1F", size = 0.7, alpha = 1) +
  geom_point(data= sumdfClimate, aes(x = WalktoHaifavisGGmap, 
                                     y = leakymalescountnohormone / (sumMnoH)), color = "#5A99A5", size = 0.7, alpha = 1) +
  geom_point(data= sumdfClimate, aes(x = WalktoHaifavisGGmap, 
                                     y = leakyfemalescountnohormone /(sumFnoH)), color = "#BFA144", size = 0.7, alpha = 1) +
  geom_function(fun = \(x) plogis((-3.228e+00 +5.136e-04 *x)), color = "darkgoldenrod4")+
  geom_function(fun = \(x) plogis((-7.0331717 +0.0010635*x)), color = "#376270")+
  xlab(label = 'Land distance from core (km)') +
  ylab(label = 'Incidence of leakiness') + theme_bw() 

ggplot(data = sumdfClimate, aes(x = AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular)) +
  geom_point(aes(y = WalktoHaifavisGGmap), color = "black", size = 1)  +
  xlab(label = 'Geodesic distance to core (km)') +
  ylab(label = 'Land distance from core (km)') + theme_bw() 

forcor <- data.frame(cbind(sumdfClimate$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular, 
                           sumdfClimate$WalktoHaifavisGGmap,
                           sumdfClimate$PC1,
                           sumdfClimate$PC2))
colnames(forcor) <- c("GeoD","LandD", "PC1", "PC2")
forcor2 <- data.frame(round(cor(forcor),3)) 
pheatmap(forcor2, main = "", display_numbers = T, 
         cluster_rows = F, cluster_cols = F, fontsize_number = 12 ,
         fontface_number = 1 ,number_color = "#000000", 
         color = colorRampPalette(c("cornflowerblue", "white", "salmon2"))(50),
         show_rownames = T, na_col = "white")

glmmIsLeakyM <- glmmTMB(is.leaky  ~ Walk  + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ Walk  + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
sweight <- glmmTMB(weight1seed  ~ Walk + PC1 + PC2, 
                   data = dfSClimate)

summary(glmmIsLeakyM)
summary(glmmIsLeakyF)
summary(sweight)
check_collinearity(glmmIsLeakyM)
check_collinearity(glmmIsLeakyF)
check_collinearity(sweight)

varpartitionwohormoneWALK <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Land Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Land Dist.&PC1" = round(max(0,ab), 3),
    "Land Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Land Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- euler(venn_data, shape = "ellipse")
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("#FDBF6F", "#1F78B4", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionlmWALK <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Land Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Land Dist.&PC1" = round(max(0,ab), 3),
    "Land Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Land Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- euler(venn_data, shape = "ellipse")
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("#FDBF6F", "#1F78B4", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionwohormoneWALKnum <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  # Generate and plot the diagram
  x <- c(residua,ua,ub,uc,bc, 0)
  return(x)
} 
varpartitionlmWALKnum <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  x <- c(residua,ua,ub ,uc,bc, 0)
  return(x)
} 

varpartitionwohormoneWALK <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Land Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Land Dist.&PC1" = round(max(0,ab), 3),
    "Land Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Land Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("#FDBF6F", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionlmWALK <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . -PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I3mod <- update(fullmod, formula=. ~ . -Walk-PC1)
  
  I12mod <- update(fullmod, formula=. ~ . -PC2)
  I13mod <- update(fullmod, formula=. ~ . -PC1)
  I23mod <- update(fullmod, formula=. ~ . - Walk)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted
  
  RI12mod <- r2(I12mod)$R2_adjusted
  RI13mod <- r2(I13mod)$R2_adjusted
  RI23mod <- r2(I23mod)$R2_adjusted
  
  # Individual sets
  ua <- Rfullmod - RI23mod
  ub <- Rfullmod - RI13mod
  uc <- Rfullmod - RI12mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI3mod - ua - ub
  ac <- Rfullmod - RI2mod - ua - uc
  bc <- Rfullmod - RI1mod - ub - uc
  
  # Triplet intersections
  abc <- Rfullmod - ua - ub - uc- ab-ac-bc
  
  va <- sum(ua+ub+uc+
              ab+ac+bc+
              abc)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Land Dist." = round(max(0,ua), 3),
    "PC1" = round(max(0,ub), 3),
    "PC2" = round(max(0,uc), 3),
    
    # Pairwise intersections
    "Land Dist.&PC1" = round(max(0,ab), 3),
    "Land Dist.&PC2" = round(max(0,ac), 3),
    "PC1&PC2" = round(max(0,bc), 3),
    
    # Triplet intersections
    "Land Dist.&PC1&PC2" = round(max(0,abc), 3)
    
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("#FDBF6F", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 

varpartitionwohormoneWALK(glmmIsLeakyF)
varpartitionwohormoneWALK(glmmIsLeakyM)
varpartitionlmWALK(sweight)

LFW <- varpartitionwohormoneWALKnum(glmmIsLeakyF)
LMW <- varpartitionwohormoneWALKnum(glmmIsLeakyM)
SWW <- varpartitionlmWALKnum(sweight)

varfull <- data.frame(LFW,LMW, SWW )
colnames(varfull) <- c("Leakiness probablity F","Leakiness probablity M", "Seed weight"  
)
varfull[6,] <- 1- colSums(varfull)
rownames(varfull) <- c("residuals","ua","ub", "uc","bc", "Others_joint")
varfull[varfull < 0] <- 0
varfull$ID <- c("residuals","Land Dist.", "PC1","PC2", "PC joint", "Others joint")
varfull2 <- varfull[-1,]
varfull2freq <- apply(varfull2[,-4],2,function(x){x/sum(x)})
varfull2freq <- data.frame(varfull2freq)
varfull2freq$ID <- c("Land Dist." , "PC1","PC2", "PC joint", "Others joint")
varfullstack <- melt(varfull2freq, id.vars=c('ID'), variable.name='Response_Trait', value.name='Explained_Marginal_Variance')
ggplot(varfullstack, aes(fill=factor(ID, levels=c("Others joint" ,
                                                  "PC joint",
                                                  "PC2",
                                                  "PC1",
                                                  "Land Dist.")), y=Explained_Marginal_Variance, x=Response_Trait)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#999999", "#33A1C9", "#A6CEE3", "cadetblue3", "#FDBF6F")) +
  labs(fill = "Predictors") + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())

######### APPENDIX A3

### FULL MODELS SEP SEX AIR WALK PCs

glmmIsLeakyM <- glmmTMB(is.leaky  ~ Air  + Walk + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ Air  + Walk + PC1 + PC2 +
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmREM <- glmmTMB(.Male.flowers  ~  Air + Walk + PC1 + PC2 +
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  Air + Walk + PC1 + PC2 +
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)
sweight <- glmmTMB(weight1seed  ~ Air  + Walk + PC1 + PC2, 
                   data = dfSClimate)

varpartitionwohormone <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . - Walk-PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC1-PC2)
  I3mod <- update(fullmod, formula=. ~ . - Walk-Air-PC2)
  I4mod <- update(fullmod, formula=. ~ . - Walk-PC1-Air)
  
  I12mod <- update(fullmod, formula=. ~ . - PC1-PC2)
  I13mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I14mod <- update(fullmod, formula=. ~ . - Walk-PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I24mod <- update(fullmod, formula=. ~ . - Air-PC1)
  I34mod <- update(fullmod, formula=. ~ . - Walk-Air)
  
  I123mod <- update(fullmod, formula=. ~ . - PC2)
  I124mod <- update(fullmod, formula=. ~ . - PC1)
  I134mod <- update(fullmod, formula=. ~ . - Walk)
  I234mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  RI4mod <- r2(I4mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI14mod <- r2(I14mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  RI24mod <- r2(I24mod)$R2_marginal
  RI34mod <- r2(I34mod)$R2_marginal
  
  RI123mod <- r2(I123mod)$R2_marginal
  RI124mod <- r2(I124mod)$R2_marginal
  RI134mod <- r2(I134mod)$R2_marginal
  RI234mod <- r2(I234mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI234mod
  ub <- Rfullmod - RI134mod
  uc <- Rfullmod - RI124mod
  ud <- Rfullmod - RI123mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI34mod - ua - ub
  ac <- Rfullmod - RI24mod - ua - uc
  ad <- Rfullmod - RI23mod - ua - ud
  bc <- Rfullmod - RI14mod - ub - uc
  bd <- Rfullmod - RI13mod - ub - ud
  cd <- Rfullmod - RI12mod - uc - ud
  
  # Triplet intersections
  abc <- Rfullmod - RI4mod - ua - ub - uc- ab-ac-bc
  abd <- Rfullmod - RI3mod - ua - ub - ud- ab-ad-bd
  acd <- Rfullmod - RI2mod - ua - uc - ud- ac-ad-cd
  bcd <- Rfullmod - RI1mod - ub - uc - ud- bc-bd-cd
  
  # 4-way intersection
  abcd <- RI1mod + RI2mod - RI12mod -ab- abc-abd
  
  va <- sum(ua+ub+uc+ud+
              ab+ac+ad+bc+bd+cd+
              abc+abd+acd+bcd+
              abcd)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "Land Dist." = round(max(0,ub), 3),
    "PC1" = round(max(0,uc), 3),
    "PC2" = round(max(0,ud), 3),
    
    # Pairwise intersections
    "Geo Dist.&Land Dist." = round(max(0,ab), 3),
    "Geo Dist.&PC1" = round(max(0,ac), 3),
    "Geo Dist.&PC2" = round(max(0,ad), 3),
    "Land Dist.&PC1" = round(max(0,bc), 3),
    "Land Dist.&PC2" = round(max(0,bd), 3),
    "PC1&PC2" = round(max(0,cd), 3),
    
    # Triplet intersections
    "Geo Dist.&Land Dist.&PC1" = round(max(0,abc), 3),
    "Geo Dist.&Land Dist.&PC2" = round(max(0,abd), 3),
    "Geo Dist.&PC1&PC2" = round(max(0,acd), 3),
    "Land Dist.&PC1&PC2" = round(max(0,bcd), 3),
    
    # 4-way intersection
    "Geo Dist.&Land Dist.&PC1&PC2" = round(max(0,abcd), 3)
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("indianred1", "#FDBF6F", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = TRUE,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 
varpartitionlm <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . - Walk-PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC1-PC2)
  I3mod <- update(fullmod, formula=. ~ . - Walk-Air-PC2)
  I4mod <- update(fullmod, formula=. ~ . - Walk-PC1-Air)
  
  I12mod <- update(fullmod, formula=. ~ . - PC1-PC2)
  I13mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I14mod <- update(fullmod, formula=. ~ . - Walk-PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I24mod <- update(fullmod, formula=. ~ . - Air-PC1)
  I34mod <- update(fullmod, formula=. ~ . - Walk-Air)
  
  I123mod <- update(fullmod, formula=. ~ . - PC2)
  I124mod <- update(fullmod, formula=. ~ . - PC1)
  I134mod <- update(fullmod, formula=. ~ . - Walk)
  I234mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted 
  RI4mod <- r2(I4mod)$R2_adjusted 
  
  RI12mod <- r2(I12mod)$R2_adjusted 
  RI13mod <- r2(I13mod)$R2_adjusted 
  RI14mod <- r2(I14mod)$R2_adjusted 
  RI23mod <- r2(I23mod)$R2_adjusted 
  RI24mod <- r2(I24mod)$R2_adjusted 
  RI34mod <- r2(I34mod)$R2_adjusted 
  
  RI123mod <- r2(I123mod)$R2_adjusted 
  RI124mod <- r2(I124mod)$R2_adjusted 
  RI134mod <- r2(I134mod)$R2_adjusted 
  RI234mod <- r2(I234mod)$R2_adjusted 
  
  # Individual sets
  ua <- Rfullmod - RI234mod
  ub <- Rfullmod - RI134mod
  uc <- Rfullmod - RI124mod
  ud <- Rfullmod - RI123mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI34mod - ua - ub
  ac <- Rfullmod - RI24mod - ua - uc
  ad <- Rfullmod - RI23mod - ua - ud
  bc <- Rfullmod - RI14mod - ub - uc
  bd <- Rfullmod - RI13mod - ub - ud
  cd <- Rfullmod - RI12mod - uc - ud
  
  # Triplet intersections
  abc <- Rfullmod - RI4mod - ua - ub - uc- ab-ac-bc
  abd <- Rfullmod - RI3mod - ua - ub - ud- ab-ad-bd
  acd <- Rfullmod - RI2mod - ua - uc - ud- ac-ad-cd
  bcd <- Rfullmod - RI1mod - ub - uc - ud- bc-bd-cd
  
  # 4-way intersection
  abcd <- RI1mod + RI2mod - RI12mod -ab- abc-abd
  
  va <- sum(ua+ub+uc+ud+
              ab+ac+ad+bc+bd+cd+
              abc+abd+acd+bcd+
              abcd)
  residua <- 1- Rfullmod
  va+residua
  venn_data <- c(
    # Individual sets
    "Geo Dist." = round(max(0,ua), 3),
    "Land Dist." = round(max(0,ub), 3),
    "PC1" = round(max(0,uc), 3),
    "PC2" = round(max(0,ud), 3),
    
    # Pairwise intersections
    "Geo Dist.&Land Dist." = round(max(0,ab), 3),
    "Geo Dist.&PC1" = round(max(0,ac), 3),
    "Geo Dist.&PC2" = round(max(0,ad), 3),
    "Land Dist.&PC1" = round(max(0,bc), 3),
    "Land Dist.&PC2" = round(max(0,bd), 3),
    "PC1&PC2" = round(max(0,cd), 3),
    
    # Triplet intersections
    "Geo Dist.&Land Dist.&PC1" = round(max(0,abc), 3),
    "Geo Dist.&Land Dist.&PC2" = round(max(0,abd), 3),
    "Geo Dist.&PC1&PC2" = round(max(0,acd), 3),
    "Land Dist.&PC1&PC2" = round(max(0,bcd), 3),
    
    # 4-way intersection
    "Geo Dist.&Land Dist.&PC1&PC2" = round(max(0,abcd), 3)
  )
  # Generate and plot the diagram
  fit <- venn(venn_data)
  p <- plot(fit,
            quantities = TRUE,     # show counts
            fills = list(fill = c("indianred1", "#FDBF6F", "cadetblue3", "#A6CEE3")),
            legend = F,
            #main = "5-Set Venn Diagram",
            labels = F, #list(font = 2),  # bold labels
            adjust_labels = T,
            shape = "ellipse"
  )
  
  for(i in seq_along(p$children$canvas.grob$children$diagram.grob.1$children$tags$children)){
    o <- p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]]
    if(!is.null(o)){
      if ( o$children[[paste0("tag.quantity.",i)]]$label == 0){ 
        o$children[[paste0("tag.quantity.",i)]]$label <- " "
        p$children$canvas.grob$children$diagram.grob.1$children$tags$children[[paste0("tag.number.", i)]] <- o
      }
    }
  }
  p
} 

varpartitionwohormone(glmmIsLeakyF)
varpartitionwohormone(glmmIsLeakyM)
varpartitionwohormone(glmmREF)
varpartitionwohormone(glmmREM)
varpartitionlm(sweight)

varpartitionwohormone <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_marginal
  
  I1mod <- update(fullmod, formula= . ~ . - Walk-PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC1-PC2)
  I3mod <- update(fullmod, formula=. ~ . - Walk-Air-PC2)
  I4mod <- update(fullmod, formula=. ~ . - Walk-PC1-Air)
  
  I12mod <- update(fullmod, formula=. ~ . - PC1-PC2)
  I13mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I14mod <- update(fullmod, formula=. ~ . - Walk-PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I24mod <- update(fullmod, formula=. ~ . - Air-PC1)
  I34mod <- update(fullmod, formula=. ~ . - Walk-Air)
  
  I123mod <- update(fullmod, formula=. ~ . - PC2)
  I124mod <- update(fullmod, formula=. ~ . - PC1)
  I134mod <- update(fullmod, formula=. ~ . - Walk)
  I234mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_marginal
  RI2mod <- r2(I2mod)$R2_marginal
  RI3mod <- r2(I3mod)$R2_marginal
  RI4mod <- r2(I4mod)$R2_marginal
  
  RI12mod <- r2(I12mod)$R2_marginal
  RI13mod <- r2(I13mod)$R2_marginal
  RI14mod <- r2(I14mod)$R2_marginal
  RI23mod <- r2(I23mod)$R2_marginal
  RI24mod <- r2(I24mod)$R2_marginal
  RI34mod <- r2(I34mod)$R2_marginal
  
  RI123mod <- r2(I123mod)$R2_marginal
  RI124mod <- r2(I124mod)$R2_marginal
  RI134mod <- r2(I134mod)$R2_marginal
  RI234mod <- r2(I234mod)$R2_marginal
  
  # Individual sets
  ua <- Rfullmod - RI234mod
  ub <- Rfullmod - RI134mod
  uc <- Rfullmod - RI124mod
  ud <- Rfullmod - RI123mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI34mod - ua - ub
  ac <- Rfullmod - RI24mod - ua - uc
  ad <- Rfullmod - RI23mod - ua - ud
  bc <- Rfullmod - RI14mod - ub - uc
  bd <- Rfullmod - RI13mod - ub - ud
  cd <- Rfullmod - RI12mod - uc - ud
  
  # Triplet intersections
  abc <- Rfullmod - RI4mod - ua - ub - uc- ab-ac-bc
  abd <- Rfullmod - RI3mod - ua - ub - ud- ab-ad-bd
  acd <- Rfullmod - RI2mod - ua - uc - ud- ac-ad-cd
  bcd <- Rfullmod - RI1mod - ub - uc - ud- bc-bd-cd
  
  # 4-way intersection
  abcd <- RI1mod + RI2mod - RI12mod -ab- abc-abd
  
  va <- sum(ua+ub+uc+ud+
              ab+ac+ad+bc+bd+cd+
              abc+abd+acd+bcd+
              abcd)
  residua <- 1- Rfullmod
  va+residua
  # Generate and plot the diagram
  x <- c(residua,ua,ub, ab ,uc,ud, cd, 0)
  return(x)
} 
varpartitionlm <- function(fullmodel){
  fullmod <- fullmodel
  Rfullmod <- r2(fullmod)$R2_adjusted
  
  I1mod <- update(fullmod, formula= . ~ . - Walk-PC1-PC2)
  I2mod <- update(fullmod, formula=. ~ . - Air-PC1-PC2)
  I3mod <- update(fullmod, formula=. ~ . - Walk-Air-PC2)
  I4mod <- update(fullmod, formula=. ~ . - Walk-PC1-Air)
  
  I12mod <- update(fullmod, formula=. ~ . - PC1-PC2)
  I13mod <- update(fullmod, formula=. ~ . - Walk-PC2)
  I14mod <- update(fullmod, formula=. ~ . - Walk-PC1)
  I23mod <- update(fullmod, formula=. ~ . - Air-PC2)
  I24mod <- update(fullmod, formula=. ~ . - Air-PC1)
  I34mod <- update(fullmod, formula=. ~ . - Walk-Air)
  
  I123mod <- update(fullmod, formula=. ~ . - PC2)
  I124mod <- update(fullmod, formula=. ~ . - PC1)
  I134mod <- update(fullmod, formula=. ~ . - Walk)
  I234mod <- update(fullmod, formula=. ~ . - Air)
  
  RI1mod <- r2(I1mod)$R2_adjusted
  RI2mod <- r2(I2mod)$R2_adjusted
  RI3mod <- r2(I3mod)$R2_adjusted 
  RI4mod <- r2(I4mod)$R2_adjusted 
  
  RI12mod <- r2(I12mod)$R2_adjusted 
  RI13mod <- r2(I13mod)$R2_adjusted 
  RI14mod <- r2(I14mod)$R2_adjusted 
  RI23mod <- r2(I23mod)$R2_adjusted 
  RI24mod <- r2(I24mod)$R2_adjusted 
  RI34mod <- r2(I34mod)$R2_adjusted 
  
  RI123mod <- r2(I123mod)$R2_adjusted 
  RI124mod <- r2(I124mod)$R2_adjusted 
  RI134mod <- r2(I134mod)$R2_adjusted 
  RI234mod <- r2(I234mod)$R2_adjusted 
  
  # Individual sets
  ua <- Rfullmod - RI234mod
  ub <- Rfullmod - RI134mod
  uc <- Rfullmod - RI124mod
  ud <- Rfullmod - RI123mod
  
  # Pairwise intersections
  ab <- Rfullmod - RI34mod - ua - ub
  ac <- Rfullmod - RI24mod - ua - uc
  ad <- Rfullmod - RI23mod - ua - ud
  bc <- Rfullmod - RI14mod - ub - uc
  bd <- Rfullmod - RI13mod - ub - ud
  cd <- Rfullmod - RI12mod - uc - ud
  
  # Triplet intersections
  abc <- Rfullmod - RI4mod - ua - ub - uc- ab-ac-bc
  abd <- Rfullmod - RI3mod - ua - ub - ud- ab-ad-bd
  acd <- Rfullmod - RI2mod - ua - uc - ud- ac-ad-cd
  bcd <- Rfullmod - RI1mod - ub - uc - ud- bc-bd-cd
  
  # 4-way intersection
  abcd <- RI1mod + RI2mod - RI12mod -ab- abc-abd
  
  va <- sum(ua+ub+uc+ud+
              ab+ac+ad+bc+bd+cd+
              abc+abd+acd+bcd+
              abcd)
  residua <- 1- Rfullmod
  va+residua
  x <- c(residua,ua,ub, ab ,uc,ud, cd, 0)
  return(x)
} 

LNF <- print(varpartitionwohormone(glmmIsLeakyF))
LNM <- print(varpartitionwohormone(glmmIsLeakyM))
REF <- print(varpartitionwohormone(glmmREF))
REM <- print(varpartitionwohormone(glmmREM))
SW <- print(varpartitionlm(sweight))

### VAR SUMMARY PLOTS
varfull <- data.frame(LNF,LNM, REF, REM, SW )
colnames(varfull) <- c("Leakiness probablity F","Leakiness probablity M", "Reproductive effort F" ,"Reproductive effort M" , 
                       "Seed weight")
varfull[8,] <- 1- colSums(varfull)
rownames(varfull) <- c("residuals","ua","ub", "ab" , "uc","ud", "cd", "Others_joint")
varfull[varfull < 0] <- 0
varfull$ID <- c("residuals","Geo Dist." , "Land Dist.", "Dist. joint", "PC1","PC2", "PC joint", "Others joint")
varfull2 <- varfull[-1,]
varfull2freq <- apply(varfull2[,-6],2,function(x){x/sum(x)})
varfull2freq <- data.frame(varfull2freq)
varfull2freq$ID <- c("Geo Dist." , "Land Dist.", "Dist. joint", "PC1","PC2", "PC joint", "Others joint")
varfullstack <- melt(varfull2freq, id.vars=c('ID'), variable.name='Response_Trait', value.name='Explained_Marginal_Variance')
ggplot(varfullstack, aes(fill=factor(ID, levels=c("Others joint" ,
                                                  "PC joint",
                                                  "PC2",
                                                  "PC1",
                                                  "Dist. joint",
                                                  "Land Dist.",
                                                  "Geo Dist.")), y=Explained_Marginal_Variance, x=Response_Trait)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#999999", "#33A1C9", "#A6CEE3", "cadetblue3", "#FF7F00", "#FDBF6F", "indianred1")) +
  labs(fill = "Predictors") + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank())

### PC1

glmmIsLeakyM <- glmmTMB(is.leaky  ~ PC1 * Hormone.Tag  + 
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ PC1 * Hormone.Tag + 
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmLDM <- glmmTMB(.Female.flowers.or.fruits ~ PC1 * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dfmalesClimate, 
                   family = nbinom2, ziformula = ~ PC1)
glmmLDF <- glmmTMB(.Male.flowers  ~ PC1 * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dffemalesClimate, 
                   family = nbinom2, ziformula = ~ PC1)
glmmREM <- glmmTMB(.Male.flowers  ~  PC1 * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  PC1 * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)
glmmBiomassF <- glmmTMB(biomassSC  ~ PC1 * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
glmmBiomassM <- glmmTMB(biomassSC  ~ PC1 * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightM <- glmmTMB(Height..cm.  ~ PC1 * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightF <- glmmTMB(Height..cm.  ~ PC1 * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
sweight <- glmmTMB(weight1seed  ~ PC1, 
                   data = dfSClimate)
sarea <- glmmTMB(Area  ~ 
                   PC1 +
                   (1|ID)
                 , 
                 data = dfClimate)

### PC2 

glmmIsLeakyM <- glmmTMB(is.leaky  ~ PC2 * Hormone.Tag  + 
                          (1|Population) + (1|Table)  , data = dfmalesClimate, family = binomial)
glmmIsLeakyF <- glmmTMB(is.leaky  ~ PC2 * Hormone.Tag + 
                          (1|Population) + (1|Table)  , data = dffemalesClimate, family = binomial)
glmmLDM <- glmmTMB(.Female.flowers.or.fruits ~ PC2 * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dfmalesClimate, 
                   family = nbinom2, ziformula = ~ PC2)
glmmLDF <- glmmTMB(.Male.flowers  ~ PC2 * Hormone.Tag + 
                     (1|Population) + (1|Table), 
                   offset = log(biomass.overRE.common), 
                   data = dffemalesClimate, 
                   family = nbinom2, ziformula = ~ PC2)
glmmREM <- glmmTMB(.Male.flowers  ~  PC2 * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dfmalesphenotypingClimate , offset = log(biomass), family = nbinom2)
glmmREF <- glmmTMB(.Female.flowers.or.fruits  ~  PC2 * Hormone.Tag + 
                     (1|Population)+ (1|Table), data = dffemalesphenotypingClimate, offset = log(biomass), family = nbinom2)
glmmBiomassF <- glmmTMB(biomassSC  ~ PC2 * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
glmmBiomassM <- glmmTMB(biomassSC  ~ PC2 * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightM <- glmmTMB(Height..cm.  ~ PC2 * Hormone.Tag + (1|Population) + (1|Table), data = dfmales)
glmmHeightF <- glmmTMB(Height..cm.  ~ PC2 * Hormone.Tag + (1|Population) + (1|Table), data = dffemales)
sweight <- glmmTMB(weight1seed  ~ PC2, 
                   data = dfSClimate)
sarea <- glmmTMB(Area  ~ 
                   PC2 +
                   (1|ID)
                 , 
                 data = dfClimate)

######### MAKE THE MAP ###################

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(maps)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

data <- read.csv("GPS26withdistancesMINFERRY.csv", header=TRUE, sep=",") 

my_sf <- st_as_sf(data, coords = c('E', 'N'))

my_sf <- st_set_crs(my_sf, crs = 4326)

(data <- st_as_sf(data, coords = c("E", "N"), 
                  crs = 4326, agr = "constant"))

ggplot(data = world) +
  geom_sf(fill = "white", color = "darkgrey") +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 55), expand = FALSE) +
  geom_sf(data = data, size = 5, 
          shape = 20, color = "#D55E00") +
  geom_point(x = 35.013214, y = 32.781618, size = 5, 
             shape = 20, color = "black") +
  theme_linedraw() +
  theme(panel.grid.major = element_line(colour = "transparent"))+
  theme(panel.background = element_rect(fill = "lightblue", colour = "grey50")) + theme(panel.grid.major = element_blank(),
                                                                                        panel.grid.minor = element_blank())+
  coord_sf(xlim = c(-10, 40), ylim = c(30, 55), expand = FALSE) +
  #geom_sf_text(data = data, aes(label = ID), size = 2.5, position = position_nudge(x = 0.5, y = 0.4))+
  ggrepel::geom_label_repel(
    data = data,
    aes(label = ID, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0
  )

ggplot(nc) +
  geom_sf() +
  ggrepel::geom_label_repel(
    data = head(nc),
    aes(label = NAME, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0
  )
