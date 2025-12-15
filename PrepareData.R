library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyverse)

df <- read_xlsx("phenobiomassfull.xlsx", .name_repair = "universal")
dfS <- read_xlsx("Seed phenotype.xlsx", .name_repair = "universal")
gps <- read.csv("GPS26withdistancesMINFERRY.csv")
file.list <- list.files(path = "seeds pictures",pattern='*.csv')
dftest <- read.csv(paste0("seeds pictures/",file.list[2]))
gps$ID <- as.character(gps$ID)
popID <- gps$ID
siteClimate <- read.csv("siteClimatePCAnew02.csv") ## OLD ONE is file siteClimatePCAnew.csv

df$Height..cm. <- as.character(sub("\\,",".", df$Height..cm.))
df$Sick <- as.character(sub("\\,",".", df$Sick))
df$Broken <- as.character(sub("\\,",".", df$Broken))
df$Yellow.Green <- as.character(sub("\\,",".", df$Yellow.Green))
df[df$ID.Tag == "16.22",]$Sick <- "1"
df$Height..cm. <- as.numeric(df$Height..cm.)
df$Sick <- as.numeric(df$Sick)
df$Broken <- as.numeric(df$Broken)
df$Yellow.Green <- as.numeric(df$Yellow.Green)
df$Air <- NA
df$Walk <- NA
df[df$ID.Tag == "6.87",]$biomass <- NA
df$biomassSC <- round(scale(df$biomass))
for (i in 1:nrow(df)) {
  df[i,]$Air <- gps[gps$ID == df[i,]$Population,]$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular
  df[i,]$Walk <- gps[gps$ID == df[i,]$Population,]$WalktoHaifavisGGmap
}

mistakensex <- c("22.67", "20.5", "27.76", "27.106", "22.80", "27.5", "19.3", "21.11", "19.10") #phenotyping inconsistent with exp design

for (i in mistakensex) {
  df[df$ID.Tag == i,]$.Male.flowers <- NA
  df[df$ID.Tag == i,]$.Female.flowers.or.fruits <- NA
}

df <- df  %>% drop_na(Hormone.Tag) 
df$.Male.flowers <- as.numeric(df$.Male.flowers)
df$.Female.flowers.or.fruits <- as.numeric(df$.Female.flowers.or.fruits)
dfmales <- df[df$Greenhouse == 1,]
dffemales <- df[df$Greenhouse == 2,]
dffemalesphenotyping <- dffemales[dffemales$Phenotyping.tag == "x" & dffemales$.Male.flowers == 0,]
dfmalesphenotyping <- dfmales[dfmales$Phenotyping.tag == "x" & dfmales$.Female.flowers.or.fruits == 0,]
dffemalesphenotyping$femaleRE <- dffemalesphenotyping$.Female.flowers.or.fruits / dffemalesphenotyping$biomass
dfmalesphenotyping$maleRE <- dfmalesphenotyping$.Male.flowers / dfmalesphenotyping$biomass
dffemalesphenotypingwohormone <- dffemalesphenotyping[dffemalesphenotyping$Hormone.Tag == "-",]
dfmalesphenotypingwohormone <- dfmalesphenotyping[dfmalesphenotyping$Hormone.Tag == "-",]
dffemalesphenotypingwhormone <- dffemalesphenotyping[dffemalesphenotyping$Hormone.Tag == "x",]
dfmalesphenotypingwhormone <- dfmalesphenotyping[dfmalesphenotyping$Hormone.Tag == "x",]
dffemalesphenotypingwohormone <- dffemalesphenotypingwohormone  %>% drop_na(femaleRE) 
dfmalesphenotypingwohormone <- dfmalesphenotypingwohormone  %>% drop_na(maleRE) 
dffemalesphenotypingwhormone <- dffemalesphenotypingwhormone %>% drop_na(femaleRE) 
dfmalesphenotypingwhormone <- dfmalesphenotypingwhormone %>% drop_na(maleRE) 
dffemalesphenotyping <- dffemalesphenotyping  %>% drop_na(femaleRE) 
dfmalesphenotyping <- dfmalesphenotyping  %>% drop_na(maleRE) 

REperpop <- gps[,c(1,7,8)]
REperpop$femaleREwoH <- NA
REperpop$femaleREwH <- NA
REperpop$maleREwoH <- NA
REperpop$maleREwH <- NA
REperpop$femaleREboth <- NA
REperpop$maleREboth <- NA


for (i in popID) {
  REperpop[REperpop$ID == i,]$femaleREwoH <- mean(dffemalesphenotypingwohormone[dffemalesphenotypingwohormone$Population == i,]$femaleRE)
  REperpop[REperpop$ID == i,]$femaleREwH <- mean(dffemalesphenotypingwhormone[dffemalesphenotypingwhormone$Population == i,]$femaleRE)
  REperpop[REperpop$ID == i,]$maleREwoH <- mean(dfmalesphenotypingwohormone[dfmalesphenotypingwohormone$Population == i,]$maleRE)
  REperpop[REperpop$ID == i,]$maleREwH <- mean(dfmalesphenotypingwhormone[dfmalesphenotypingwhormone$Population == i,]$maleRE)
  REperpop[REperpop$ID == i,]$femaleREboth <- mean(dffemalesphenotyping[dffemalesphenotyping$Population == i,]$femaleRE)
  REperpop[REperpop$ID == i,]$maleREboth <- mean(dfmalesphenotyping[dfmalesphenotyping$Population == i,]$maleRE)
}

df$leaky.degree <- NA

for (i in 1:nrow(df)) {
  if (df[i,]$Greenhouse == "1"){
    if (df[i,]$Hormone.Tag == "x"){
      df[i,]$leaky.degree <- (df[i,]$.Female.flowers.or.fruits / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$femaleREwH
    } else {
      df[i,]$leaky.degree <- (df[i,]$.Female.flowers.or.fruits / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$femaleREwoH
    }
  } else{
    if (df[i,]$Hormone.Tag == "x"){
      df[i,]$leaky.degree <- (df[i,]$.Male.flowers / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$maleREwH
    } else {
      df[i,]$leaky.degree <- (df[i,]$.Male.flowers / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$maleREwoH
    }
  }
}

REperpop$avleaky.degreeMwH <- NA 
REperpop$avleaky.degreeMwoH <- NA 
REperpop$avleaky.degreeFwH <- NA 
REperpop$avleaky.degreeFwoH <- NA 

df <- df  %>% drop_na(leaky.degree) 
for (i in popID) {
  REperpop[REperpop$ID == i,]$avleaky.degreeMwH <- mean(df[df$Population == i & df$Greenhouse == "1" & df$Hormone.Tag == "x",]$leaky.degree) 
  REperpop[REperpop$ID == i,]$avleaky.degreeMwoH <- mean(df[df$Population == i & df$Greenhouse == "1" & df$Hormone.Tag == "-",]$leaky.degree)  
  REperpop[REperpop$ID == i,]$avleaky.degreeFwH <- mean(df[df$Population == i & df$Greenhouse == "2" & df$Hormone.Tag == "x",]$leaky.degree) 
  REperpop[REperpop$ID == i,]$avleaky.degreeFwoH <- mean(df[df$Population == i & df$Greenhouse == "2" & df$Hormone.Tag == "-",]$leaky.degree) 
}

### Combine hormone with no hormone

df$leaky.degree.comhormone <- NA

for (i in 1:nrow(df)) {
  if (df[i,]$Greenhouse == "1"){
    df[i,]$leaky.degree.comhormone <- (df[i,]$.Female.flowers.or.fruits / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$femaleREwoH
  } else{
    df[i,]$leaky.degree.comhormone <- (df[i,]$.Male.flowers / df[i,]$biomass) / REperpop[REperpop$ID == df[i,]$Population,]$maleREwoH
  }
}

REperpop$avleaky.degreeM <- NA 
REperpop$avleaky.degreeF <- NA 
df <- df  %>% drop_na(leaky.degree.comhormone) 
for (i in popID) {
  REperpop[REperpop$ID == i,]$avleaky.degreeM <- mean(df[df$Population == i & df$Greenhouse == "1",]$leaky.degree.comhormone) 
  REperpop[REperpop$ID == i,]$avleaky.degreeF <- mean(df[df$Population == i & df$Greenhouse == "2",]$leaky.degree.comhormone) 
}


### Leaky count number

allmales <- dfmales
allfemales <- dffemales
sumdf <- gps[,c(1,7,8)]
sumdf$leakymalescountnohormone <- 0
sumdf$leakyfemalescountnohormone <- 0
sumdf$leakymalescountwhormone <- 0
sumdf$leakyfemalescountwhormone <- 0
sumdf$sumMnoH <- 0
sumdf$sumMwH <- 0
sumdf$sumFnoH <- 0
sumdf$sumFwH <- 0

allmalecountleaky <- allmales
allmalecountleaky$countleaky <- replace(allmalecountleaky$.Female.flowers.or.fruits, allmalecountleaky$.Female.flowers.or.fruits > 0, 1) 
allfemalecountleaky <- allfemales
allfemalecountleaky$countleaky <- replace(allfemalecountleaky$.Male.flowers, allfemalecountleaky$.Male.flowers > 0, 1) 
library(tidyr)
allmalecountleaky <- allmalecountleaky %>% drop_na(countleaky)
allfemalecountleaky <- allfemalecountleaky %>% drop_na(countleaky)

for (i in popID) {
  newdfMnoH <- allmalecountleaky[allmalecountleaky$Population == i & 
                                   allmalecountleaky$Hormone.Tag == "-",]
  newdfMwH <- allmalecountleaky[allmalecountleaky$Population == i & 
                                  allmalecountleaky$Hormone.Tag == "x",]
  newdfFnoH <- allfemalecountleaky[allfemalecountleaky$Population == i & 
                                     allfemalecountleaky$Hormone.Tag == "-",]
  newdfFwH <- allfemalecountleaky[allfemalecountleaky$Population == i & 
                                    allfemalecountleaky$Hormone.Tag == "x",]
  sumdf[sumdf$ID == i,]$leakymalescountnohormone <- sum(
    newdfMnoH$countleaky) 
  sumdf[sumdf$ID == i,]$leakymalescountwhormone <- sum(
    newdfMwH$countleaky) 
  sumdf[sumdf$ID == i,]$leakyfemalescountnohormone <- sum(
    newdfFnoH$countleaky) 
  sumdf[sumdf$ID == i,]$leakyfemalescountwhormone <- sum(
    newdfFwH$countleaky) 
  
  sumdf[sumdf$ID == i,]$sumMnoH <- length(
    newdfMnoH$countleaky) 
  sumdf[sumdf$ID == i,]$sumMwH <- length(
    newdfMwH$countleaky) 
  sumdf[sumdf$ID == i,]$sumFnoH <- length(
    newdfFnoH$countleaky) 
  sumdf[sumdf$ID == i,]$sumFwH <- length(
    newdfFwH$countleaky) 
  
}
sumdf$leakymalescount <- sumdf$leakymalescountnohormone + sumdf$leakymalescountwhormone
sumdf$leakyfemalescount <- sumdf$leakyfemalescountnohormone + sumdf$leakyfemalescountwhormone


df$is.leaky <- 0
df[df$leaky.degree >0, ]$is.leaky <- 1

## Deal with degree 
df$biomass.overRE.common <- NA
for (i in 1:nrow(df)) {
  if (df[i,]$Greenhouse == "1"){
    df[i,]$biomass.overRE.common <- (df[i,]$biomass) * REperpop[REperpop$ID == df[i,]$Population,]$femaleREwoH
  } else{
    df[i,]$biomass.overRE.common <- (df[i,]$biomass) * REperpop[REperpop$ID == df[i,]$Population,]$maleREwoH
  }
}
dffemales <- df[df$Greenhouse == "2",]
dfmales <- df[df$Greenhouse == "1",]

sumdf$avleakydegreecommhormoneM <- 0
sumdf$avleakydegreecommhormoneF <- 0

for (i in popID) {
  newdfM <- dfmales[dfmales$Population == i,]
  newdfF <- dffemales[dffemales$Population == i,]
  sumdf[sumdf$ID == i,]$avleakydegreecommhormoneM <- mean(
    newdfM$leaky.degree.comhormone) 
  sumdf[sumdf$ID == i,]$avleakydegreecommhormoneF <-  mean(
    newdfF$leaky.degree.comhormone) 
}

sumdf$avLDMminusF <- sumdf$avleakydegreecommhormoneM-sumdf$avleakydegreecommhormoneF

sumdf$avleakydegreecommhormoneMallL <- 0
sumdf$avleakydegreecommhormoneFallL <- 0
sumdf$avLDMminusFallL <- 0

for (i in popID) {
  newdfM <- dfmales[dfmales$Population == i & dfmales$is.leaky == 1,]
  newdfF <- dffemales[dffemales$Population == i & dffemales$is.leaky == 1,]
  sumdf[sumdf$ID == i,]$avleakydegreecommhormoneMallL <- mean(
    newdfM$leaky.degree.comhormone) 
  sumdf[sumdf$ID == i,]$avleakydegreecommhormoneFallL <-  mean(
    newdfF$leaky.degree.comhormone) 
}

sumdf[is.na(sumdf)] <- 0
sumdf$avLDMminusFallL <- sumdf$avleakydegreecommhormoneMallL-sumdf$avleakydegreecommhormoneFallL
sumdf$difLNFtoLNM <- (sumdf$leakyfemalescount / (sumdf$sumFwH + sumdf$sumFnoH)) - (sumdf$leakymalescount / (sumdf$sumMwH + sumdf$sumMnoH))

### REPRODUCTIVE EFFORT
dffemalesphenotyping$femaleRE <- dffemalesphenotyping$.Female.flowers.or.fruits / dffemalesphenotyping$biomass
dfmalesphenotyping$maleRE <- dfmalesphenotyping$.Male.flowers / dfmalesphenotyping$biomass

#### SEEDSSSSS

# ---- Seed weight ----
gps$ID <- as.character(gps$ID)
popID <- gps$ID
dfS$Air <- NA
dfS$Walk <- NA
for (i in 1:nrow(dfS)) {
  dfS[i,]$Air <- gps[gps$ID == dfS[i,]$ID,]$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular
  dfS[i,]$Walk <- gps[gps$ID == dfS[i,]$ID,]$WalktoHaifavisGGmap
}
dfS <- dfS  %>% drop_na(Weight.of.selected.seeds..g.) 

dfS$Weight.of.selected.seeds..g. <- as.numeric(dfS$Weight.of.selected.seeds..g.)
dfS$Selection.seeds.Number..max.50. <- as.numeric(dfS$Selection.seeds.Number..max.50.)
dfS$weight1seed <- dfS$Weight.of.selected.seeds..g./dfS$Selection.seeds.Number..max.50.

# ---- Seed size ----
df <- dftest[1,]
df$ID <- NA
for (i in file.list) {
  dfnew <- read.csv(paste0("seeds pictures/", i))
  dfnew$ID <- str_match(i, "pop\\s*(.*?)\\s*.csv")[,2]
  df <- rbind(df, dfnew)
}

df <- df[-1,]

df$Air <- NA
df$Walk <- NA
for (i in 1:nrow(df)) {
  df[i,]$Air <- gps[gps$ID == df[i,]$ID,]$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular
  df[i,]$Walk <- gps[gps$ID == df[i,]$ID,]$WalktoHaifavisGGmap
}

### CLIMATE
siteClimate$leakynumberM <- sumdf$leakymalescount / (sumdf$sumMwH + sumdf$sumMnoH)
siteClimate$leakynumberF <- sumdf$leakyfemalescount / (sumdf$sumFwH + sumdf$sumFnoH)
siteClimate$leakydegM <- sumdf$avleakydegreecommhormoneM
siteClimate$leakydegF <- sumdf$avleakydegreecommhormoneF
siteClimate$REM <- REperpop$maleREboth
siteClimate$REF <- REperpop$femaleREboth

dfmalesClimate <- dfmales
dffemalesClimate <- dffemales
dfmalesphenotypingClimate <- dfmalesphenotyping
dffemalesphenotypingClimate <- dffemalesphenotyping
sumdfClimate <- sumdf
dfSClimate <- dfS
dfClimate <- df

dfmalesClimate$PC1 <- NA
dffemalesClimate$PC1 <- NA
dfmalesphenotypingClimate$PC1 <- NA
dffemalesphenotypingClimate$PC1 <- NA
sumdfClimate$PC1 <- NA
dfSClimate$PC1 <- NA
dfClimate$PC1 <- NA

dfmalesClimate$PC2 <- NA
dffemalesClimate$PC2 <- NA
dfmalesphenotypingClimate$PC2 <- NA
dffemalesphenotypingClimate$PC2 <- NA
sumdfClimate$PC2 <- NA
dfSClimate$PC2 <- NA
dfClimate$PC2 <- NA

i <- 1
for (i in 1:nrow(dfmalesClimate)) {
  dfmalesClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dfmalesClimate[i,]$Population,]$PC1
  dfmalesClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dfmalesClimate[i,]$Population,]$PC2
} 

i <- 1
for (i in 1:nrow(dffemalesClimate)) {
  dffemalesClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dffemalesClimate[i,]$Population,]$PC1
  dffemalesClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dffemalesClimate[i,]$Population,]$PC2
} 

i <- 1
for (i in 1:nrow(dfmalesphenotypingClimate)) {
  dfmalesphenotypingClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dfmalesphenotypingClimate[i,]$Population,]$PC1
  dfmalesphenotypingClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dfmalesphenotypingClimate[i,]$Population,]$PC2
} 

i <- 1
for (i in 1:nrow(dffemalesphenotypingClimate)) {
  dffemalesphenotypingClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dffemalesphenotypingClimate[i,]$Population,]$PC1
  dffemalesphenotypingClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dffemalesphenotypingClimate[i,]$Population,]$PC2
} 

i <- 1
for (i in 1:nrow(sumdfClimate)) {
  sumdfClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == sumdfClimate[i,]$ID,]$PC1
  sumdfClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == sumdfClimate[i,]$ID,]$PC2
} 

i <- 1
for (i in 1:nrow(dfSClimate)) {
  dfSClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dfSClimate[i,]$ID,]$PC1
  dfSClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dfSClimate[i,]$ID,]$PC2
} 

i <- 1
for (i in 1:nrow(dfClimate)) {
  dfClimate[i,]$PC1 <- siteClimate[siteClimate$Site_ID == dfClimate[i,]$ID,]$PC1
  dfClimate[i,]$PC2 <- siteClimate[siteClimate$Site_ID == dfClimate[i,]$ID,]$PC2
} 

saveRDS(dfmalesClimate, file = "dfmalesClimate.rds")
saveRDS(dffemalesClimate, file = "dffemalesClimate.rds")
saveRDS(dfmalesphenotypingClimate, file = "dfmalesphenotypingClimate.rds")
saveRDS(dffemalesphenotypingClimate, file = "dffemalesphenotypingClimate.rds")
saveRDS(sumdfClimate, file = "sumdfClimate.rds")
saveRDS(dfSClimate, file = "dfSClimate.rds")
saveRDS(dfClimate, file = "dfClimate.rds")
