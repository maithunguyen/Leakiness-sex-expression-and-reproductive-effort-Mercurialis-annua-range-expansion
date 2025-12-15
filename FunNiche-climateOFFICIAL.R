#*********************************************************************************
# FunNiche
# climate
# Thu Nguyen, adapted from Shengman Lyu
# Date updated: 31.1.2025
#*********************************************************************************
# Load packages
library(readxl)
library(viridis) # V.0.6.4
library(tidyverse)
library(terra) # V.1.7-55
library(ggplot2) # V.3.4.4
library(readxl)
library(gtools)
library(corrplot)
library(pheatmap)
library(ggpubr)
library('corrr')
library(ggcorrplot)
#*****************************************************
# ---- Data ----
#*****************************************************
# FunNiche site information
site <- read.csv("GPS26withdistancesMINFERRY.csv")

# CHELSA data
raster_files <- mixedsort(list.files(path = "CHELSA", pattern = "",
                                     full.names = TRUE))
# gather all raster layers

#####
# gather all raster layers
rlist <- lapply(raster_files, rast)

# Check extents
sapply(rlist, ext)

#check to see if they are all spatial rasters
sapply(rlist, inherits, what = "SpatRaster")

#identify messed up rasters
ref <- rlist[[1]]
bad_idx <- which(!sapply(rlist, function(r) ext(r) == ext(ref)))
bad_idx #has different extension, need separated extraction

rlistclean <- rlist[-23:-20]

bioclim <- terra::rast(rlistclean)
bioclim2 <- terra::rast(rlist[20:23])

# check the coordinate system (CRS)
cat(crs(bioclim)) # cat() allows to clearly read the string from crs()

# Change the names of the rasters
names(bioclim) = sub("CHELSA_","",names(bioclim))
names(bioclim) = sub("_1981-2010_V.2.1","",names(bioclim))

names(bioclim2) = sub("CHELSA_","",names(bioclim2))
names(bioclim2) = sub("_1981-2010_V.2.1","",names(bioclim2))

#*****************************************************
# ---- Extract climate data ----
#*****************************************************
# Extract climatic data get the bioclimatic raster files
coord <- site %>% select(E, N)
plot(subset(bioclim, 1), xlim=c(min(site$E) - 1 ,max(site$E) + 1), ylim=c(min(site$N) - 1 ,max(site$N) + 1))
points(site$E, site$N)

## Extract bioclimatic values

site.climate2 <- na.omit(data.frame(terra::extract(bioclim2, coord)))
site.climate2 <- (data.frame(terra::extract(bioclim2, coord)))
site.climate2 <- site.climate2 %>%
  select_if(~ !any(is.na(.)))

site.climate <- na.omit(data.frame(terra::extract(bioclim, coord)))
site.climate <- (data.frame(terra::extract(bioclim, coord)))
site.climate <- site.climate %>%
  select_if(~ !any(is.na(.)))
site.climate <- cbind(site.climate, site.climate2[2:5]) # combining both types of extension

site.climate$Site_ID = site$ID
site.climate$walk = site$WalktoHaifavisGGmap
site.climate$air = site$AirKmToHaifa.32.781618..35.013214.via.distancefromto.netAnddistance.to.VINCENTYFomular


#*****************************************************
# ---- Climate selection PCA ----
#*****************************************************
dim(site.climate)
site.climate.representatives <- site.climate[,2:66]
data_normalized <- scale(site.climate.representatives)
colMeans(data_normalized)
data.pca <- prcomp(data_normalized)
summary(data.pca)


# Select 2 PCs 
chosenPCA <- data.pca$rotation[, 1:2]
site.climate.PCAed <- as.matrix(data_normalized) %*% chosenPCA
site.climate.representatives$PC1 <- site.climate.PCAed[,1]
site.climate.representatives$PC2 <- site.climate.PCAed[,2]
site.climate.representatives$Site_ID <- site.climate$Site_ID
site.climate.representatives$walk <- site.climate$walk
site.climate.representatives$air <- site.climate$air
write.csv(site.climate.representatives, "siteClimatePCAnew02.csv")
