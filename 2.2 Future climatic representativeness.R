library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)

# Load data
climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/ENVIREM/", ".tif", full.names = T))

library(stringr)
nombres <- list.files("T:/MODCLIM_R_DATA/ENVIREM/", ".tif", full.names = T)


study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))
polygon <- read_sf(list.files("T:/MODCLIM_R_DATA/LIC/", "\\.shp$", full.names = T))
polygon <- read_sf("T:/MODCLIM_R_DATA/LIC/LIC_SELECT.shp")
polygon <- read_sf("T:/MODCLIM_R_DATA/THIC/THIC.shp")

study_area <- st_transform(study_area, crs(projection(climatic_variables)))
polygon <- st_transform(polygon, crs(projection(climatic_variables)))


# Crop raster to study area and polygon(centroid)
beginCluster()
climatic_variables <-  mask(crop(climatic_variables, study_area), study_area)
endCluster()

# PCA
beginCluster()
PCA <- RStoolbox::rasterPCA(climatic_variables, spca = TRUE)
endCluster()

# PCA results
round(PCA$model$loadings[,1:5],3)
kk <- summary(PCA$model)
kk <- kk$loadings
kk <- as.data.frame(kk)
PCA$model

# Create raster with 3 first PCA factors
raster_present <- PCA$map[[1:5]] # n of factors


# Raster to datafame

data_present <- raster::as.data.frame(raster_present, xy = TRUE)
data_present <- na.omit(data_present)

#nrow(polygon)
names <- polygon$THIC

mh_f <- data.frame(matrix(1,    # Create empty data frame
                  nrow = nrow(data_present),
                  ncol = length(names)))

names(mh_f) <- names
beginCluster()
for (i in 1:nrow(polygon)){
  pol <- polygon[i,]
  raster_polygon <- mask(crop(raster_present, pol), pol)
  data_polygon <- raster::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  # Mahalanobis distance
  
  mh <- mahalanobis(data_present[,3:length(data_present)], 
                    colMeans(data_polygon[,3:length(data_polygon)]), 
                    cov(data_present[,3:length(data_present)]), 
                    inverted = F)
  
  # Agregar información espacial al reusltado de mh
  mh_f[,i] <- mh
}

# Agregar información espacial al reusltado de mh
mh <- cbind(data_present[,1:2], mh_f)


#Creamos raster
mh_present <- raster::brick()

for(j in 3:length(mh)){
  aa <- rasterFromXYZ(mh[, c(1:2,j)])
  mh_present <- raster::stack(mh_present, aa)
}

endCluster()

plot(mh_present)

# Umbral para establecer el corte percentil 90

-
plot(mh_present[[1]])
plot(mh_polygon)

names_THIC <- readxl::read_xlsx("T:/MODCLIM_R_DATA/THIC/THIC_names.xlsx")
plot(mh_present_bin_f[[1]])
names(mh_present_bin_f)

monthly_pcp