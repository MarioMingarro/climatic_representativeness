library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)

# Load data
climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/ENVIREM/", ".tif", full.names = T))

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))
polygon <- read_sf(list.files("T:/MODCLIM_R_DATA/shp", "\\.shp$", full.names = T))

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
summary(PCA$model)

# Create raster with 3 first PCA factors
raster_present <- PCA$map[[1:5]] # n of factors
raster_polygon <- mask(crop(raster_present, polygon), polygon)

# Raster to datafame

data_present <- raster::as.data.frame(raster_present, xy = TRUE)
data_present <- na.omit(data_present)

data_polygon <- raster::as.data.frame(raster_polygon, xy = TRUE)
data_polygon <- na.omit(data_polygon)


# Mahalanobis distance

mh <- mahalanobis(data_present[,3:length(data_present)], 
                  colMeans(data_polygon[,3:length(data_polygon)]), 
                  cov(data_present[,3:length(data_present)]), 
                  inverted = F)

# Agregar información espacial al reusltado de mh
mh <- cbind(data_present, mh)

mh_present <- mh[,c(1:2,8)]


#Creamos raster
mh_present <- rasterFromXYZ(mh_present)


plot(mh_present)


# Umbral para establecer el corte percentil 90
mh_polygon <- mask(crop(mh_present, polygon), polygon)
mh_polygon <- raster::as.data.frame(mh_polygon, xy = T)
mh_polygon <- quantile(na.omit(mh_polygon[,3]), probs = c(.95))
mh_polygon


# Selección de esas distancias en PI
# present ----
mh_present_bin <- reclassify(mh_present, c(mh_polygon,Inf,0))
plot(mh_present_bin)
