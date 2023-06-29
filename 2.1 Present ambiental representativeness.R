library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/analysis/", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/THIC/THIC.shp")


# Reference system
reference_system <- "+init=epsg:4326"
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))
polygon <- subset(polygon, polygon$THIC == 9260)

plot(polygon)

# Climatic representativeness -----
climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/analysis/climatic", "\\.tif$", full.names = T))
climatic_variables <- projectRaster(climatic_variables, crs = reference_system)

# Crop raster to study area
climatic_variables <-  raster::mask(crop(climatic_variables, study_area), study_area)

# Geomorphological representativeness ----
geomorphologic_variables <-  raster::stack(list.files("T:/MODCLIM_R_DATA/analysis/geomophological/", ".tif", full.names = T))
geomorphologic_variables <- projectRaster(geomorphologic_variables, crs = reference_system)

# Crop raster to study area
geomorphologic_variables <-  mask(crop(geomorphologic_variables, study_area), study_area)


# Climate change variables
#climatic_change_variables <-
  
  # PCA ----
# Climatic representativeness
PCA_climatic_variables <- RStoolbox::rasterPCA(climatic_variables, spca = TRUE)

# Geomorphological representativeness 
PCA_geomorphologic_variables <- RStoolbox::rasterPCA(geomorphologic_variables, spca = TRUE)

# Climate change variables
PCA <- RStoolbox::rasterPCA(climatic_change_variables, spca = TRUE)

# PCA results
summary(PCA_climatic_variables$model)
summary(PCA_geomorphologic_variables$model)


round(PCA_climatic_variables$model$loadings[,1:6],3)
round(PCA_geomorphologic_variables$model$loadings[,1:3],3)

selected_PCA_climate_variable <- PCA_climatic_variables$map[[1:6]]
selected_PCA_geomorphological_variable <- PCA_geomorphologic_variables$map[[1:3]]
# Create raster with 3 first PCA factors

selected_PCA_geomorphological_variable <- resample(selected_PCA_geomorphological_variable, selected_PCA_climate_variable, method='bilinear')
raster_present <- raster::stack(selected_PCA_climate_variable, selected_PCA_geomorphological_variable) # n of factors



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

mh_f <- rasterFromXYZ(mh)

endCluster()

plot(mh_f)


# Umbral para establecer el corte percentil 90
mh_poligono <- mask(crop(mh_f, polygon), polygon)
mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.90))
mh_poligono

writeRaster(mh_presente_bin, "T:/MODCLIM_R_DATA/analysis/result/castaños_nuevo.tif")
writeRaster(mh_f, "T:/MODCLIM_R_DATA/analysis/result/castaños_mh_nuevo.tif")


# Selección de esas distancias en PI
# Presente ----
mh_presente_bin <- reclassify(mh_f, c(mh_poligono,Inf,NA))
plot(mh_presente_bin)
plot(polygon[1], add = T)
