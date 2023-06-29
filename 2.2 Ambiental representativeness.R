library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))
study_area <- raster("T:/MODCLIM_R_DATA/analysis/mask2.tif")
polygon <- read_sf(list.files("T:/MODCLIM_R_DATA/LIC/", "\\.shp$", full.names = T))
polygon <- read_sf("T:/MODCLIM_R_DATA/LIC/LIC_SELECT.shp")
polygon <- read_sf("T:/MODCLIM_R_DATA/THIC/THIC.shp")


# Reference system
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_transform(study_area, crs(reference_system))
study_area <- projectRaster(study_area, crs = reference_system)
polygon <- st_transform(polygon, crs(reference_system))

polygon <- polygon$THIC[8]

nombres <- list.files("T:/MODCLIM_R_DATA/ENVIREM/", ".tif", full.names = T)


# Climatic representativeness -----
climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/ENVIREM/", ".tif", full.names = T))
# Crop raster to study area
climatic_variables <-  mask(crop(climatic_variables, study_area), study_area)


# Geomorphological representativeness ----
geomorphologic_variables<-  raster::stack(list.files("T:/MODCLIM_R_DATA/analysis/geomophological/", ".tif", full.names = T))

# Crop raster to study area
geomorphologic_variables <-  mask(crop(geomorphologic_variables, climatic_variables[[1]]), climatic_variables[[1]])

# Climate change variables
climatic_change_variables <-
  
  # PCA ----
# Climatic representativeness
PCA_climatic_variables <- RStoolbox::rasterPCA(climatic_variables, spca = TRUE)

# Geomorphological representativeness 
PCA_geomorphologic_variables <- RStoolbox::rasterPCA(geomorphologic_variables, spca = TRUE)

# Climate change variables
PCA <- RStoolbox::rasterPCA(climatic_change_variables, spca = TRUE)

writeRaster(PCA_climatic_variables$map[[1]], "T:/MODCLIM_R_DATA/kk/climatic.tif" )
writeRaster(PCA_geomorphologic_variables$map[[1]], "T:/MODCLIM_R_DATA/kk/geomorf.tif" )

# PCA results
summary(PCA_climatic_variables$model)
summary(PCA_geomorphologic_variables$model)


round(PCA_climatic_variables$model$loadings[,1:5],3)
round(PCA_geomorphologic_variables$model$loadings[,1:5],3)

selected_PCA_climate_variable <- PCA_climatic_variables$map[[1:4]]
selected_PCA_geomorphological_variable <- PCA_geomorphologic_variables$map[[1:4]]
# Create raster with 3 first PCA factors
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