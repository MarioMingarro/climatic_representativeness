library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/natural_parks_2.shp")


# Reference system
reference_system <- "+init=epsg:4326"
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))


plot(polygon)

# Climatic representativeness -----
present_climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/ANALISIS/CLIMA/PRESENT", "\\.tif$", full.names = T))
present_climatic_variables <- projectRaster(climatic_variables, crs = reference_system)

# Climatic representativeness -----
future_climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/ANALISIS/CLIMA/FUTURO", "\\.tif$", full.names = T))
future_climatic_variables <- projectRaster(climatic_variables, crs = reference_system)

# Crop raster to study area
present_climatic_variables <-  raster::mask(crop(present_climatic_variables, study_area), study_area)
future_climatic_variables <-  raster::mask(crop(future_climatic_variables, study_area), study_area)

# Geomorphological representativeness ----
geomorphologic_variables <-  raster::stack(list.files("T:/MODCLIM_R_DATA/analysis/geomophological/", ".tif", full.names = T))
geomorphologic_variables <- projectRaster(geomorphologic_variables, crs = reference_system)

# Crop raster to study area
geomorphologic_variables <-  mask(crop(geomorphologic_variables, study_area), study_area)


# Climate change variables
#climatic_change_variables <-

# PCA ----
# Climatic representativeness
PCA_present_climatic_variables <- RStoolbox::rasterPCA(present_climatic_variables, spca = TRUE)
PCA_future_climatic_variables <- RStoolbox::rasterPCA(future_climatic_variables, spca = TRUE)

# Geomorphological representativeness 
PCA_geomorphologic_variables <- RStoolbox::rasterPCA(geomorphologic_variables, spca = TRUE)



# PCA results
summary(PCA_present_climatic_variables$model)
summary(PCA_future_climatic_variables$model)



round(PCA_present_climatic_variables$model$loadings[,1:6],3)
round(PCA_future_climatic_variables$model$loadings[,1:6],3)

selected_PCA_PCA_present_climatic_variables <- PCA_present_climatic_variables$map[[1:5]]
selected_PCA_PCA_PCA_future_climatic_variables<- PCA_future_climatic_variables$map[[1:5]]
selected_PCA_geomorphological_variable <- geomorphologic_PCA$map[[1:3]]

plot(selected_PCA_PCA_present_climatic_variables)


writeRaster(selected_PCA_PCA_present_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
writeRaster(selected_PCA_PCA_PCA_future_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
writeRaster(selected_PCA_geomorphological_variable, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")


###########################################################################################################
####################################################################################

PCA_present <- raster("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
PCA_future <- raster("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
PCA_geo <- raster("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")
PCA_geo <- resample(PCA_geo, PCA_present, method='bilinear')

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/natural_parks_2.shp")



PCA_present <- raster::as.data.frame(PCA_present, xy = TRUE)
PCA_present <- na.omit(PCA_present)
PCA_present <- mutate(PCA_present, Periodo = c("Presente"))

PCA_future <- raster::as.data.frame(PCA_future, xy = TRUE)
PCA_future <- na.omit(PCA_future)
PCA_future <- mutate(PCA_future, Periodo = c("Futuro"))

PCA_geo <- raster::as.data.frame(PCA_geo, xy = TRUE)
PCA_geo <- na.omit(PCA_geo)
PCA_geo <- mutate(PCA_geo, Periodo = c("Presente"))

data <- left_join(PCA_geo, PCA_present, by = c("x", "y"))
data <- left_join(data, PCA_future, by = c("x", "y"))
data <- na.omit(data)
colnames(data) <- c("x", "y", "a", "a", "a", "a", "a", "a")

data2 <- data[,1:4]
data2 <- rbind(data2, data[,c(1, 2, 5, 6)] )
data2 <- rbind(data2, data[,c(1, 2, 7, 8)] )
# Create raster with 3 first PCA factors

raster_present <- raster::stack(selected_PCA_PCA_present_climatic_variables, selected_PCA_PCA_PCA_future_climatic_variables, selected_PCA_geomorphological_variable) # n of factors



# Raster to datafame
data_present <- raster::as.data.frame(raster_present, xy = TRUE)
data_present <- na.omit(data_present)

selected_PCA_PCA_present_climatic_variables <- raster::as.data.frame(selected_PCA_PCA_present_climatic_variables, xy = TRUE)
selected_PCA_PCA_present_climatic_variables <- na.omit(selected_PCA_PCA_present_climatic_variables)
selected_PCA_PCA_present_climatic_variables <- mutate(selected_PCA_PCA_present_climatic_variables, Periodo = c("Presente"))

selected_PCA_PCA_PCA_future_climatic_variables <- raster::as.data.frame(selected_PCA_PCA_PCA_future_climatic_variables, xy = TRUE)
selected_PCA_PCA_PCA_future_climatic_variables <- na.omit(selected_PCA_PCA_PCA_future_climatic_variables)
selected_PCA_PCA_PCA_future_climatic_variables <- mutate(selected_PCA_PCA_PCA_future_climatic_variables, Periodo = c("Futuro"))
selected_PCA_geomorphological_variables <- raster::as.data.frame(selected_PCA_geomorphological_variables, xy = TRUE)
selected_PCA_geomorphological_variables <- na.omit(selected_PCA_geomorphological_variables)
data_presente <- mutate(data_presente, Periodo = c("Presente"))


#nrow(polygon)
names <- polygon$ORIG_NAME

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

plot(mh_present)
nlayers(mh_present)

# Umbral para establecer el corte percentil 90
mh_present_umbral <- raster::brick()

for (i in 1:nlayers(mh_present)){
  aa <- polygon$geometry[i]
  mh_poligono <- mask(crop(mh_present[[i]], polygon[i,]), polygon[i,])
  mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
  mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.90))
  mh_presente_bin <- reclassify(mh_present[[i]], c(mh_poligono,Inf,NA))
  mh_present_umbral <- raster::stack(mh_present_umbral, mh_presente_bin)
}

plot(mh_poligono)
plot(mh_present)
plot(mh_present_umbral)

mh_poligono <- mask(crop(mh_present, polygon), polygon)
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
