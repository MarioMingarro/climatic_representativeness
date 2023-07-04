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

# Geomorphological representativeness ----
geomorphologic_variables <-  raster::stack(list.files("T:/MODCLIM_R_DATA/ANALISIS/GEO/IP/2/", ".tif", full.names = T))
geomorphologic_variables <- projectRaster(geomorphologic_variables, crs = reference_system)

# Crop raster to study area
present_climatic_variables <-  raster::mask(crop(present_climatic_variables, study_area), study_area)
future_climatic_variables <-  raster::mask(crop(future_climatic_variables, study_area), study_area)
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
summary(PCA_geomorphologic_variables$model)


round(PCA_present_climatic_variables$model$loadings[,1:6],3)
round(PCA_future_climatic_variables$model$loadings[,1:6],3)
round(PCA_geomorphologic_variables$model$loadings[,1:5],3)

selected_PCA_present_climatic_variables <- PCA_present_climatic_variables$map[[1:4]]
selected_PCA_future_climatic_variables <- PCA_future_climatic_variables$map[[1:4]]
selected_PCA_geomorphological_variable <- PCA_geomorphologic_variables$map[[1:4]]

plot(selected_PCA_present_climatic_variables)


writeRaster(selected_PCA_present_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
writeRaster(selected_PCA_future_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
writeRaster(selected_PCA_geomorphological_variable, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")

PCA_present <-selected_PCA_present_climatic_variables
PCA_future <- selected_PCA_future_climatic_variables 
PCA_geo <- selected_PCA_geomorphological_variable 

###########################################################################################################
####################################################################################

PCA_present <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
PCA_future <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
PCA_geo <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")


PCA_geo <- resample(PCA_geo, PCA_present, method='bilinear')

delete <- c("PCA_present", "PCA_future", "PCA_geo", "study_area", "polygon")
rm(list=(ls()[ls()!= delete]))

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

data <- cbind(PCA_geo, PCA_present)

str(PCA_future$x)
data <- left_join(PCA_geo, PCA_present, by = c("x", "y"))
data <- left_join(data, PCA_future, by = c("x", "y"))
data <- na.omit(data)

colnames(data) <- c("a", "a", "a","a", "a", "a","a", "a", "a","a", "a", "a","a", "a", "a","a", "a")
colnames(data) <- c("x", "y", "PC_1", "PC_2", "PC_3", "PC_4", "Period")


data2 <- data[,1:7]
data2 <- rbind(data2, data[,c(1,2, 8,9,10,11,12)])
data2 <- rbind(data2, data[,c(1,2, 13,14,15,16, 17)])
colnames(data2) <- c("x", "y", "PC_1", "PC_2", "PC_3", "PC_4", "Period")

data <- data2
rm(data2)
#shp
he renombrado los archivos
PCA_present <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
PCA_future <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
PCA_geo <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")

raster <- raster::stack(PCA_present, PCA_future, PCA_geo)
#nrow(polygon)
names <- polygon$ORIG_NAME

mh_f <- data.frame(matrix(1,    # Create empty data frame
                          nrow = nrow(raster),
                          ncol = length(names)))

names(mh_f) <- names
beginCluster()
for (i in 1:nrow(polygon)){
  pol <- polygon[i,]
  raster_polygon <- raster::mask(raster::crop(raster, pol), pol)
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
