library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/national_parks_2.shp")


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
present_climatic_variables <-  raster::mask(crop(present_climatic_variables, geomorphologic_variables[[1]]), geomorphologic_variables[[1]])
future_climatic_variables <-  raster::mask(crop(future_climatic_variables, geomorphologic_variables[[1]]), geomorphologic_variables[[1]])
geomorphologic_variables <-  mask(crop(geomorphologic_variables, present_climatic_variables[[1]]), present_climatic_variables[[1]])


variables <- raster::stack(present_climatic_variables, future_climatic_variables, geomorphologic_variables)
cor <-layerStats(variables,'pearson')

# Climate change variables
#climatic_change_variables <-

# PCA ----
# Climatic representativeness
PCA_variables <- RStoolbox::rasterPCA(variables, spca = TRUE)


# PCA results
summary(PCA_variables$model)


round(PCA_variables$model$loadings[,1:6],3)


selected_PCA_variables <- PCA_variables$map[[1:6]]


plot(selected_PCA_variables)

data_PCA <- raster::as.data.frame(variables, xy = TRUE)
data_PCA <- na.omit(data_PCA)

presente <- data_PCA[,1:21]
futuro <- data_PCA[,c(1:2, 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)]

geo <- data_PCA[,c(1:2, 41,42,43,44,45,46)]


presente <- mutate(presente, Periodo = c("Presente"))
futuro <- mutate(futuro, Periodo = c("Futuro"))

presente <- left_join(presente, geo, by = c("x", "y"))
futuro <- left_join(futuro, geo, by = c("x", "y"))

colnames(futuro) <- colnames(presente)
data <- rbind(presente, futuro)




writeRaster(selected_PCA_present_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")




###########################################################################################################
####################################################################################

PCA_present <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
PCA_future <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
PCA_geo <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")


#PCA_geo <- resample(PCA_geo, PCA_present, method='bilinear')
#
#delete <- c("PCA_present", "PCA_future", "PCA_geo", "study_area", "polygon")
#rm(list=(ls()[ls()!= delete]))

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/national_parks.shp")



data_PCA <- raster::as.data.frame(selected_PCA_variables, xy = TRUE)
data_PCA_present <- na.omit(data_PCA_present)
data_PCA_present <- mutate(data_PCA_present, Periodo = c("Presente"))
data_PCA_present <- mutate(data_PCA_present, Var = c("climatic"))

data_PCA_future <- raster::as.data.frame(PCA_future, xy = TRUE)
data_PCA_future <- na.omit(data_PCA_future)
data_PCA_future <- mutate(data_PCA_future, Periodo = c("Futuro"))
data_PCA_future <- mutate(data_PCA_future, Var = c("climatic"))

data_PCA_geo <- raster::as.data.frame(PCA_geo, xy = TRUE)
data_PCA_geo <- na.omit(data_PCA_geo)
data_PCA_geo <- mutate(data_PCA_geo, Periodo = c("Presente"))
data_PCA_geo <- mutate(data_PCA_geo, Var = c("geomor"))


colnames(data_PCA_present) <- c("x", "y", "PC_1", "PC_2", "PC_3", "PC_4", "Period", "Var")
colnames(data_PCA_geo) <- colnames(data_PCA_present)
colnames(data_PCA_future) <- colnames(data_PCA_present)


data <- rbind(data_PCA_present, data_PCA_geo)
data <- rbind(data, data_PCA_future)




raster <- raster::stack(PCA_present, PCA_future, PCA_geo)
#nrow(polygon)
names <- polygon$ORIG_NAME


mh_f <- data.frame(matrix(1,    # Create empty data frame
                          nrow = nrow(data),
                          ncol = length(names)))

names(mh_f) <- names

for (i in 1:nrow(polygon)){
  pol <- polygon[i,]
  raster_polygon <- raster::mask(raster::crop(raster, pol), pol)
  data_polygon <- raster::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  colnames(data_polygon) <- c("a", "a", "a","a", "a", "a","a", "a", "a","a", "a", "a","a", "a")
  data_polygon_1 <- rbind(data_polygon[,c(1:6)], data_polygon[,c(1:2, 7,8,9,10)])
  data_polygon_1 <- rbind(data_polygon_1, data_polygon[,c(1:2, 11,12,13,14)])
  data_polygon <- data_polygon_1
  colnames(data_polygon) <- c("x", "y", "PC_1", "PC_2", "PC_3", "PC_4")
  
  # Mahalanobis distance
  
  mh <- mahalanobis(data[,3:6], 
                    colMeans(data_polygon[,3:6]), 
                    cov(data[,3:6]), 
                    inverted = F)
  
  # Agregar información espacial al reusltado de mh
  mh_f[,i] <- mh
}

# Agregar informacion espacial al reusltado de mh
mh <- cbind(data[,c(1:2, 7,8)], mh_f)


#Creamos raster
mh_presente<- raster::brick()

for(j in 5:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Presente")
  mh_f <- dplyr::filter(mh_f, mh_f$Var == "climatic")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_presente <- raster::stack(mh_presente, mh_f)
}

mh_futuro <- raster::brick()

for(j in 5:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Futuro")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_futuro <- raster::stack(mh_futuro, mh_f)
}

writeRaster(mh_presente, "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_presente_NP.tif")
writeRaster(mh_futuro, "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_future_NP.tif")
plot(mh_futuro[[1]])

endCluster()


plot(mh_presente)
nlayers(mh_presente)

# Umbral para establecer el corte percentil 90
mh_present_umbral <- raster::brick()
mh_futuro_umbral <- raster::brick()

for (i in 1:nlayers(mh_presente)){
  mh_poligono <- mask(crop(mh_presente[[i]], polygon[i,]), polygon[i,])
  mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
  mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.90))
  mh_presente_bin <- reclassify(mh_presente[[i]], c(mh_poligono,Inf,NA))
  mh_present_umbral <- raster::stack(mh_present_umbral, mh_presente_bin)
  mh_futuro_bin <- reclassify(mh_futuro[[i]], c(mh_poligono,Inf,NA))
  mh_futuro_umbral <- raster::stack(mh_futuro_umbral, mh_futuro_bin)
}
plot(mh_present_umbral)
plot(mh_futuro_umbral)
mh_futuro_umbral <- raster::brick()

for (i in 1:nlayers(mh_futuro_umbral)){
  mh_poligono <- mask(crop(mh_futuro[[i]], polygon[i,]), polygon[i,])
  mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
  mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.90))
  mh_futuro_bin <- reclassify(mh_futuro_umbral[[i]], c(mh_poligono,Inf,NA))
  mh_futuro_umbral <- raster::stack(mh_futuro_umbral, mh_futuro_bin)
}

writeRaster(mh_present_umbral, "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_present_umbral_NP.tif")
writeRaster(mh_futuro_umbral, "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_futuro_umbral_NP.tif")

plot(mh_present_umbral)
plot(mh_futuro_umbral)

# Selección de esas distancias en PI
# Presente ----
mh_presente_bin <- reclassify(mh_f, c(mh_poligono,Inf,NA))
plot(mh_presente_bin)
plot(polygon[1], add = T)
