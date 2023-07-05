library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)


PCA_presente <- raster::stack(list.files("T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Presente", "\\.TIF$", full.names = T))
PCA_futuro <- raster::stack(list.files("T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Futuro", "\\.TIF$", full.names = T))

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))
polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/national_parks_2.shp")



data_PCA_presente <- raster::as.data.frame(PCA_presente, xy = TRUE)
data_PCA_presente <- na.omit(data_PCA_presente)
data_PCA_presente <- mutate(data_PCA_presente, Periodo = c("Presente"))
data_PCA_presente <- mutate(data_PCA_presente, Var = c("climatic"))

data_PCA_futuro <- raster::as.data.frame(PCA_futuro, xy = TRUE)
data_PCA_futuro <- na.omit(data_PCA_futuro)
data_PCA_futuro <- mutate(data_PCA_futuro, Periodo = c("Futuro"))
data_PCA_futuro <- mutate(data_PCA_futuro, Var = c("climatic"))


colnames(data_PCA_futuro) <- colnames(data_PCA_presente)

data <- rbind(data_PCA_presente, data_PCA_futuro)

PCA_futuro <- resample(PCA_futuro, PCA_presente, method='bilinear')
raster <- raster::stack(PCA_presente, PCA_futuro)
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
  colnames(data_polygon) <- c("a","a","a","a","a","a","a","a","a","a", "a", "a")
  data_polygon_1 <- rbind(data_polygon[,c(1:7)], data_polygon[,c(1:2, 8,9,10,11,12)])
  data_polygon <- data_polygon_1
  colnames(data_polygon) <- colnames(data_PCA_presente[1:7])
  
  # Mahalanobis distance
  
  mh <- mahalanobis(data[,3:6], 
                    colMeans(data_polygon[,3:6]), 
                    cov(data[,3:6]), 
                    inverted = F)
  
  # Agregar información espacial al reusltado de mh
  mh_f[,i] <- mh
}

# Agregar informacion espacial al reusltado de mh
mh <- cbind(data[,c(1:2, 8,9)], mh_f)


#Creamos raster
mh_presente<- raster::brick()

for(j in 5:length(mh)){
  mh_f <- dplyr::filter(mh, Periodo == "Presente")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_presente <- raster::stack(mh_presente, mh_f)
}

mh_futuro <- raster::brick()

for(j in 5:length(mh)){
  mh_f <- dplyr::filter(mh, Periodo == "Futuro")
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
plot(mh_presente[[8]])
plot(mh_futuro_umbral)

writeRaster(mh_presente[[8]], "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_present_gua.tif")
writeRaster(mh_present_umbral[[8]], "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_present_umbral_gua.tif")
writeRaster(mh_futuro_umbral[[8]], "T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/mh_futuro_umbral_gua.tif")

plot(mh_present_umbral)
plot(mh_futuro_umbral)

# Selección de esas distancias en PI
# Presente ----
mh_presente_bin <- reclassify(mh_f, c(mh_poligono,Inf,NA))
plot(mh_presente_bin)
plot(polygon[1], add = T)




################################################################
library(readxl)
datos <- read_excel("T:/MODCLIM_R_DATA/datos_test_JML.xlsx")
dat <- datos[,5:9]
centroide <- filter(datos, datos$Centroide == 1)
centroide <- centroide[,5:9]

mh <- mahalanobis(dat, 
                  colMeans(centroide), 
                  cov(dat), 
                  inverted = F)
resultado <- cbind(datos, mh)

cor(resultado$mh, resultado$MH_Jorge)
