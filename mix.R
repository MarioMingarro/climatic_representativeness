library(raster)
library(sf)
library(stars)
library(tidyverse)
library(tictoc)
library(geodata)
library(stringr)

geodata::worldclim_country()


ae <- shapefile(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))
lat_long <- coordinates(spTransform(ae, CRS("+proj=longlat + datum=WGS84")))
lat<- 40.41 
long <- -3.8

rm(raster_presente)

tmin <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.tmn.dat.nc"), varname = "tmn")
tmax <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.tmx.dat.nc"), varname = "tmx")
pcp <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.pre.dat.nc"), varname = "pre")


nlayers(tmin)
plot(tmin)

raster_presente <- getData('worldclim', var='bio', res=0.5, lon=long, lat=lat)
plot(raster_presente)

kk <- area_estudio %>% 
  st_transform(., "+init=epsg:4326") %>%
  st_centroid()

long <- st_bbox(kk)[1]
lat <- st_bbox(kk)[2]

kk  %>% 
  mutate(lon = map_dbl(geometry, ~st_point_on_surface(.kk)[[1]]),
         lat = map_dbl(geometry, ~st_point_on_surface(.x)[[2]]))

raster_futuro <- getData('CMIP5', var='bio', res=0.5, rcp=85, model='AC', year=70, lon=long, lat=lat)

P_url<- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip"
download.file(P_url, destfile="wc2.1_2.5m_bio.zip")

tic()
# Presente -----
## Cargar los datos -----

raster_presente <- raster::stack(list.files("T:/MODCLIM_R_DATA/Presente", ".TIF", full.names = T))
raster_futuro <- raster::stack(list.files("T:/MODCLIM_R_DATA/Futuro", ".TIF", full.names = T))

poligono <- read_sf(list.files("T:/MODCLIM_R_DATA/shp", "\\.shp$", full.names = T))
area_estudio <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))

#Ajustar los datos espaciales (mismo sistema coordenadas)

if(projection(raster_presente) == projection(raster_futuro)){
  print("Adecuado")
} else {
  print("Se procede a la transformacion del SHP")
  raster_futuro <- projectRaster(raster_futuro, crs = projection(raster_presente))
  print(paste0("Nuevo SR: ", projection(raster_futuro)))
}

if(st_crs(poligono)$proj4string == projection(raster_presente)){
  print("Adecuado")
} else {
  print("Se procede a la transformacion del SHP")
  poligono <- st_transform(poligono, crs(projection(raster_presente)))
  print(paste0("Nuevo SR: ", projection(raster_presente)))
}

if(st_crs(area_estudio)$proj4string == projection(raster_presente)){
  print("Adecuado")
} else {
  print("Se procede a la transformacion del SHP")
  area_estudio <- st_transform(area_estudio, crs(projection(raster_presente)))
  print(paste0("Nuevo SR: ", projection(raster_presente)))
}

#Extraer raster en area de estudio
raster_presente <-  mask(crop(raster_presente, area_estudio), area_estudio)
raster_futuro <- mask(crop(raster_futuro, area_estudio), area_estudio)
raster_poligono <- mask(crop(raster_presente, poligono), poligono)

# Extraer datos del raster
# Presente
data_presente <- raster::as.data.frame(x = raster_presente, xy = TRUE)
data_presente <- na.omit(data_presente)
data_presente <- mutate(data_presente, Periodo = c("Presente"))

# Futuro 
data_futuro <- raster::as.data.frame(x = raster_futuro, xy = TRUE)
data_futuro <- na.omit(data_futuro)
data_futuro <- mutate(data_futuro, Periodo = c("Futuro"))

# Poligono

data_poligono <- raster::as.data.frame(raster_poligono, xy = TRUE)
data_poligono <- na.omit(data_poligono)

# Preparar datos mahalanobis

colnames(data_futuro) <- colnames(data_presente)
data_presente_futuro <- rbind(data_presente, data_futuro)


#Calculo de mahalanobis

mh <- mahalanobis(data_presente_futuro[,3:7], colMeans(data_poligono[,3:7]), cov(data_presente_futuro[,3:7]), inverted = F)

# Agregar información espacial al reusltado de mh
mh <- cbind(data_presente_futuro, mh)

mh_presente <- dplyr::filter(mh, Periodo == "Presente")
mh_futuro <- dplyr::filter(mh, Periodo == "Futuro")

mh_presente <- mh_presente[,c(1:2,9)]
mh_futuro <- mh_futuro[,c(1:2,9)]

#Creamos raster
mh_presente <- rasterFromXYZ(mh_presente)
mh_futuro <- rasterFromXYZ(mh_futuro)



plot(mh_presente)
plot(mh_futuro)


# Umbral para establecer el corte percentil 90
mh_poligono <- mask(crop(mh_presente, poligono), poligono)
mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.90))
mh_poligono


# Selección de esas distancias en PI
# Presente ----
mh_presente_bin <- reclassify(mh_presente, c(mh_poligono,Inf,0))
plot(mh_presente_bin)

# Futuro ----
mh_futuro_bin <- reclassify(mh_futuro, c(mh_poligono, Inf,0))
plot(mh_futuro_bin)

toc()


devtools::install_github("achubaty/grainscape")
