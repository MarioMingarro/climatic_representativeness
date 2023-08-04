library(raster)
library(sf)
library(tidyverse)
library(rgeos)
library(terra)

pre <- raster("D:/REPRESENATTIVENESS/TEST/mh_present_IPSL_2040_2070_SSP85Parque Nacional de Aiguestortes i Estany de Sant Maurici_T.tif")
fut <- raster("D:/REPRESENATTIVENESS/TEST/mh_future_IPSL_2040_2070_SSP85Parque Nacional de Aiguestortes i Estany de Sant Maurici_T.tif")

study_area <- read_sf("D:/REPRESENATTIVENESS/IP/IP.shp")
polygon <- read_sf("D:/REPRESENATTIVENESS/ENP/RAMSAR/Selected_Ramsar.shp")

reference_system <- projection("+init=epsg:25828") # "+proj=longlat +datum=WGS84 +no_defs"   


pre <- projectRaster(pre, crs = reference_system)
fut <- projectRaster(fut, crs = reference_system)
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))

polygon <- polygon[5,]

pre_in <- raster::mask(pre, polygon)
n_pre_in <- length(na.omit(pre_in@data@values))
pre_out <- raster::mask(pre, polygon, inverse = T)
fut_out <- raster::mask(fut, polygon, inverse = T)
pre_out <- rasterToPoints(pre_out)
fut_out <- rasterToPoints(fut_out)

class(pre_out)

kk <- SpatialPoints(pre_out[,1:2])

polygon <- as(polygon, 'Spatial')
pre_out <- SpatialPoints(pre_out[,1:2])
fut_out <- SpatialPoints(fut_out[,1:2])

distance_pre <- sort(rgeos::gDistance(pre_out, polygon, byid=TRUE))
distance_fut <- sort(rgeos::gDistance(fut_out, polygon, byid=TRUE))

mean(distance_pre[1:n_pre_in])
mean(distance_fut[1:n_pre_in])

######################○

reference_system <- projection("+init=epsg:25828") # "+proj=longlat +datum=WGS84 +no_defs"   


mh_present_umbral <- projectRaster(mh_present_umbral, crs = reference_system)
mh_future_umbral <- projectRaster(mh_future_umbral, crs = reference_system)
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))

isolation <- data.frame(name=character(), 
                        present_dist=numeric(), 
                        future_dist=numeric())



for (i in 1:nlayers(mh_present_umbral)){
  aa <- data.frame(name="a", 
                   present_dist=2, 
                   future_dist=3)
  pol <- polygon[i,]
  pre_in <- raster::mask(mh_present_umbral, pol)
  n_pre_in <- length(na.omit(pre_in@data@values))
  pre_out <- raster::mask(mh_present_umbral, pol, inverse = T)
  fut_out <- raster::mask(mh_future_umbral, pol, inverse = T)
  pre_out <- rasterToPoints(pre_out)
  fut_out <- rasterToPoints(fut_out)
  
  polygon <- as(polygon, 'Spatial')
  pre_out <- SpatialPoints(pre_out[,1:2])
  fut_out <- SpatialPoints(fut_out[,1:2])
  
  distance_pre <- sort(rgeos::gDistance(pre_out, pol, byid=TRUE))
  distance_fut <- sort(rgeos::gDistance(fut_out, pol, byid=TRUE))
  
  
  aa$name <- pol$NAME
  aa$present_dist <-  mean(distance_pre[1:n_pre_in])
  aa$future_dist <-  mean(distance_fut[1:n_pre_in])
  isolation <- rbind(isolation, aa)
}

i=1
