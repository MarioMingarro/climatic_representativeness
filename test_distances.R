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

plot(polygon)

polygon <- polygon[5,]

plot(polygon)

aa <- raster::mask(pre, polygon)
kk <- length(na.omit(aa@data@values))
pre_out <- raster::mask(pre, polygon, inverse = T)
fut_out <- raster::mask(fut, polygon, inverse = T)
pre_out <- rasterToPoints(pre_out)
fut_out <- rasterToPoints(fut_out)

plot(pre_out)
plot(fut_out, col = "red", add= T)

bb <- gDistance(pre_out, polygon,  byid=T)
