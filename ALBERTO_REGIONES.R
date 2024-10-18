library(terra)
library(sf)
library(tidyverse)
library(stringr)


regiones <- read_sf("A:/REGIONES_KARST/geokarst_reg_final_def.shp")

# Área total ----
regiones$area_m2 <- st_area(regiones)

# Latitud y longitud máxima y mínima ----
bbox_coords <- function(regiones) {
  bbox <- st_bbox(regiones)
  return(data.frame(
    lon_min = bbox["xmin"],
    lon_max = bbox["xmax"],
    lat_min = bbox["ymin"],
    lat_max = bbox["ymax"]
  ))
}

coords_df <- do.call(rbind, lapply(st_geometry(regiones), bbox_coords))
regiones <- cbind(regiones, coords_df)

# Temperatura máxima/mínima/media ----

dir_TMAX<- "B:/A_DATA/TERRACLIMATE/TMAX/"
dir_TMIN<- "B:/A_DATA/TERRACLIMATE/TMIN/"
dir_TMED<- "B:/A_DATA/TERRACLIMATE/TMED/"
dir_PCP<- "B:/A_DATA/TERRACLIMATE/PCP/"

kk <- terra::rast("B:/A_DATA/TERRACLIMATE/agg_terraclimate_tmax_1958_CurrentYear_GLOBE.nc")

library(terra)
library(sf)
mask <- st_read("A:/REGIONES_KARST/AE.shp")
mask <- st_transform(regiones, crs("+proj=longlat +datum=WGS84 +no_defs"))
for (i in 1958:2023){
  terraclimate <- raster::brick(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data/TerraClimate_ppt_",i,".nc"))
  terraclimate <- terra::mask(terra::crop(terraclimate,mask), mask)
  for(j in 1:12){
    name <- names(terraclimate[[j]])
    name <- str_sub(str_sub(gsub("\\.", "_", name),  2, -1),  1, 7)
    writeRaster(terraclimate[[j]], paste0(dir_PCP, "Tmax_", name, ".tif"))
    }}

ff <- kk[[1]]

plot(ff)


TMAX <- terra::rast(list.files(dir_TMAX, "\\.tif$", full.names = T))
TMIN <- terra::rast(list.files(dir_TMIN, "\\.tif$", full.names = T))
TMED <- terra::rast(list.files(dir_TMED, "\\.tif$", full.names = T))
PCP <- terra::rast(list.files(dir_PCP, "\\.tif$", full.names = T))


names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")

names(future_climatic_variables) <- names(present_climatic_variables)

# Reference system ----
reference_system <- "EPSG:4326" 
terra::crs(present_climatic_variables)
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))

# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)
