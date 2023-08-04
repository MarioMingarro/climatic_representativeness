library(raster)
library(elevatr)
library(sf)

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Functions.R")
#--------------------

dem <- raster("T:/MODCLIM_R_DATA/ANALISIS/GEO/DEM_SRTM30.tif")

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))


mdt <- get_elev_raster(study_area, z = 7)
res(mdt)

dem <-  raster::mask(raster::crop(dem, study_area), study_area)

# Rescalar
dem <- resample(mdt, dem) 
#sistema ref
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
dem = projectRaster(dem, crs = reference_system)


geomorphologic <- terrain(dem, opt=c("slope", "aspect", "TPI", "TRI", "roughness", "flowdir"), unit='degrees')

plot(geomorphologic)

tri <-        geomorphologic[[1]]
tpi <-        geomorphologic[[2]]
roughness <-  geomorphologic[[3]]
slope <-      geomorphologic[[4]] 
aspect <-     geomorphologic[[5]] 
flowdir <-    geomorphologic[[6]] 


writeRaster(geomorphologic[[1]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/tri.tif")
writeRaster(geomorphologic[[2]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/tpi.tif")
writeRaster(geomorphologic[[3]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/roughness.tif")
writeRaster(geomorphologic[[4]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/slope.tif")
writeRaster(geomorphologic[[5]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/aspect.tif")
writeRaster(geomorphologic[[6]], "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/flowdir.tif")

unstack(geomorphologic)

# SOIL
library(geodata)
ph <- soil_world_vsi(var="phh2o", depth=5)
ph <- soil_world_vsi(var="phh2o", depth=5)
ph <- soil_world_vsi(var="phh2o", depth=5)
writeRaster(ph, "T:/MODCLIM_R_DATA/ANALISIS/GEO/VARIABLES/pH.tif")



geomorphologic <- raster::stack(list.files("T:/MODCLIM_R_DATA/ANALISIS/GEO/IP/", "\\.tif$", full.names = T))

rm(PCA)
geomorphologic_PCA <- RStoolbox::rasterPCA(geomorphologic, spca = TRUE)
summary(PCA$model)
PCA$model
round(PCA$model$loadings[,1:5],5)

writeRaster(geomorphologic, "T:/MODCLIM_R_DATA/analysis/geomorphologic.tif")

kk <- raster(geomorphologic, "T:/MODCLIM_R_DATA/analysis/bioclim_ordenado.tif")
nlayers(kk)
plot(kk)
