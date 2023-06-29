library(raster)
library(elevatr)

dem <- raster("T:/MODCLIM_R_DATA/DEM/MDT_PI_100M.tif")


mdt <- get_elev_raster(study_area, z = 7)
res(mdt)

mdt <-  mask(crop(mdt, study_area), study_area)

mdt <- resample(mdt, tmin)

reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
dem = projectRaster(dem, crs = reference_system)

geomorphologic <- terrain(mdt, opt=c("slope", "aspect", "TPI", "TRI", "roughness", "flowdir"), unit='degrees')
plot(geomorphologic)

tri <-        geomorphologic[[1]]
tpi <-        geomorphologic[[2]]
roughness <-  geomorphologic[[3]]
slope <-      geomorphologic[[4]] 
aspect <-     geomorphologic[[5]] 
flowdir <-    geomorphologic[[6]] 


writeRaster(geomorphologic[[1]], "T:/MODCLIM_R_DATA/analysis/geomophological/tri.tif")
writeRaster(geomorphologic[[2]], "T:/MODCLIM_R_DATA/analysis/geomophological/tpi.tif")
writeRaster(geomorphologic[[3]], "T:/MODCLIM_R_DATA/analysis/geomophological/roughness.tif")
writeRaster(geomorphologic[[4]], "T:/MODCLIM_R_DATA/analysis/geomophological/slope.tif")
writeRaster(geomorphologic[[5]], "T:/MODCLIM_R_DATA/analysis/geomophological/aspect.tif")
writeRaster(geomorphologic[[6]], "T:/MODCLIM_R_DATA/analysis/geomophological/flowdir.tif")

unstack(geomorphologic)

PCA <- RStoolbox::rasterPCA(geomorphologic, spca = TRUE)
summary(PCA$model)
PCA$model
round(PCA$model$loadings[,1:5],3)

writeRaster(geomorphologic, "T:/MODCLIM_R_DATA/analysis/geomorphologic.tif")

kk <- raster(geomorphologic, "T:/MODCLIM_R_DATA/analysis/bioclim_ordenado.tif")
nlayers(kk)
plot(kk)
