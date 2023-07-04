library(raster)
library(tidyverse)
library(sf)
library(dismo)

# Load raster data
tmin <- raster::stack(list.files("T:/MODCLIM_R_DATA/PRESENT/TMIN/", pattern = ".tif", full.names = TRUE))
tmax <- raster::stack(list.files("T:/MODCLIM_R_DATA/PRESENT/TMAX/", pattern = ".tif", full.names = TRUE))
pcp <-  raster::stack(list.files("T:/MODCLIM_R_DATA/PRESENT/PCP/", pattern = ".tif", full.names = TRUE))

paste0("The raster resolution is", res(tmin))

# Load study area shp
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))

# Reference system
projection(tmin)
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_transform(study_area, crs(reference_system))
rm(reference_system)

# Crop raster to study area
#Skip if var were obtained whit function
beginCluster()

tmin <- mask(crop(tmin, study_area), study_area)
tmax <- mask(crop(tmax, study_area), study_area)
pcp <-  mask(crop(pcp, study_area), study_area)

endCluster()


# Period dates
starting_year <- 1976
finishing_year <- 2016

year <- seq(starting_year, finishing_year, 1)

rm(starting_year)
rm(finishing_year)



### Creating 19 bioclimatic variables (Worldclim)----

beginCluster()

bioclim_all <- raster::stack()

for (i in 1:length(year)){
  monthly_tmax <- raster::subset(tmax, 
                                 grep(paste0(year[i]), 
                                      names(tmax), 
                                      value = T))
  monthly_tmax <- stack(monthly_tmax[[1]],
                 monthly_tmax[[5]],
                 monthly_tmax[[6]],
                 monthly_tmax[[7]],
                 monthly_tmax[[8]],
                 monthly_tmax[[9]], 
                 monthly_tmax[[10]],
                 monthly_tmax[[11]],
                 monthly_tmax[[12]],
                 monthly_tmax[[2]],
                 monthly_tmax[[3]],
                 monthly_tmax[[4]])
  
  monthly_tmin <- raster::subset(tmin, 
                                 grep(paste0(year[i]), 
                                      names(tmin), 
                                      value = T))
  monthly_tmin <- stack(monthly_tmin[[1]],
                        monthly_tmin[[5]],
                        monthly_tmin[[6]],
                        monthly_tmin[[7]],
                        monthly_tmin[[8]],
                        monthly_tmin[[9]], 
                        monthly_tmin[[10]],
                        monthly_tmin[[11]],
                        monthly_tmin[[12]],
                        monthly_tmin[[2]],
                        monthly_tmin[[3]],
                        monthly_tmin[[4]])
  
  monthly_pcp <- raster::subset(pcp, 
                                 grep(paste0(year[i]), 
                                      names(pcp), 
                                      value = T))
  monthly_pcp  <- stack(monthly_pcp[[1]],
                        monthly_pcp[[5]],
                        monthly_pcp[[6]],
                        monthly_pcp[[7]],
                        monthly_pcp[[8]],
                        monthly_pcp[[9]], 
                        monthly_pcp[[10]],
                        monthly_pcp[[11]],
                        monthly_pcp[[12]],
                        monthly_pcp[[2]],
                        monthly_pcp[[3]],
                        monthly_pcp[[4]])
  bioclim <- dismo::biovars(monthly_pcp, monthly_tmin, monthly_tmax) 
  bioclim_all <- raster::stack(bioclim_all, bioclim)
}
endCluster()

writeRaster(bioclim, "T:/MODCLIM_R_DATA/bioclim_ordenado.tif" )

plot(raster)
plot(bioclim)
nlayers(bioclim_all)

plot(monthly_pcp[[1:12]])


###
library(doParallel)
registerDoParallel(cl <- makeCluster(4))

bioclim_all <- raster::stack()

foreach::foreach(i = 1:2) %dopar% {
  monthly_tmax <- raster::subset(tmax, 
                                 grep(paste0(year[i]), 
                                      names(tmax), 
                                      value = T))
  monthly_tmax <- stack(monthly_tmax[[1]],
                        monthly_tmax[[5]],
                        monthly_tmax[[6]],
                        monthly_tmax[[7]],
                        monthly_tmax[[8]],
                        monthly_tmax[[9]], 
                        monthly_tmax[[10]],
                        monthly_tmax[[11]],
                        monthly_tmax[[12]],
                        monthly_tmax[[2]],
                        monthly_tmax[[3]],
                        monthly_tmax[[4]])
  
  monthly_tmin <- raster::subset(tmin, 
                                 grep(paste0(year[i]), 
                                      names(tmin), 
                                      value = T))
  monthly_tmin <- stack(monthly_tmin[[1]],
                        monthly_tmin[[5]],
                        monthly_tmin[[6]],
                        monthly_tmin[[7]],
                        monthly_tmin[[8]],
                        monthly_tmin[[9]], 
                        monthly_tmin[[10]],
                        monthly_tmin[[11]],
                        monthly_tmin[[12]],
                        monthly_tmin[[2]],
                        monthly_tmin[[3]],
                        monthly_tmin[[4]])
  
  monthly_pcp <- raster::subset(pcp, 
                                grep(paste0(year[i]), 
                                     names(pcp), 
                                     value = T))
  monthly_pcp  <- stack(monthly_pcp[[1]],
                        monthly_pcp[[5]],
                        monthly_pcp[[6]],
                        monthly_pcp[[7]],
                        monthly_pcp[[8]],
                        monthly_pcp[[9]], 
                        monthly_pcp[[10]],
                        monthly_pcp[[11]],
                        monthly_pcp[[12]],
                        monthly_pcp[[2]],
                        monthly_pcp[[3]],
                        monthly_pcp[[4]])
  bioclim <- dismo::biovars(monthly_pcp, monthly_tmin, monthly_tmax) 
  bioclim_all <- raster::stack(bioclim_all, bioclim)
}

stopImplicitCluster()
stopCluster(cl)

######################################################
# Climate velocity
