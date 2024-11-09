library(terra)
library(sf)
library(tidyverse)
library(caret)
library(tictoc)
library(spdep)
library(corrplot)
library(doParallel)
library(foreach)


closeAllConnections()
gc(reset = T)

dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURO/GFDL/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_AREA_RED/RES/"

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_AREA_RED/MURCIA.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_AREA_RED/KBA_MURCIA_TEST2.shp")

# Create name object
names <- polygon$NatName

tic()

# CLIMATE ----
## Load data ----
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = T))
present_climatic_variables <- present_climatic_variables[[1:17]]

future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = T))
future_climatic_variables <- future_climatic_variables[[1:17]]

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

# Raster to df
data_present_climatic_variables <- terra::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- terra::as.data.frame(future_climatic_variables, xy = TRUE)

# Delete NA
data_present_climatic_variables <- na.omit(data_present_climatic_variables)
data_future_climatic_variables <- na.omit(data_future_climatic_variables)

# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:length(data_present_climatic_variables)])

# Select variables less correlated 
drop_1  <-  caret::findCorrelation(cor, cutoff = .8)
drop  <-  names(data_present_climatic_variables[,3:length(data_present_climatic_variables)])[drop_1]
data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]

present_climatic_variables <- terra::subset(present_climatic_variables, !names(present_climatic_variables) %in% drop)
future_climatic_variables <- terra::subset(future_climatic_variables, !names(future_climatic_variables) %in% drop)

# Add field period 
data_present_climatic_variables <- dplyr::mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- dplyr::mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)




# Función para procesar los datos de la serie presente
process_present <- function(j) {
  # Presente
  data <- data_present_climatic_variables
  
  # Mahalanobis distance calculation ----
  
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(raster::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data[,4:length(data)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data[,4:length(data)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh <- cbind(data[,c(1:3)], mh_p)
  
  mh_raster <- dplyr::filter(mh, Period == "Present")
  mh_raster <- terra::rast(mh_raster[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster) <- colnames(mh[j+3])
  writeRaster(mh_raster, paste0(dir_result, "PRE_", names[j], ".tif"), overwrite = TRUE)
  plot(mh_raster)
  
  pol <- polygon[j, ]
  puntos_todos <- terra::as.points(mh_raster)
  puntos_todos <- sf::st_as_sf(puntos_todos)
  colnames(puntos_todos) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos, pol)
  puntos_todos$inv_mh <- 1 / puntos_todos$mh
  max_mh_p <- max(puntos_dentro$mh)
  
  coords <- st_coordinates(puntos_todos)
  puntos_todos <- cbind(puntos_todos, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos$inv_mh, lw)
  LM_df <- as.data.frame(LM)
  
  puntos_todos$SA <- LM_df$Ii
  puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))`
  
  puntos_todos$lisa_cluster <- case_when(
    puntos_todos$SA_sig >= 0.01 ~ 0,
    puntos_todos$mh < max_mh_p & puntos_todos$SA > 0 ~ 1,
    puntos_todos$mh < max_mh_p & puntos_todos$SA < 0 ~ 0,
    puntos_todos$mh > max_mh_p & puntos_todos$SA > 0 ~ 0,
    puntos_todos$mh > max_mh_p & puntos_todos$SA < 0 ~ 0
  )
  
  puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)
  
  res <- res(mh_raster)
  bbox <- ext(mh_raster)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_todos <- puntos_todos %>% mutate(aa = puntos_todos$inv_mh * puntos_todos$SA)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "PRE_sig_", names[j], ".tif"), overwrite = TRUE)
  
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "PRE_SA_", names[j], ".tif"), overwrite = TRUE)
  
  # FUTURO 
  data <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  # Mahalanobis distance calculation ----
  
  mh_f <- data.frame(matrix(1,    
                            nrow = nrow(data),
                            ncol = length(names)))
  
  names(mh_f) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(raster::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data[,4:length(data)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data[,4:length(data)]), 
                      inverted = F)
    
    mh_f[,i] <- mh
  }
  
  mh <- cbind(data[,c(1:3)], mh_f)
  
  mh_raster <- dplyr::filter(mh, Period == "Future")
  mh_raster <- terra::rast(mh_raster[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster) <- colnames(mh[j+3])
  writeRaster(mh_raster, paste0(dir_result, "SSP585_2070_", names[j], ".tif"), overwrite = TRUE)
  
  pol <- polygon[j, ]
  puntos_todos <- terra::as.points(mh_raster)
  puntos_todos <- sf::st_as_sf(puntos_todos)
  colnames(puntos_todos) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos, pol)
  puntos_todos$inv_mh <- 1 / puntos_todos$mh
  
  coords <- st_coordinates(puntos_todos)
  puntos_todos <- cbind(puntos_todos, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos$inv_mh, lw)
  LM_df <- as.data.frame(LM)
  
  puntos_todos$SA <- LM_df$Ii
  puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))`
  
  puntos_todos$lisa_cluster <- case_when(
    puntos_todos$SA_sig >= 0.01 ~ 0,
    puntos_todos$mh < max_mh_p & puntos_todos$SA > 0 ~ 1,
    puntos_todos$mh < max_mh_p & puntos_todos$SA < 0 ~ 0,
    puntos_todos$mh > max_mh_p & puntos_todos$SA > 0 ~ 0,
    puntos_todos$mh > max_mh_p & puntos_todos$SA < 0 ~ 0
  )
  
  puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)
  
  res <- res(mh_raster)
  bbox <- ext(mh_raster)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_todos <- puntos_todos %>% mutate(aa = puntos_todos$inv_mh * puntos_todos$SA)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "SSP585_2070_sig_", names[j], ".tif"), overwrite = TRUE)
  
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "SSP585_2070_SA_", names[j], ".tif"), overwrite = TRUE)
}




# PRESENTE ----
closeAllConnections()
num_cores <- 5  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

foreach(j = 1:nrow(polygon), .packages = c("terra", "sf", "dplyr", "spdep")) %dopar% {
  process_present(j)
}


stopCluster(cl)
toc()








process_future <- function(j) {
  
}

























