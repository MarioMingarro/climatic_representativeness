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
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_AREA_RED/KBA_MURCIA_TEST3.shp")

# Create name object
names <- polygon$NatName
year <- "2070"
model = "GFDL"
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
present_climatic_variables <-  terra::mask (terra::crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(terra::crop(future_climatic_variables,  study_area), study_area)

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
representativeness <- function(j) {
  # PRESENT
  # Presente
  data_p <- data_present_climatic_variables
  
  # Mahalanobis distance calculation ----
  
  mh_p <- data.frame(matrix(1,    
                            nrow = nrow(data_p),
                            ncol = length(names)))
  
  names(mh_p) <- names
  
  for (i in 1:nrow(polygon)){
    pol <- polygon[i,]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p[,4:length(data_p)], 
                      colMeans(data_polygon[,3:length(data_polygon)]), 
                      cov(data_p[,4:length(data_p)]), 
                      inverted = F)
    
    mh_p[,i] <- mh
  }
  
  mh_p <- cbind(data_p[,c(1:3)], mh_p)
  mh_raster_p <- dplyr::filter(mh_p, Period == "Present")
  mh_raster_p <- terra::rast(mh_raster_p[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_p) <- colnames(mh_p[j+3])
  writeRaster(mh_raster_p, paste0(dir_result, "PRE_", names[j], ".tif"), overwrite = TRUE)
  
  pol <- polygon[j, ]
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos_p, pol)
  puntos_todos_p$inv_mh <- 1 / puntos_todos_p$mh
  max_mh_p <- max(puntos_dentro$mh)
  coords <- st_coordinates(puntos_todos_p)
  puntos_todos_p <- cbind(puntos_todos_p, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos_p$inv_mh, lw)
  LM_df <- as.data.frame(LM)
  
  puntos_todos_p$SA <- LM_df$Ii
  puntos_todos_p$SA_sig <- LM_df$`Pr(z != E(Ii))`
  
  puntos_todos_p$lisa_cluster <- case_when(
    puntos_todos_p$SA_sig >= 0.01 ~ 0,
    puntos_todos_p$mh < max_mh_p & puntos_todos_p$SA > 0 ~ 1,
    puntos_todos_p$mh < max_mh_p & puntos_todos_p$SA < 0 ~ 0,
    puntos_todos_p$mh > max_mh_p & puntos_todos_p$SA > 0 ~ 0,
    puntos_todos_p$mh > max_mh_p & puntos_todos_p$SA < 0 ~ 0
  )
  
  puntos_todos_p$lisa_cluster <- as.numeric(puntos_todos_p$lisa_cluster)
  
  res <- res(mh_raster_p)
  bbox <- ext(mh_raster_p)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_p <- vect(puntos_todos_p)
  
  raster <- terra::rasterize(puntos_vect_p, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster_p)
  writeRaster(raster, paste0(dir_result, "PRE_sig_", names[j], ".tif"), overwrite = TRUE)
  
  raster <- terra::rasterize(puntos_vect_p, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster_p)
  writeRaster(raster, paste0(dir_result, "PRE_SA_", names[j], ".tif"), overwrite = TRUE)
  
  
  
  # FUTURO
  data_p_f <- rbind(data_present_climatic_variables, data_future_climatic_variables)
  
  mh_f <- data.frame(matrix(1,    
                            nrow = nrow(data_p_f),
                            ncol = length(names)))
  
  names(mh_f) <- names
  
  for (i in 1:nrow(polygon)) {
    pol <- polygon[i, ]
    raster_polygon <- terra::mask(terra::crop(present_climatic_variables, pol), pol)
    data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
    data_polygon <- na.omit(data_polygon)
    
    mh <- mahalanobis(data_p_f[, 4:length(data_p_f)], colMeans(data_polygon[, 3:length(data_polygon)]), cov(data_p_f[, 4:length(data_p_f)]), inverted = F)
    
    mh_f[, i] <- mh
  }
  
  mh_f <- cbind(data_p_f[, c(1:3)], mh_f)
  
  mh_raster_f <- dplyr::filter(mh_f, Period == "Future")
  mh_raster_f <- terra::rast(mh_raster_f[, c(1:2, j+3)], crs = reference_system)
  names(mh_raster_f) <- colnames(mh_f[j+3])
  writeRaster(mh_raster_f,
              paste0(dir_result, model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  pol <- polygon[j, ]
  puntos_todos_f <- terra::as.points(mh_raster_f)
  puntos_todos_f <- sf::st_as_sf(puntos_todos_f)
  colnames(puntos_todos_f) <- c("mh", "geometry")
  puntos_todos_f$inv_mh <- 1 / puntos_todos_f$mh
  
  coords <- st_coordinates(puntos_todos_f)
  puntos_todos_f <- cbind(puntos_todos_f, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos_f$inv_mh, lw)
  LM_df <- as.data.frame(LM)
  
  puntos_todos_f$SA <- LM_df$Ii
  puntos_todos_f$SA_sig <- LM_df$`Pr(z != E(Ii))`
  
  puntos_todos_f$lisa_cluster <- case_when(
    puntos_todos_f$SA_sig >= 0.01 ~ 0,
    puntos_todos_f$mh < max_mh_p &
      puntos_todos_f$SA > 0 ~ 1,
    puntos_todos_f$mh < max_mh_p &
      puntos_todos_f$SA < 0 ~ 0,
    puntos_todos_f$mh > max_mh_p &
      puntos_todos_f$SA > 0 ~ 0,
    puntos_todos_f$mh > max_mh_p &
      puntos_todos_f$SA < 0 ~ 0
  )
  
  puntos_todos_f$lisa_cluster <- as.numeric(puntos_todos_f$lisa_cluster)
  
  res <- res(mh_raster_f)
  bbox <- ext(mh_raster_f)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_f <- vect(puntos_todos_f)
  
  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster_f)
  writeRaster(raster,
              paste0(dir_result, "Sig", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)
  
  raster <- terra::rasterize(puntos_vect_f, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster_f)
  writeRaster(raster,
              paste0(dir_result, "SA", model, "_", year, "_", names[j], ".tif"),
              overwrite = TRUE)

}
  

for(j in 1:length(names)){
  representativeness(j)
}


# PRESENTE ----
closeAllConnections()
num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

foreach(j = 1:length(names),
        .packages = c("terra", "sf", "dplyr", "spdep")) %dopar% {
          representativeness(j)
        }


stopCluster(cl)
toc()




rm(list = setdiff(ls(), c("data_present_climatic_variables", "data_future_climatic_variables", "dir_future_climate_data",
                          "dir_present_climate_data", "dir_result", "future_climatic_variables", "max_mh_p", "polygon", 
                          "present_climatic_variables", "reference_system", "study_area", "names", "j", "model", "year"
)))  


rm(list = setdiff(ls(), c("dir_future_climate_data","dir_present_climate_data", "dir_result", "future_climatic_variables", "polygon", 
                          "present_climatic_variables", "reference_system", "study_area", "names", "model", "year", "data_present_climatic_variables",
                          "data_future_climatic_variables", "j" , "max_mh_p"
)))






