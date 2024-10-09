library(terra)
library(sf)
library(tidyverse)
library(caret)
library(tictoc)
library(spdep)
library(corrplot)

# library(RStoolbox)
# library(stringr)

gc(reset = T)

dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA_OLD/PRESENTE/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA_OLD/RCP85_2050/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_OLD/"


# Cargar área de estudio y polígonos
study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST/parques_nacionales.shp")

# Iniciar temporizador
tic()

# Cargar datos climáticos
present_climatic_variables <- terra::rast(list.files(dir_present_climate_data, "\\.tif$", full.names = TRUE))
future_climatic_variables <- terra::rast(list.files(dir_future_climate_data, "\\.tif$", full.names = TRUE))

# Asignar nombres a las variables
names(present_climatic_variables) <-  c("isotermal", "prec_meshum","rango_temp_anual", "temp_max_mescal", "temp_med_anual")# Excluye las que no están
names(future_climatic_variables) <- names(present_climatic_variables)

reference_system <-"EPSG:4326" 
present_climatic_variables <- terra::project(present_climatic_variables, reference_system)
future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))


# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)
# Raster to df
data_present_climatic_variables <- terra::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- terra::as.data.frame(future_climatic_variables, xy = TRUE)
# Delete NA
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)

# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:length(data_present_climatic_variables)])

# Select variables less correlated 
drop_1  <-  caret::findCorrelation(cor, cutoff = .8)
drop  <-  names(data_present_climatic_variables[,3:length(data_present_climatic_variables)])[drop_1]
data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]

present_climatic_variables <- terra::subset(present_climatic_variables, !names(present_climatic_variables) %in% drop)
future_climatic_variables <- terra::subset(future_climatic_variables, !names(future_climatic_variables) %in% drop)


#corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
#          method = "number",
#          type = "upper")
#corrplot(cor(data_future_climatic_variables[3:length(data_present_climatic_variables)]),
#          method = "number",
#          type = "upper")


# Add field period 
data_present_climatic_variables <- dplyr::mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- dplyr::mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)
data <- rbind(data_present_climatic_variables, data_future_climatic_variables)

# Create name object
names <- polygon$NAME

# Mahalanobis distance ----
# Create empty data frame
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

# MH raster
for(j in 4:length(mh)){
  mh_raster <- dplyr::filter(mh, Period == "Present")
  mh_raster <- terra::rast(mh_raster[, c(1:2, j)], crs = reference_system)
  names(mh_raster) <- colnames(mh[j])
  writeRaster(mh_raster, paste0(dir_result, "PRE_", names[j-3], ".tif"), overwrite=TRUE)
  plot(mh_raster)
  
  pol <- polygon[j - 3, ]
  puntos_todos <- terra::as.points(mh_raster)
  puntos_todos <- sf::st_as_sf(puntos_todos)
  colnames(puntos_todos) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos, pol)
  puntos_todos$inv_mh <-  1 / puntos_todos$mh
  # Autocorrelacion Local ----
  coords <- st_coordinates(puntos_todos)
  puntos_todos <- cbind(puntos_todos, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos$inv_mh , lw)
  
  LM_df <- as.data.frame(LM)
  puntos_todos$SA <- LM_df$Ii  # Índice local de Moran (Ii)
  puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))` # P-valor del índice local de Moran
  
  puntos_todos$lisa_cluster <- case_when(
    puntos_todos$SA_sig >= 0.05 ~ NA,
    # No significativo
    puntos_todos$mh < max(puntos_dentro$mh) &
      puntos_todos$SA > 0 ~ 1,
    # Punto alto rodeado de puntos altos
    puntos_todos$mh < max(puntos_dentro$mh) &
      puntos_todos$SA < 0 ~ NA,
    # Punto bajo rodeado de puntos bajos
    puntos_todos$mh > max(puntos_dentro$mh) &
      puntos_todos$SA > 0 ~ NA,
    # Punto alto rodeado de puntos bajos
    puntos_todos$mh > max(puntos_dentro$mh) &
      puntos_todos$SA < 0 ~ NA
  )
  
  puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)
  
  
  res <- res(mh_raster)  # Resolución del raster existente
  bbox <- ext(mh_raster)         # Extensión del raster existente
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])  # Número de filas
  ncols <- round((bbox[2] - bbox[1]) / res[1])  # Número de columnas
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_todos <- puntos_todos %>% mutate(aa = puntos_todos$inv_mh * puntos_todos$SA)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "PRE_sig_5_", names[j - 3], ".tif"), overwrite =
                TRUE)
  parches <- terra::patches(raster, directions = 8, zeroAsNA=T)
  crs(parches) <- crs(raster)
  parche_pt <- terra::extract(parches, puntos_todos)
  puntos_todos <- puntos_todos %>% 
    mutate(parche = parche_pt$patches)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "parche")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "PRE_patch_", names[j - 3], ".tif"), overwrite =
                TRUE)
  plot(raster)
  raster <- terra::rasterize(puntos_vect, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "PRE_SA_", names[j - 3], ".tif"), overwrite =
                TRUE)
}

for(j in 4:length(mh)){
  mh_raster <- dplyr::filter(mh, Period == "Future")
  mh_raster <- terra::rast(mh_raster[, c(1:2, j)], crs = reference_system)
  names(mh_raster) <- colnames(mh[j])
  writeRaster(mh_raster, paste0(dir_result, "RCP85_2050_", names[j-3], ".tif"), overwrite=TRUE)
  plot(mh_raster)
  
  pol <- polygon[j - 3, ]
  puntos_todos <- terra::as.points(mh_raster)
  puntos_todos <- sf::st_as_sf(puntos_todos)
  colnames(puntos_todos) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos, pol)
  puntos_todos$inv_mh <-  1 / puntos_todos$mh
  # Autocorrelacion Local ----
  coords <- st_coordinates(puntos_todos)
  puntos_todos <- cbind(puntos_todos, coords)
  kk <- puntos_todos[1:10, ]
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb,  zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos$inv_mh , lw)
  
  LM_df <- as.data.frame(LM)
  puntos_todos$SA <- LM_df$Ii  # Índice local de Moran (Ii)
  puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))` # P-valor del índice local de Moran
  
  puntos_todos$lisa_cluster <- case_when(
    puntos_todos$SA_sig >= 0.01 ~ NA, # No significativo
    puntos_todos$mh < max(puntos_dentro$mh) &
      puntos_todos$SA > 0 ~ 1,
    puntos_todos$mh < max(puntos_dentro$mh) &
      puntos_todos$SA < 0 ~ NA,
    puntos_todos$mh > max(puntos_dentro$mh) &
      puntos_todos$SA > 0 ~ NA,
    puntos_todos$mh > max(puntos_dentro$mh) &
      puntos_todos$SA < 0 ~ NA)
  
  puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)
  
  
  res <- res(mh_raster)  # Resolución del raster existente
  bbox <- ext(mh_raster)         # Extensión del raster existente
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])  # Número de filas
  ncols <- round((bbox[2] - bbox[1]) / res[1])  # Número de columnas
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_todos <- puntos_todos %>% mutate(aa = puntos_todos$inv_mh * puntos_todos$SA)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "RCP85_2050_sig_", names[j - 3], ".tif"), overwrite =
                TRUE)
  parches <- terra::patches(raster, directions = 8, zeroAsNA=T)
  crs(parches) <- crs(raster)
  parche_pt <- terra::extract(parches, puntos_todos)
  puntos_todos <- puntos_todos %>% 
    mutate(parche = parche_pt$patches)
  puntos_vect <- vect(puntos_todos)
  
  raster <- terra::rasterize(puntos_vect, raster_template, field = "parche")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "RCP85_2050_patch_", names[j - 3], ".tif"), overwrite =
                TRUE)
  plot(raster)
  raster <- terra::rasterize(puntos_vect, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster)
  writeRaster(raster, paste0(dir_result, "RCP85_2050_SA_", names[j - 3], ".tif"), overwrite =
                TRUE)
}

toc()
