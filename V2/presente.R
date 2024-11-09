library(terra)
library(sf)
library(tidyverse)
library(tictoc)

library(caret)

library(RStoolbox)
library(stringr)
library(corrplot)

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
dir_present_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/"
dir_future_climate_data <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURO/IPSL/"
dir_result <- "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_PEDRIZA/TEST/"

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")
polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/TEST_PEDRIZA/KBA_pedriza.shp")


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

reference_system <-"EPSG:4326" 
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
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)




vif_filter <- function(x, th = 10) {
  
  # Función para calcular el VIF de todas las variables
  calc_vif <- function(df) {
    vif_values <- sapply(1:ncol(df), function(i) {
      formula <- as.formula(paste(names(df)[i], "~ ."))  # Formula de lm con todas las demás variables
      model <- lm(formula, data = df)
      return(1 / (1 - summary(model)$r.squared))
    })
    names(vif_values) <- colnames(df)
    return(vif_values)
  }
  
  # Convertir el objeto Raster a un data frame
  x <- as.data.frame(x, na.rm = TRUE)
  
  # Eliminar variables multicolineales iterativamente
  exc <- character(0)  # Lista de variables excluidas
  while (TRUE) {
    v <- calc_vif(x)  # Calcular el VIF
    if (max(v) < th) break  # Terminar si no hay VIF por encima del umbral
    ex <- names(v)[which.max(v)]  # Variable con mayor VIF
    exc <- c(exc, ex)  # Agregar a la lista de excluidos
    x <- x[, !(colnames(x) %in% ex)]  # Remover variable
  }
  
  # Crear lista de resultados
  result <- list(
    variables = colnames(x),
    excluded = exc,
    corMatrix = cor(x, method = "pearson"),
    results = data.frame(Variables = names(v), VIF = v)
  )
  
  return(result)
}


kk <- vif_filter(present_climatic_variables)
drop <- kk$excluded

data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]

present_climatic_variables <- terra::subset(present_climatic_variables, !names(present_climatic_variables) %in% drop)
future_climatic_variables <- terra::subset(future_climatic_variables, !names(future_climatic_variables) %in% drop)


corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number",
         type = "upper")
corrplot(cor(data_future_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number",
         type = "upper")



# Add field period 
data_present_climatic_variables <- dplyr::mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- dplyr::mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)
data <- rbind(data_present_climatic_variables, data_future_climatic_variables)


# Create name object
names <- polygon$NatName


# Mahalanobis distance ----

# Función para procesar los datos de la serie presente
representativeness <- function(j, th = .95) {
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
  plot(mh_raster_p)
  #writeRaster(mh_raster_p, paste0(dir_present, "PRE_", names[j], ".tif"), overwrite = TRUE)
  
  pol <- polygon[j, ]
  puntos_todos_p <- terra::as.points(mh_raster_p)
  puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
  colnames(puntos_todos_p) <- c("mh", "geometry")
  puntos_dentro <- sf::st_intersection(puntos_todos_p, pol)
  puntos_todos_p$inv_mh <- 1 / puntos_todos_p$mh
  
  th_mh_p <- quantile(na.omit(puntos_dentro$mh), probs = th)
  
  coords <- st_coordinates(puntos_todos_p)
  puntos_todos_p <- cbind(puntos_todos_p, coords)
  nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  
  LM <- spdep::localmoran(puntos_todos_p$inv_mh, lw)
  LM_df <- as.data.frame(LM)
  
  puntos_todos_p$SA <- LM_df$Ii
  puntos_todos_p$SA_sig <- LM_df$`Pr(z != E(Ii))`
  
  puntos_todos_p$th <- case_when(
    puntos_todos_p$mh > th_mh_p ~ 0,
    puntos_todos_p$mh <= th_mh_p  ~ 1
  )
  
  puntos_todos_p$th <- as.numeric(puntos_todos_p$th)
  
  res <- res(mh_raster_p)
  bbox <- ext(mh_raster_p)
  
  nrows <- round((bbox[4] - bbox[3]) / res[2])
  ncols <- round((bbox[2] - bbox[1]) / res[1])
  raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
  puntos_vect_p <- vect(puntos_todos_p)
  
  raster <- terra::rasterize(puntos_vect_p, raster_template, field = "th")
  crs(raster) <- crs(mh_raster_p)
  plot(raster)
  #writeRaster(raster, paste0(dir_present, "Pre_sig_", names[j], ".tif"), overwrite = TRUE)
  
  raster <- terra::rasterize(puntos_vect_p, raster_template, field = "SA")
  crs(raster) <- crs(mh_raster_p)
  #writeRaster(raster, paste0(dir_present, "Pre_SA_", names[j], ".tif"), overwrite = TRUE)
  
}

tic()
for(j in 1:length(names)){
  representativeness(j)
}
toc()



rm(list = setdiff(ls(), c("data_present_climatic_variables", "data_future_climatic_variables", "dir_future_climate_data",
                          "dir_present_climate_data", "dir_result", "future_climatic_variables", "polygon", 
                          "present_climatic_variables", "reference_system", "study_area", "names", "model", "year"
)))  



