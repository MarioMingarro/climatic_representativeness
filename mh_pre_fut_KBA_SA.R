library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)
library(corrplot)
library(caret)
library(tictoc)
gc(reset = T)


tic()
# CLIMATE ----
## Load data ----
## Climatic representativeness -----
present_climatic_variables <- raster::stack(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENTE/", "\\.tif$", full.names = T))

future_climatic_variables <- raster::stack(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/RCP85_2050", "\\.tif$", full.names = T))

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")

polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/KBA/KBA_test.shp")


# Reference system
reference_system <- projection(study_area) # "+proj=longlat +datum=WGS84 +no_defs"

if(compareCRS(present_climatic_variables, future_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  future_climatic_variables <- projectRaster(future_climatic_variables, crs = reference_system)
}

if(compareCRS(study_area, present_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  study_area <- st_transform(study_area, crs(reference_system))
}

if(compareCRS(polygon, present_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  polygon <- st_transform(polygon, crs(reference_system))
  print(paste0("Reference System was modified to ", reference_system))
}


# Crop raster to study area
present_climatic_variables <-  raster::mask(crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  raster::mask(crop(future_climatic_variables,  study_area), study_area)

# Extract raster data
data_present_climatic_variables <- raster::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- raster::as.data.frame(future_climatic_variables, xy = TRUE)

# Delete NA
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)

# Rename columns
colnames(data_present_climatic_variables) <- c("x","y","isotermal", "prec_meshum","rango_temp_anual", "temp_max_mescal", "temp_med_anual")

colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)



# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:7])

# Select variables less correlated 
drop_1  <-  findCorrelation(cor, cutoff = .8)
drop  <-  names(data_present_climatic_variables[,3:7])[drop_1]
data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]


present_climatic_variables <- dropLayer(present_climatic_variables, names(present_climatic_variables)[drop_1])
future_climatic_variables <- dropLayer(future_climatic_variables, names(future_climatic_variables)[drop_1])


corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number",
         type = "upper")
corrplot(cor(data_future_climatic_variables[3:length(data_present_climatic_variables)]),
         method = "number",
         type = "upper")



# Add field period 
data_present_climatic_variables <- mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)
data <- rbind(data_present_climatic_variables, data_future_climatic_variables)



# Create name object
names <- polygon$NatName


# Mahalanobis distance ----
# Create empty data frame
mh_f <- data.frame(matrix(1,    
                          nrow = nrow(data),
                          ncol = length(names)))

names(mh_f) <- names

for (i in 1:nrow(polygon)){
  pol <- polygon[i,]
  raster_polygon <- raster::mask(raster::crop(present_climatic_variables, pol), pol)
  data_polygon <- raster::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  mh <- mahalanobis(data[,4:length(data)], 
                    colMeans(data_polygon[,3:length(data_polygon)]), 
                    cov(data[,4:length(data)]), 
                    inverted = F)
  
  # Agregar información espacial al reusltado de mh
  mh_f[,i] <- mh
}

# Agregar informacion espacial al resultado de mh
mh <- cbind(data[,c(1:3)], mh_f)


#Creamos raster
#mh_present <- raster::brick()

for(j in 4:length(mh)){
  mh_present <- raster::brick()
  mh_f <- dplyr::filter(mh, Period == "Present")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_present <- raster::stack(mh_present, mh_f)
  crs(mh_present) <- reference_system
  writeRaster(mh_present, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PR_", names[j-3], ".tif"), overwrite=TRUE)
  plot(mh_present)
}




for(j in 4:length(mh)){
  mh_future <- raster::brick()
  mh_f <- dplyr::filter(mh, Period == "Future")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_future <- raster::stack(mh_future, mh_f)
  crs(mh_future) <- reference_system
  writeRaster( mh_future, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/RCP85_2050_", names[j-3],   ".tif"), overwrite=TRUE)
}

toc()

library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(tictoc)

tic()
# Cargar archivos
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PR_Sierra de Gádor.tif")
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/RCP85_2070_Sierra de Gádor.tif")
area_protegida <- sf::st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/KBA/KBA_Gador.shp")


# Convertimos el shapefile a la misma proyección que el raster
area_protegida <- sf::st_transform(area_protegida, crs = terra::crs(mh_raster))

# Extraer los puntos del raster usando terra
puntos_todos <- terra::as.points(mh_raster)

# Convertir los puntos del raster a un objeto 'sf'
puntos_todos <- sf::st_as_sf(puntos_todos)
colnames(puntos_todos) <- c("mh", "geometry")  

# Encontrar los puntos que están dentro del área protegida usando st_intersection
puntos_dentro <- sf::st_intersection(puntos_todos, area_protegida)

#Calcular invers mahalanobis
puntos_todos$inv_mh <-  1 / puntos_todos$mh



# Distancia euclidea ----

# # Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
# dist <- sf::st_distance(puntos_todos, puntos_dentro)
# 
# # Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
# dist <- apply(dist, 1, min)
# 
# # Agregar las distancias mínimas como un atributo a los puntos del raster
# puntos_todos$dist <- dist
# 
# # Convertir las distancias a kilómetros y redondear
# puntos_todos$dist <- round(puntos_todos$dist / 1000, 0)
# puntos_todos$dist[puntos_todos$dist == 0] <- NA




# Autocorrelacion Local ----

# Crear la matriz de pesos espaciales
coords <- st_coordinates(puntos_todos)
puntos_todos <- cbind(puntos_todos, coords)
nb <- dnearneigh(coords, 0, 1000)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

## Moran ----

LM <- localmoran(puntos_todos$inv_mh , lw)

LM_df <- as.data.frame(LM)
puntos_todos$SA <- LM_df$Ii  # Índice local de Moran (Ii)
puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))` # P-valor del índice local de Moran


# Crear una nueva columna con el inverso de mh_guadarrama

puntos_todos$lisa_cluster <- case_when(
  puntos_todos$SA_sig >= 0.05 ~ NA,  # No significativo
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ 1,  # Punto alto rodeado de puntos altos
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA,    # Punto bajo rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ NA,   # Punto alto rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA)

puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)


# Ponderar ----
#resultado <- (scale(puntos_todos$mh)*0.69)+(scale(puntos_todos$dist)*0.2) +
  #(puntos_todos$lisa_cluster *0.1)

#puntos_todos <- cbind(puntos_todos, ponderacion =resultado)


resolution <- 1000
bbox <- st_bbox(puntos_todos)
raster_template <- rast(ext(bbox), nrows = 869, ncols = 1083)
puntos_vect <- vect(puntos_todos)
raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")
crs(raster) <- "EPSG:25830"

# Parches ----
parches <- patches(raster, directions = 8, zeroAsNA=T)
crs(parches) <- crs(raster)
parche_pt <- terra::extract(parches, puntos_todos)
puntos_todos <- puntos_todos %>% 
  mutate(parche = parche_pt$patches)
toc()

# Crear raster----
resolution <- 1000
bbox <- st_bbox(puntos_todos)
raster_template <- rast(ext(bbox), nrows = 869, ncols = 1083)
puntos_vect <- vect(puntos_todos)
raster <- terra::rasterize(puntos_vect, raster_template, field = "SA" )# SA parche)
crs(raster) <- "EPSG:25830"

writeRaster(raster, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PRE_Gador_", "SA",   ".tif"), overwrite=TRUE)


writeRaster(raster, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/RCP85_2070_Gador_", "SA",  ".tif"), overwrite=TRUE)

