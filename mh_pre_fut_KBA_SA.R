library(terra)
library(sf)
library(tidyverse)
library(caret)

library(RStoolbox)
library(stringr)
library(corrplot)

library(tictoc)
gc(reset = T)


tic()
# CLIMATE ----
## Load data ----
## Climatic representativeness -----
present_climatic_variables <- terra::rast(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENT/", "\\.tif$", full.names = T))
present_climatic_variables <- present_climatic_variables[[1:17]]

future_climatic_variables <- terra::rast(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/FUTURO/GFDL/", "\\.tif$", full.names = T))
future_climatic_variables <- future_climatic_variables[[1:17]]

names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")
names(future_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                      "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                      "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7")


#present_climatic_variables <- raster::stack(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/PRESENTE/", "\\.tif$", full.names = T))

#future_climatic_variables <- raster::stack(list.files("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/CLIMA/RCP85_2050", "\\.tif$", full.names = T))

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")

polygon <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/KBA/KBA_test.shp")


# Reference system
reference_system <-terra::crs(present_climatic_variables) # "+proj=longlat +datum=WGS84 +no_defs"

if(terra::same.crs(present_climatic_variables, future_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
}

if(terra::same.crs(study_area, present_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  study_area <- st_transform(study_area, crs(reference_system))
  print(paste0("Reference System was modified to ", reference_system))
}

if(terra::same.crs(polygon, present_climatic_variables) == TRUE) {
  print("Same Reference System")
} else {
  polygon <- st_transform(polygon, crs(reference_system))
  print(paste0("Reference System was modified to ", reference_system))
}

# Crop raster to study area
present_climatic_variables <-  terra::mask (crop(present_climatic_variables, study_area), study_area)
future_climatic_variables  <-  terra::mask(crop(future_climatic_variables,  study_area), study_area)
# a ETRS8930N
# reference_system <- "EPSG:25830"
# present_climatic_variables <- terra::project(present_climatic_variables, reference_system)
# future_climatic_variables <- terra::project(future_climatic_variables, reference_system)
# polygon <- sf::st_transform(polygon, crs(reference_system))
# study_area <- sf::st_transform(study_area, crs(reference_system))
# Extract raster data
data_present_climatic_variables <- terra::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- terra::as.data.frame(future_climatic_variables, xy = TRUE)

# Delete NA
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)

# Rename columns
#colnames(data_present_climatic_variables) <- c("x","y","isotermal", "prec_meshum","rango_temp_anual", "temp_max_mescal", "temp_med_anual")

#colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)





# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:length(data_present_climatic_variables)])

# Select variables less correlated 
drop_1  <-  caret::findCorrelation(cor, cutoff = .8)
drop  <-  names(data_present_climatic_variables[,3:length(data_present_climatic_variables)])[drop_1]
data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]

present_climatic_variables <- terra::subset(present_climatic_variables, !names(present_climatic_variables) %in% drop)
future_climatic_variables <- terra::subset(future_climatic_variables, !names(future_climatic_variables) %in% drop)


# corrplot(cor(data_present_climatic_variables[3:length(data_present_climatic_variables)]),
#          method = "number",
#          type = "upper")
# corrplot(cor(data_future_climatic_variables[3:length(data_present_climatic_variables)]),
#          method = "number",
#          type = "upper")



# Add field period 
data_present_climatic_variables <- dplyr::mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- dplyr::mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

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
  raster_polygon <- terra::mask(raster::crop(present_climatic_variables, pol), pol)
  data_polygon <- terra::as.data.frame(raster_polygon, xy = TRUE)
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
  mh_present <- dplyr::filter(mh, Period == "Present")
  mh_present <- terra::rast(mh_present[, c(1:2, j)], crs = reference_system)  # Asignar el CRS aquí
  names(mh_present) <- colnames(mh[j])
  writeRaster(mh_present, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PR_WGS84_", names[j-3], ".tif"), overwrite=TRUE)
  plot(mh_present)
}




for(j in 4:length(mh)){
  mh_future <- terra::rast()
  mh_future <- dplyr::filter(mh, Period == "Present")
  mh_p_rast <- terra::rast(mh_f[, c(1:2, j)], crs = reference_system)  # Asignar el CRS aquí
  names(mh_p_rast) <- colnames(mh[j])
  mh_present <- terra::app(c(mh_present, mh_p_rast), fun = sum, na.rm = TRUE)
  writeRaster(mh_future, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/RCP85_2050_", names[j-3],   ".tif"), overwrite=TRUE)
}

toc()

library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(tictoc)

tic()
# Cargar archivos
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PR_WGS84_Sierra de Gádor.tif")#RCP85_2070_Sierra de Gádor.tif

area_protegida <- sf::st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/KBA/KBA_Gador.shp")

study_area <- read_sf("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/Peninsula_Iberica_89.shp")




# Reference system
reference_system <- "+proj=longlat +datum=WGS84 +no_defs"

if(terra::same.crs(study_area, mh_raster) == TRUE) {
  print("Same Reference System")
} else {
  study_area <- st_transform(study_area, crs(reference_system))
  print(paste0("Reference System was modified to ", reference_system))
}
if(terra::same.crs(area_protegida, mh_raster) == TRUE) {
  print("Same Reference System")
} else {
  area_protegida <- st_transform(area_protegida, crs(reference_system))
  print(paste0("Reference System was modified to ", reference_system))
}


  
  
# Convertimos el shapefile a la misma proyección que el raster
area_protegida <- sf::st_transform(area_protegida, crs(reference_system))

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

# Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
dist <- sf::st_distance(puntos_todos, puntos_dentro)

# Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
dist <- apply(dist, 1, min)

# Agregar las distancias mínimas como un atributo a los puntos del raster
puntos_todos$dist <- dist

# Convertir las distancias a kilómetros y redondear
puntos_todos$dist <- round(puntos_todos$dist / 1000, 0)
puntos_todos$dist[puntos_todos$dist == 0] <- NA


# Crear la matriz de pesos espaciales

## TEST vecinos

# Supongamos que 'puntos_todos' es tu conjunto de datos
kk <- puntos_todos[1:20, ]
coords <- st_coordinates(kk)

# Convertir a dataframe y renombrar columnas
coords_df <- as.data.frame(coords)
colnames(coords_df) <- c("x", "y")

# Calcular vecinos
nb <- spdep::dnearneigh(coords, d1=0, d2=0.015, use_s2=TRUE)

# Crear una lista de líneas de los vecinos
lines_list <- lapply(1:length(nb), function(i) {
  neighbors <- nb[[i]]
  if (length(neighbors) > 0) {
    lapply(neighbors, function(j) {
      if (i != j) {
        return(rbind(coords[i, ], coords[j, ]))
      } else {
        return(NULL)
      }
    })
  } else {
    return(NULL)
  }
})

# Filtrar líneas vacías y crear geometrías válidas
lines_list <- unlist(lines_list, recursive = FALSE)
lines_list <- Filter(Negate(is.null), lines_list)

# Crear un objeto sf solo si hay líneas válidas
if (length(lines_list) > 0) {
  lines_sf <- st_sf(geometry = st_sfc(lapply(lines_list, function(x) {
    if (nrow(x) == 2) {
      st_linestring(x)
    } else {
      NULL
    }
  })))
  
  # Filtrar geometrías nulas
  lines_sf <- lines_sf[!sapply(st_geometry(lines_sf), is.null), ]
  
  # Crear el gráfico solo si hay geometrías válidas
  if (nrow(lines_sf) > 0) {
    ggplot() +
      geom_point(data = coords_df, aes(x = x, y = y), color = "blue", size = 3) +
      geom_sf(data = lines_sf, color = "red") +
      labs(title = "Mapa de Grafos de Vecindades") +
      theme_minimal()
  } else {
    print("No hay líneas válidas para graficar.")
  }
} else {
  print("No hay líneas para graficar.")
}






# Autocorrelacion Local ----

# Crear la matriz de pesos espaciales
coords <- st_coordinates(puntos_todos)
puntos_todos <- cbind(puntos_todos, coords)
kk <- puntos_todos[1:10,]
nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2=TRUE, dwithin=TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)


## Moran ----

LM <- localmoran(puntos_todos$inv_mh , lw)

LM_df <- as.data.frame(LM)
puntos_todos$SA <- LM_df$Ii  # Índice local de Moran (Ii)
puntos_todos$SA_sig <- LM_df$`Pr(z != E(Ii))` # P-valor del índice local de Moran


# Crear una nueva columna con el inverso de mh_guadarrama

puntos_todos$lisa_cluster <- case_when(
  puntos_todos$SA_sig >= 0.01 ~ NA,  # No significativo
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ 1,  # Punto alto rodeado de puntos altos
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA,    # Punto bajo rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ NA,   # Punto alto rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA)

puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)


# Ponderar ----
#resultado <- (scale(puntos_todos$mh)*0.69)+(scale(puntos_todos$dist)*0.2) +
  #(puntos_todos$lisa_cluster *0.1)

#puntos_todos <- cbind(puntos_todos, ponderacion =resultado)



res <- res(mh_raster)  # Resolución del raster existente
bbox <- ext(mh_raster)         # Extensión del raster existente

# Crear un nuevo raster utilizando la extensión y resolución del raster existente
nrows <- round((bbox[4] - bbox[3]) / res[2])  # Número de filas
ncols <- round((bbox[2] - bbox[1]) / res[1])  # Número de columnas
raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)

# Convertir tus puntos a un objeto SpatVector
puntos_vect <- vect(puntos_todos)

# Rasterizar los puntos en el nuevo raster usando el campo deseado
raster <- terra::rasterize(puntos_vect, raster_template, field = "lisa_cluster")

# Establecer el sistema de referencia de coordenadas
crs(raster) <- crs(mh_raster)
plot(raster)




# Parches ----
parches <- patches(raster, directions = 8, zeroAsNA=T)
crs(parches) <- crs(raster)
parche_pt <- terra::extract(parches, puntos_todos)
puntos_todos <- puntos_todos %>% 
  mutate(parche = parche_pt$patches)
toc()

# Crear raster----
# Convertir tus puntos a un objeto SpatVector
puntos_vect <- vect(puntos_todos)

# Rasterizar los puntos en el nuevo raster usando el campo deseado
raster <- terra::rasterize(puntos_vect, raster_template, field = "SA") # parche SA

# Establecer el sistema de referencia de coordenadas
crs(raster) <- crs(mh_raster)



writeRaster(raster, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/PRE_Gador_", "SA_KKKK",   ".tif"), overwrite=TRUE)


writeRaster(raster, paste0("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/MH_RES/RCP85_2070_Gador_", "SA",  ".tif"), overwrite=TRUE)

