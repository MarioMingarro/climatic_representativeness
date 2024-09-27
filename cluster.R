library(raster)   # Para manejar archivos raster
library(sf)       # Para manejar shapefiles
library(gstat)    # Para generar el variograma
library(sp)       # Para trabajar con objetos espaciales
library(RColorBrewer)


# 1. Cargar el archivo raster de elevación con terra
library(terra)
library(sf)
library(tidyverse)

# Cargar el archivo raster usando terra
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/mh_guadarrama.tif")

# 2. Cargar el shapefile del área protegida (usando sf)
area_protegida <- sf::st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/guadarrama_PN.shp")

# 3. Convertimos el shapefile a la misma proyección que el raster (si es necesario)
area_protegida <- sf::st_transform(area_protegida, crs = terra::crs(mh_raster))

# 4. Extraer los puntos del raster usando terra
puntos_todos <- terra::as.points(mh_raster)

# Convertir los puntos del raster a un objeto 'sf'
puntos_todos <- sf::st_as_sf(puntos_todos)

# 5. Encontrar los puntos que están dentro del área protegida usando st_intersection
puntos_dentro <- sf::st_intersection(puntos_todos, area_protegida)

# 6. Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
# Convertir los puntos del polígono a 'sf' si es necesario
distancias_minimas <- sf::st_distance(puntos_todos, puntos_dentro)

# Para obtener la distancia mínima por cada punto, tomamos el valor mínimo de cada fila
distancias_minimas <- apply(distancias_minimas, 1, min)

# 7. Agregar las distancias mínimas como un atributo a los puntos del raster
puntos_todos$distancia_minima <- distancias_minimas

# 8. Convertir las distancias a kilómetros y redondear
puntos_todos$distancia_minima <- round(puntos_todos$distancia_minima / 1000, 0)

# 9. Crear un dataframe con la distancia mínima y el valor de elevación (MH)
bb <- data.frame(
  Distancia = puntos_todos$distancia_minima,
  MH = puntos_todos$mh_guadarrama
)


# Definir un rango de valores K
k_values <- 5:10
sse_values <- numeric(length(k_values))

# Ejecutar K-means para cada K y guardar la SSE
for (k in k_values) {
  tictoc::tic()
  kmeans_result <- kmeans(bb, centers = 10, nstart = 50)
  print(k)
  tictoc::toc()
  sse_values[k] <- kmeans_result$tot.withinss  # SSE
}

# Crear un data frame para graficar
elbow_data <- data.frame(K = k_values, SSE = sse_values)

ggplot(elbow_data, aes(x = K, y = SSE)) +
  geom_line() +
  geom_point() +
  labs(title = "Diagrama del Codo", x = "Número de Clusters (K)", y = "Suma de Errores Cuadráticos (SSE)") +
  theme_minimal()

# Kmeans----
kk <- kmeans(x = bb, centers = 10, nstart = 50)
bb <- bb %>% mutate(cluster = kk$cluster)

bb <- filter(bb, bb$cluster ==1)
min(bb$V1)
max(bb$V1)

bb <- data.frame(
  Distancia = puntos_todos$distancia_minima,
  MH = puntos_todos$mh_guadarrama
)
kk <- kmeans(x = bb, centers = 2, nstart = 50)
bb <- bb %>% mutate(cluster = kk$cluster)

ggplot2::ggplot(bb, aes(x=bb$Distancia, y= bb$MH, 
                        color =  as.factor(bb$cluster)))+
  ggplot2::geom_point()+
  scale_color_brewer(palette = "RdBu") +
  geom_smooth()






# 6. Convertir a un objeto 'SpatialPointsDataFrame' para el semivariograma
puntos_sp <- as(cc, "Spatial")
puntos_sp <- st_as_sf(puntos_sp)

st_write(puntos_sp, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/my_shapefile.shp")




