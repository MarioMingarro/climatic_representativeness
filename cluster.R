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


# Número de clusters ----
# Calcular la TSS (Suma total de cuadrados)
tss <- sum((bb - colMeans(bb))^2)

# Valores de K a probar
k_values <- 1:10
sse_values <- numeric(length(k_values))
varianza_explicada <- numeric(length(k_values))

# Ejecutar K-means para cada K y guardar la SSE
for (k in k_values) {
  tictoc::tic()
  kmeans_result <- kmeans(bb, centers = k, nstart = 50)
  print(k)
  tictoc::toc()
  
  sse_values[k] <- kmeans_result$tot.withinss  # SSE
  varianza_explicada[k] <- 1 - (sse_values[k] / tss)  # Porcentaje de varianza explicada
}

# Crear un data frame para graficar
elbow_data <- data.frame(K = k_values, VarianzaExplicada = varianza_explicada,  SSE = sse_values)

## Varianza ----
# Graficar el diagrama del codo
ggplot(elbow_data, aes(x = K, y = VarianzaExplicada)) +
  geom_line() +
  geom_point() +
  labs(title = "Diagrama del Codo - Varianza Explicada", 
       x = "Número de Clusters (K)", 
       y = "Porcentaje de Varianza Explicada") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent)

## SSE ----
ggplot(elbow_data, aes(x = K, y = SEE)) +
  geom_line() +
  geom_point() +
  labs(title = "Diagrama del Codo - Varianza Explicada", 
       x = "Número de Clusters (K)", 
       y = "Porcentaje de Varianza Explicada") +
  theme_minimal() +
  scale_y_continuous()

# Kmeans----
kk <- kmeans(x = bb, centers = 6, nstart = 50)
bb <- bb %>% mutate(cluster = kk$cluster)

bb <- data.frame(
  Distancia = puntos_todos$distancia_minima,
  MH = puntos_todos$mh_guadarrama
)


# Ver clusters----
ggplot2::ggplot(bb, aes(x=bb$Distancia, y= bb$MH, 
                        color =  as.factor(bb$cluster)))+
  ggplot2::geom_point()+
  scale_color_brewer(palette = "RdBu") +
  geom_smooth()

puntos_sp <- puntos_todos %>% mutate(bb)

class(puntos_sp)

st_write(puntos_sp, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/my_shapefile.shp")

dentro <- filter(puntos_sp, puntos_sp$Distancia == 0)
max(dentro$MH)
aa <- filter(puntos_sp, puntos_sp$cluster == 1 & puntos_sp$MH <= max(dentro$MH))



# Visualizar el objeto sf llamado 'aa' sobre un mapa base de relieve (Stamen Terrain)
ggplot() +
  
  # Mostrar los datos de 'aa' (polígonos o puntos)
  geom_sf(data = aa, aes(geometry = geometry), color = "blue", fill = NA, size = 1) +
  # Título y ajustes de visualización
  labs(title = "aaa", 
       subtitle = "Visualización de datos sf sobre mapa") +
  theme_minimal() +
  coord_sf()

library(tmap)

# Cambiar el modo de tmap a interactivo
tmap_mode("view")

# Visualizar las CCAA con un mapa de relieve interactivo
tm_shape(aa) +
  tm_borders(col = "blue") +  # Mostrar los bordes de las CCAA
  tm_fill(col = "lightblue", alpha = 0.5) +  # Rellenar las CCAA con color
  tm_basemap() +  # Agregar mapa base de relieve
  tm_layout(title = "Ckkkkk")
