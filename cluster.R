library(raster)   # Para manejar archivos raster
library(sf)       # Para manejar shapefiles
library(gstat)    # Para generar el variograma
library(sp)       # Para trabajar con objetos espaciales
library(dplyr)
library(RColorBrewer)
library(ggplot2)


# 1. Cargar el archivo raster de elevación
mh_raster <- raster("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/mh_guadarrama.tif")

# 2. Cargar el shapefile del área protegida (usando sf)
area_protegida <- st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/guadarrama_PN.shp")

# Convertimos el shapefile a la misma proyección que el raster (si es necesario)
area_protegida <- st_transform(area_protegida, crs = crs(mh_raster))

# 3. Extraer los puntos del raster
puntos_todos <- rasterToPoints(mh_raster, spatial = TRUE)

# Convertir los puntos del raster a un objeto 'sf'
puntos_todos <- st_as_sf(puntos_todos)

# 4. Encontrar los puntos que están dentro del área protegida
# Utilizamos la función st_intersection para seleccionar los puntos dentro del polígono
puntos_dentro <- st_intersection(puntos_todos, area_protegida)

# 5. Convertir de vuelta a 'Spatial' para calcular distancias (si es necesario)
puntos_dentro <- as(puntos_dentro, "Spatial")

# 6. Calcular la distancia mínima desde cada punto del raster a los puntos dentro del polígono
puntos_todos <- as(puntos_todos, "Spatial")

# Para cada punto del raster, calculamos la distancia mínima a un punto dentro del polígono
distancias_minimas <- apply(spDists(puntos_todos, puntos_dentro), 1, min)

# 5. Agregar las distancias mínimas como un atributo a los puntos del raster
puntos_todos$distancia_minima <- distancias_minimas
puntos_todos$distancia_minima <- round(puntos_todos$distancia_minima/1000,0)

# numero de clusters----
bb <- as.data.frame(cbind(puntos_todos$distancia_minima, puntos_todos$mh_guadarrama))
colnames(bb) <- c("Distancia", "MH")

# Definir un rango de valores K
k_values <- 5:10
sse_values <- numeric(length(k_values))

# Ejecutar K-means para cada K y guardar la SSE
for (k in k_values) {
  tictoc::tic()
  kmeans_result <- kmeans(bb, centers = k, nstart = 50)
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
kk <- kmeans(x = bb, centers = 15, nstart = 50)
bb <- bb %>% mutate(cluster = kk$cluster)

bb <- filter(bb, bb$cluster ==1)
min(bb$V1)
max(bb$V1)



ggplot2::ggplot(cc, aes(x=cc$distancia_minima, y= cc$mh_guadarrama, 
                        color =  as.factor(cc$cluster)))+
  ggplot2::geom_point()+
  scale_color_brewer(palette = "RdBu") +
  geom_smooth()






# 6. Convertir a un objeto 'SpatialPointsDataFrame' para el semivariograma
puntos_sp <- as(cc, "Spatial")
puntos_sp <- st_as_sf(puntos_sp)

st_write(puntos_sp, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/my_shapefile.shp")



# Cargar librerías necesarias
library(ggplot2)


