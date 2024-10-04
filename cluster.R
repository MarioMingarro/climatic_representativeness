library(terra)
library(sf)
library(tidyverse)
library(spdep)
library(tictoc)

tic()
# Cargar archivos
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/mh_guadarrama.tif")
mh_raster <- terra::rast("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/mh_guadarrama_2050.tif")

area_protegida <- sf::st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/guadarrama_PN.shp")

# Convertimos el shapefile a la misma proyección que el raster
area_protegida <- sf::st_transform(area_protegida, crs = terra::crs(mh_raster))

# Extraer los puntos del raster usando terra
puntos_todos <- terra::as.points(mh_raster)

# Convertir los puntos del raster a un objeto 'sf'
puntos_todos <- sf::st_as_sf(puntos_todos)
colnames(puntos_todos) <- c("mh", "geometry")  

# Encontrar los puntos que están dentro del área protegida usando st_intersection
puntos_dentro <- sf::st_intersection(puntos_todos, area_protegida)


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
#Calcular invers mahalanobis
puntos_todos$inv_mh <-  1 / puntos_todos$mh

# Autocorrelacion Local ----
# Crear la matriz de pesos espaciales
coords <- st_coordinates(puntos_todos)
puntos_todos <- cbind(puntos_todos, coords)
nb <- dnearneigh(coords, 0, 1000)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

## Moran ----

G30 <- localmoran(puntos_todos$inv_mh , lw)

G30_df <- as.data.frame(G30)
puntos_todos$SA <- G30_df$Ii  # Índice local de Moran (Ii)
puntos_todos$SA_sig <- G30_df$`Pr(z != E(Ii))` # P-valor del índice local de Moran


# Crear una nueva columna con el inverso de mh_guadarrama

puntos_todos$lisa_cluster <- case_when(
  puntos_todos$SA_sig >= 0.01 ~ NA,  # No significativo
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ 1,  # Punto alto rodeado de puntos altos
  puntos_todos$mh < max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA,    # Punto bajo rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA > 0 ~ NA,   # Punto alto rodeado de puntos bajos
  puntos_todos$mh > max(puntos_dentro$mh) & puntos_todos$SA < 0 ~ NA)

puntos_todos$lisa_cluster <- as.numeric(puntos_todos$lisa_cluster)


# Ponderar ----
resultado <- (scale(puntos_todos$mh)*0.69)+(scale(puntos_todos$dist)*0.2) +
  (puntos_todos$lisa_cluster *0.1)
puntos_todos <- cbind(puntos_todos, ponderacion =resultado)
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
raster <- terra::rasterize(puntos_vect, raster_template, field = "parche")
crs(raster) <- "EPSG:25830"

st_write(puntos_todos, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/puntos_todos_parches_2050.shp", append = FALSE)

###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

plot(raster)

#puntos_todos <- puntos_todos %>%
#select(c(1, 2, 3, 4, 5, 6, 7, 8, 9,  11,10))

puntos_todos <-  st_read("C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/puntos_todos.shp")
colnames(puntos_todos) <- c("mh_guadarrama","distancia_minima","inv_mh","X","Y","SA","SA_sig","geometry")


writeRaster(raster, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/raster_resultado.tif")
writeRaster(raster, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/raster_resultado.tif", overwrite=TRUE)

# Plotear LISA -----
# Definir los colores para cada nivel del factor lisa_cluster
colors <- c("1" = "red",     # Color para el nivel 1
            "2" = "green",   # Color para el nivel 2
            "3" = "blue",    # Color para el nivel 3
            "4" = "yellow")  # Color para el nivel 4
# Graficar
ggplot() +
  # Mostrar los datos de 'puntos_todos'
  geom_sf(data = puntos_todos, aes(geometry = geometry, color = factor(lisa_cluster)), size = 1) +
  # Título y ajustes de visualización
  labs(title = "LISA",
       subtitle = "p > 0.05") +
  theme_minimal() +
  coord_sf() +
  # Aplicar los colores definidos
  scale_color_manual(values = colors,
                     name = "LISA Cluster",
                     labels = c("1", "2", "3", "4")) +
  theme(legend.position = "right")  # Ajustar la posición de la leyenda

df_resultados <- puntos_todos %>%
  group_by(parches) %>%  # Agrupamos por el ID del parche
  summarize(
    area_total = n(),         # Sumamos el área total de cada grupo
    ponderacion_promedio = mean(ponderacion, na.rm = TRUE)  # Calculamos el valor promedio del raster en cada grupo
  )%>%
  mutate(FIN = 1/(area_total / res_promedio))
puntos_vect <- vect(df_resultados)
raster <- terra::rasterize(puntos_vect, raster_template, field = "FIN")
fin <- terra::extract(raster, puntos_todos)
View(fin)
df_resultados$ID <- as.numeric(rownames(df_resultados))
df_resultados <- left_join(as.data.frame(df_resultados), as.data.frame(fin))
puntos_todos <- left_join(as.data.frame(puntos_todos), as.data.frame(df_resultados), by = c("patches" = "patches"))

write.csv2(puntos_todos, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/resultado.csv")

# Crear raster----
resolution <- 1000
bbox <- st_bbox(puntos_todos)
raster_template <- rast(ext(bbox), nrows = 869, ncols = 1083)
puntos_vect <- vect(puntos_todos)
raster <- terra::rasterize(puntos_vect, raster_template, field = "RES")
crs(raster) <- "EPSG:25830"


# Parches ----
parches <- patches(raster, directions = 8, zeroAsNA=T)
crs(parches) <- crs(raster)
#ii <- terra::as.points(parches)
#ii <- sf::st_as_sf(ii)

pp <- terra::extract(parches, puntos_todos)
pp <- cbind(puntos_todos, ii)

# Calculos----

df_resultados <- pp %>%
  group_by(patches) %>%  # Agrupamos por el ID del parche
  summarize(
    area_total = n(),         # Sumamos el área total de cada grupo
    res_promedio = mean(RES, na.rm = TRUE)  # Calculamos el valor promedio del raster en cada grupo
  )%>%
  mutate(FIN = 1/(area_total / res_promedio))

resolution <- 1000
bbox <- st_bbox(puntos_todos)
raster_template <- rast(ext(bbox), nrows = 869, ncols = 1083)
puntos_vect <- vect(df_resultados)
raster <- terra::rasterize(puntos_vect, raster_template, field = "FIN")
crs(raster) <- "EPSG:25830"
fin <- terra::extract(raster, puntos_todos)
plot(raster)

pp <- left_join(as.data.frame(pp), as.data.frame(df_resultados))

###############################

aa <- puntos_todos

bb <- as.data.frame(parches, xy = TRUE, na.rm = TRUE)
bb$x <- round(bb$x , 0)
bb$y <- round(bb$y , 0)

cc <- left_join(df_raster, df_resultados, by = c("patches" = "patches"))


bb <- st_join(aa, kk, by = c("X" = "x", "Y" = "y"))

##################





head(puntos_todos)
# Calcular el área de cada celda
areas_celdas <- cellSize(kkvk, unit = "m")
plot(areas_celdas)
# Convertir el raster de parches y el raster original en data frame para manipular con dplyr
df_raster <- as.data.frame(c(kkvk, raster), xy = TRUE, na.rm = TRUE)
# Agregar el área a df_raster utilizando el valor de las celdas de kkvk
df_raster$area <- areas_celdas[kkvk[]]

# Agregar el valor del raster original al data frame
df_raster$raster_value <- rast_t[kkvk[]]



colnames(puntos_todos) <- c("x", "y","patches","mh","area" , "kk")

# Ahora, agrupamos por el ID de cada parche (kkvk) y calculamos el promedio del valor del raster y el área total
df_resultados <- puntos_todos %>%
  group_by(grupo) %>%  # Agrupamos por el ID del parche
  summarize(
    area_total = n(),         # Sumamos el área total de cada grupo
    RES_medio = mean(RES, na.rm = TRUE)  # Calculamos el valor promedio del raster en cada grupo
  )

df_resultados$patches
df_resultados %>%
  mutate(FIN = valor_promedio / area_total)
pp <- as.data.frame(pp)

# Ahora puedes hacer el left_join
df_resultados <- left_join(pp, df_resultados, by = "patches")

# Ahora dividimos el valor promedio del raster entre el área total de cada grupo
df_resultados <- left_join(pp, df_resultados, by = c("patches" = "patches"))



# Ver los resultados
print(df_resultados)

# Exportar ----
writeRaster(raster, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/pedo.tif")
st_write(puntos_todos, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/resultados4.shp")
ggplot() +
  geom_sf(data = puntos_todos, aes(geometry = geometry, color = puntos_todos$RES ), size = 1) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "aaa", 
       subtitle = "Visualización de datos sf sobre mapa",
       color = "Valor de Agg") +  
  theme_minimal() +
  coord_sf()  

ggplot(pp) +
  geom_sf(aes(geometry = geometry, color = pp$patches)) +
  theme_minimal() +
  coord_sf()



##################################################
#################################################

## Getis-Ord G ----
tic()
G30 <- localG(puntos_todos$mh_guadarrama, lw)
toc()
SAlocalI <- as.data.frame(attr(G30, "internals"))
puntos_todos$SA <- SAlocalI$Gi
puntos_todos$SA_sig <- SAlocalI$`Pr(z != E(Gi))`


# Hotspots ----
hotspots <- hotspot(G30, Prname = "Pr(z != E(Gi))", cutoff = 0.05)  # Hotspots con corte de 0.05
puntos_todos$hotspot <- as.character(hotspots) 
levels(hotspots)
unique(droplevels(hotspots))
puntos_todos$hotspot_binary <- ifelse(puntos_todos$hotspot == "High", 1, 0)

# Crear variables estandarizadas de mh_guadarrama (z-scores)
puntos_todos$mh_guadarrama_z <- as.numeric(scale(puntos_todos$mh_guadarrama))

# Crear la variable de lag espacial de mh_guadarrama_z usando pesos espaciales
puntos_todos$mh_guadarrama_zlag <- lag.listw(lw, puntos_todos$mh_guadarrama_z)

# Agrupar en categorías "High-High", "Low-Low", etc. usando el valor estandarizado y la autocorrelación local (SA)
puntos_todos$lisa_cluster <- case_when(
  puntos_todos$SA_sig >= 0.05 ~ 0,  # No significativo si p > 0.05
  puntos_todos$mh_guadarrama_z > 0 & puntos_todos$SA > 0 ~ 1,  # Alto valor y autocorrelación alta
  puntos_todos$mh_guadarrama_z > 0 & puntos_todos$SA < 0 ~ 0,   # Alto valor y autocorrelación baja
  puntos_todos$mh_guadarrama_z < 0 & puntos_todos$SA > 0 ~ 0,    # Bajo valor y autocorrelación alta
  puntos_todos$mh_guadarrama_z < 0 & puntos_todos$SA < 0 ~ 0    # Bajo valor y autocorrelación baja
)
puntos_todos$lisa_cluster <- ifelse(puntos_todos$lisa_cluster  == "High-high", 1, 0)






#################################################################################

bb <- data.frame(
  Distancia = puntos_todos$distancia_minima,
  MH = puntos_todos$mh_guadarrama,
  SA = puntos_todos$SAlocalI)

ggplot(bb, aes(x = Distancia, y= MH, color = agg)) + 
  geom_point() + 
  scale_color_gradient(low = "blue", high = "red")



# 7. Convertir el raster a un data frame para ggplot
raster_df <- as.data.frame(raster_localI, xy = TRUE, na.rm = TRUE)
# Cambiamos el nombre de la columna con los valores de 'localI'
colnames(raster_df)[3] <- "RES"  # La tercera columna es donde están los valores
# 8. Graficar los resultados con ggplot2
ggplot() +
geom_raster(data = raster_df, aes(x = x, y = y, fill = RES)) +
scale_fill_viridis_c() +  # Escala de color continua
theme_minimal() +
labs(title = "Índice de Moran Bivariado Local",
x = "Coordenada X", y = "Coordenada Y",
fill = "Índice de Moran Local") +
coord_equal()



head(puntos_todos$hotspot)

# 9. Crear un dataframe con la distancia mínima y el valor de elevación (MH)
bb <- data.frame(
  Distancia = puntos_todos$distancia_minima,
  MH = puntos_todos$mh_guadarrama,
  SA = puntos_todos$hotspot_binary)



# Coordenadas geográficas
#######################################################################################
# Instalar los paquetes necesarios si no los tienes
install.packages(c("sf", "spdep", "terra", "ggplot2", "stars", "tidyverse"))

# Cargar las librerías
library(sf)
library(spdep)
library(terra)
library(ggplot2)
library(stars)
library(dplyr)
library(tictoc)
# 1. Cargar los datos
puntos <- puntos_todos

# 2. Crear matriz de pesos espaciales (Vecinos y pesos)
coords <- st_coordinates(puntos)
nb30 <- dnearneigh(coords, 0, 1000)  # Vecinos dentro de un radio de 1000 metros
lw <- nb2listw(nb30, style = "B", zero.policy = TRUE)  # Creamos la lista de pesos espaciales

# 3. Calcular el índice local de Getis-Ord G usando la variable mh_guadarrama
# Esto devuelve tanto los valores G como los valores p (significancia) usando permutaciones
tic()
G30 <- localG(puntos$mh_guadarrama,  nb2listw(nb30, style = "B", zero.policy = TRUE))  # nsim define el número de permutaciones
toc()
# 4. Extraer los valores G y las probabilidades p de significancia
G30_res <- as.data.frame(attr(G30, "internals"))

puntos$G_local <- G30_res$Gi  # Índice G local basado en la variable mh_guadarrama
puntos$p_value <- G30_res$`Pr(z != E(Gi))`  # Valores p

# 5. Identificar hotspots (clústeres) basados en los valores de p
# Usamos la función hotspot y especificamos el nombre de la columna que tiene los valores p
hotspots <- hotspot(G30, Prname = "Pr(z != E(Gi))", cutoff = 0.05)  # Hotspots con corte de 0.05

# 6. Agregar los resultados de hotspots a los datos espaciales
puntos$hotspot_G <- hotspots  # Agregar la clasificación de hotspots como nueva columna

# 7. Crear un raster de 1000 metros de resolución en el sistema de coordenadas actual (ETRS89 / UTM zone 30N)
resolution <- 1000  # 1000 metros

# Definir el bounding box del área de estudio
bbox <- st_bbox(puntos)

# Crear un raster vacío con la resolución adecuada usando 'terra'
raster_template <- rast(ext(bbox), res = c(resolution, resolution))

# 8. Rasterizar los puntos usando las categorías de clúster (hotspots)
# Convertimos el objeto sf a SpatVector para usar con terra
puntos_vect <- vect(puntos)

# Rasterizar la categoría sobre el raster template
raster_hotspot <- rasterize(puntos_vect, raster_template, field = "hotspot_G")

# 9. Convertir el raster a un data frame para ggplot
raster_df <- as.data.frame(raster_hotspot, xy = TRUE, na.rm = TRUE)

# Asegurarnos de que la columna 'hotspot_G' sea de tipo factor
raster_df$hotspot_G <- as.factor(raster_df$hotspot_G)

# 10. Graficar los resultados con ggplot, utilizando la clasificación de hotspots
ggplot() +
  geom_tile(data = raster_df, aes(x = x, y = y, fill = hotspot_G)) +
  scale_fill_manual(values = c("1" = "red", "0" = "blue"), labels = c("Hotspot", "No Hotspot")) +
  theme_minimal() +
  labs(title = "Hotspots de Autocorrelación Local (Getis-Ord G)",
       x = "Coordenada X (UTM)", y = "Coordenada Y (UTM)",
       fill = "Clúster G") +
  coord_equal() +
  theme(legend.position = "right")
writeRaster(raster_hotspot, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/hotspot.tif")






# 2. Crear la matriz de pesos espaciales
coords <- st_coordinates(puntos)
nb <- spdep::dnearneigh(coords, 0, 10000)  # Vecinos dentro de un radio de 1000 metros
listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

# 3. Variables para el análisis
var1 <- puntos$mh_guadarrama  # Primera variable cuantitativa
var2 <- puntos$distancia_minima  # Segunda variable cuantitativa

# 4. Cálculo del Índice de Moran Bivariado Local (manual)
# Normalizamos las variables (opcional pero recomendado)
var1_norm <- (var1 - mean(var1)) / sd(var1)
var2_norm <- (var2 - mean(var2)) / sd(var2)

# Producto cruzado entre las dos variables
cross_product <- var1_norm * lag.listw(listw, var2_norm, zero.policy = TRUE)

# Cálculo del Índice de Moran Bivariado Local para cada punto
local_bivariate_I <- cross_product

# Añadir los resultados al objeto sf
puntos$localI <- local_bivariate_I

# 5. Clasificar las autocorrelaciones locales en HH, LL, HL, LH
# Calcular el valor promedio ponderado en los vecinos
var1_lag <- lag.listw(listw, var1_norm, zero.policy = TRUE)
var2_lag <- lag.listw(listw, var2_norm, zero.policy = TRUE)

# Clasificar los puntos según su relación con sus vecinos
puntos$category <- case_when(
  var1_norm > 0 & var1_lag > 0 ~ "HH",  # High-High
  var1_norm > 0 & var1_lag < 0 ~ "HL",  # High-Low
  var1_norm < 0 & var1_lag < 0 ~ "LL",  # Low-Low
  var1_norm < 0 & var1_lag > 0 ~ "LH"   # Low-High
)


ggplot(data = puntos) +
  geom_point(aes(x = st_coordinates(geometry)[,1], 
                 y = st_coordinates(geometry)[,2], 
                 color = category), size = 0.5) +
  scale_color_manual(values = c("HH" = "red", "HL" = "blue", "LL" = "green", "LH" = "yellow")) +
  theme_minimal() +
  labs(title = "Clasificación de Autocorrelación Local (HH, HL, LL, LH)",
       x = "Coordenada X (UTM)", y = "Coordenada Y (UTM)",
       color = "Categoría") +
  coord_equal() +
  theme(legend.position = "right")

ggplot(data = puntos) +
  geom_point(aes(x = st_coordinates(geometry)[,1], 
                 y = st_coordinates(geometry)[,2], 
                 color = puntos$G_local), size = 0.5) +
  theme_minimal() +
  labs(title = "Clasificación de Autocorrelación Local (HH, HL, LL, LH)",
       x = "Coordenada X (UTM)", y = "Coordenada Y (UTM)",
       color = "Categoría") +
  coord_equal() +
  theme(legend.position = "right")

# 6. Crear un raster de 1000 metros de resolución
# Definir la resolución deseada (1000 metros)
resolution <- 1000

# Definir el bounding box del área de estudio
bbox <- st_bbox(puntos)

# Crear un raster vacío con la resolución y bounding box adecuado usando 'terra'
raster_template <- rast(ext(bbox), resolution, resolution)

# 7. Rasterizar los puntos usando las categorías HH, HL, LL, LH
# Convertimos el objeto sf a SpatVector para usar con terra
puntos_vect <- vect(puntos)

# Rasterizar la categoría sobre el raster template
raster_category <- rasterize(puntos_vect, raster_template, field = "category")

# 8. Convertir el raster a un data frame para ggplot
raster_df <- as.data.frame(raster_category, xy = TRUE, na.rm = TRUE)

# Cambiamos el nombre de la columna con los valores de categoría
# Asegurarnos de que la columna 'category' sea de tipo factor
raster_df$category <- as.factor(raster_df$category)

# 9. Graficar los resultados con geom_tile() para ver las categorías
ggplot() +
  geom_tile(data = raster_df, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = c("HH" = "red", "HL" = "blue", "LL" = "green", "LH" = "yellow")) +
  theme_minimal() +
  labs(title = "Clasificación de Autocorrelación Local (HH, HL, LL, LH)",
       x = "Coordenada X", y = "Coordenada Y",
       fill = "Categoría") +
  coord_equal()

# 5. Crear un raster de 1000 metros de resolución
# Definir la resolución deseada (1000 metros)
resolution <- 1000

# Definir el bounding box del área de estudio
bbox <- st_bbox(puntos)

# Crear un raster vacío con la resolución y bounding box adecuado usando 'terra'
raster_template <- rast(ext(bbox), resolution, resolution)

# 6. Rasterizar los puntos usando los valores del índice de Moran Local
# Convertimos el objeto sf a SpatVector para usar con terra
puntos_vect <- vect(puntos)

# Rasterizar el índice de Moran Local sobre el raster template
raster_localI <- rasterize(puntos_vect, raster_template, field = "localI")

# 7. Convertir el raster a un data frame para ggplot
raster_df <- as.data.frame(raster_localI, xy = TRUE, na.rm = TRUE)

# Cambiamos el nombre de la columna con los valores de 'localI'
colnames(raster_df)[3] <- "localI"  # La tercera columna es donde están los valores

# 8. Graficar los resultados con ggplot2
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = localI)) +
  scale_fill_viridis_c() +  # Escala de color continua
  theme_minimal() +
  labs(title = "Índice de Moran Bivariado Local",
       x = "Coordenada X", y = "Coordenada Y",
       fill = "Índice de Moran Local") +
  coord_equal()
writeRaster(raster_localI, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/kk.tif")
plot(raster_localI
     )
##############################################################################




# Número de clusters ----
bb <- data.frame(
  LiM = puntos$localI,
  MH = puntos$mh_guadarrama
)
# Calcular la TSS (Suma total de cuadrados)
tss <- sum((puntos - colMeans(bb))^2)

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
bb <- data.frame(
  SA = puntos_todos,
  MH = puntos$mh_guadarrama
)

kk <- kmeans(x = na.omit(puntos_todos$RES), centers = 6, nstart = 50)
puntos
bb <- bb %>% mutate(cluster = kk$cluster) %>% mutate(puntos)

bb <- bb[,-c(4,7)]

# Ver clusters----
ggplot2::ggplot(bb, aes(x=bb$distancia_minima , y= bb$mh_guadarrama , 
                        color =  as.factor(bb$cluster)))+
  ggplot2::geom_point()+
  scale_color_brewer(palette = "RdBu")

puntos_sp <- puntos_todos %>% mutate(bb)

class(puntos_sp)

st_write(puntos_todos, "C:/A_TRABAJO/A_GABRIEL/REPRESENTATIVIDAD/ponderados.shp")

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
