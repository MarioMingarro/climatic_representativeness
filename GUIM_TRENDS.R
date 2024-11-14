
library(terra)
library(tidyverse)
library(sf)
library(spdep)
library(tidyterra)
library(mapSpain)


mh_raster_p <- rast("C:/A_TRABAJO/GUIM/TRENDS/climate_trend.tif")
puntos_todos_p <- terra::as.points(mh_raster_p)
puntos_todos_p <- sf::st_as_sf(puntos_todos_p)
coords <- st_coordinates(puntos_todos_p)
puntos_todos_p <- cbind(puntos_todos_p, coords)
nb <- spdep::dnearneigh(coords, 0, 0.015, use_s2 = TRUE, dwithin = TRUE)
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

LM <- spdep::localmoran(puntos_todos_p$, lw)
LM_df <- as.data.frame(LM)

puntos_todos_p$SA <- LM_df$Ii
puntos_todos_p$SA_sig <- LM_df$`Pr(z != E(Ii))`

puntos_todos_p$SA_cluster <- case_when(
  puntos_todos_p$SA_sig >= 0.05 ~ NA_real_, 
  TRUE ~ puntos_todos_p$layer 
)

puntos_todos_p$SA_cluster <- as.numeric(puntos_todos_p$SA_cluster)
res <- res(mh_raster_p)
bbox <- ext(mh_raster_p)
nrows <- round((bbox[4] - bbox[3]) / res[2])
ncols <- round((bbox[2] - bbox[1]) / res[1])
raster_template <- rast(ext = bbox, nrows = nrows, ncols = ncols)
puntos_vect_p <- vect(puntos_todos_p)
raster <- terra::rasterize(puntos_vect_p, raster_template, field = "SA_cluster")
crs(raster) <- crs(mh_raster_p)
plot(raster)
writeRaster(raster, "C:/A_TRABAJO/GUIM/TRENDS/cluster_SA.tif", overwrite = TRUE)
write.csv2(puntos_todos_p, "C:/A_TRABAJO/GUIM/TRENDS/data_trends.csv")


cuartiles <- quantile(puntos_todos_p$layer, probs = c(0, 0.25, 0.5, 0.75, 1))
puntos_todos_p$SA_cluster <- case_when(
  puntos_todos_p$SA_sig >= 0.05 ~ 0,  # No significativo
  puntos_todos_p$layer <= cuartiles[1] ~ 1,  # Mínimo
  puntos_todos_p$layer <= cuartiles[2] ~ 2,  # Primer cuartil
  puntos_todos_p$layer <= cuartiles[3] ~ 3,  # Segundo cuartil (mediana)
  puntos_todos_p$layer <= cuartiles[4] ~ 4,  # Tercer cuartil
  puntos_todos_p$layer <= cuartiles[5] ~ 5,  # Cuarto cuartil
)
puntos_vect_p <- vect(puntos_todos_p)
raster <- terra::rasterize(puntos_vect_p, raster_template, field = "SA_cluster")
crs(raster) <- crs(mh_raster_p)
plot(raster)
writeRaster(raster, "C:/A_TRABAJO/GUIM/TRENDS/cluster_SA_Quart.tif", overwrite = TRUE)

prov <- esp_get_prov()
prov <- prov %>% filter(prov$nuts.prov.code != "NA")
ggplot(prov) +
  geom_spatraster(data = raster) +
  geom_spatvector(col = "grey60", fill = NA) +
  scale_fill_whitebox_c(
    palette = "muted", direction = 1,
    labels = scales::label_number(suffix = "º/year")
  ) +
  theme_minimal() +
  coord_sf(crs = 25830) +
  labs(
    fill = "Trend"
  )
ggplot(prov) +
  geom_spatraster(data = raster) +
  geom_spatvector(col = "grey60", fill = NA) +
  scale_fill_whitebox_c(
    palette = "muted", direction = 1,
    labels = scales::label_number(suffix = "º Q")
  ) +
  theme_minimal() +
  coord_sf(crs = 25830)
  

