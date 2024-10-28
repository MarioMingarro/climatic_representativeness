library(terra)
library(sf)
library(tidyverse)
library(stringr)
library(tictoc)

directorio <- "C:/A_TRABAJO/A_ALBERTO_JIM/"
directorio_result <- "C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/"

regiones <- read_sf(paste0(directorio, "geokarst_reg_final_def.shp"))

# Área total ----
regiones$area_m2 <- st_area(regiones)

# Latitud y longitud máxima y mínima ----
bbox_coords <- function(regiones) {
  bbox <- st_bbox(regiones)
  return(data.frame(
    lon_min = bbox["xmin"],
    lon_max = bbox["xmax"],
    lat_min = bbox["ymin"],
    lat_max = bbox["ymax"]
  ))
}

coords_df <- do.call(rbind, lapply(st_geometry(regiones), bbox_coords))
regiones <- cbind(regiones, coords_df)

regiones <- st_transform(regiones, crs("+proj=longlat +datum=WGS84 +no_defs"))

# aet (Actual Evapotranspiration, monthly total), units = mm
# def (Climate Water Deficit, monthly total), units = mm
# pet (Potential evapotranspiration, monthly total), units = mm
# ppt (Precipitation, monthly total), units = mm
# q (Runoff, monthly total), units = mm
# soil (Soil Moisture, total column - at end of month), units = mm
# srad (Downward surface shortwave radiation), units = W/m2
# swe (Snow water equivalent - at end of month), units = mm
# vap (Vapor pressure, average for month), units  = kPa
# ws (Wind speed, average for month), units = m/s
# vpd (Vapor Pressure Deficit, average for month), units = kpa
# PDSI (Palmer Drought Severity Index, at end of month), units = unitless
# tmax (Max Temperature, average for month), units = C
# tmin (Min Temperature, average for month), units = C
variable_name <- "ppt"
variable <- raster::stack()

tic()
for (i in 1958:2023){
  terraclimate <- raster::brick(
    paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data/TerraClimate_", variable_name,"_",i,".nc"))
  terraclimate <- terra::mask(terra::crop(terraclimate, regiones), regiones)
  if (variable_name %in% c("tmax", "tmin")){
    var <- mean(terraclimate)
  }else{
    var <- sum(terraclimate)
  }
  names(var) <- paste0("Y_", i)
  variable <- raster::stack(variable, var)
}
variable <- mean(variable)
writeRaster(variable, paste0(directorio_result, variable_name, "1958_2023.tif"), overwrite = TRUE)
variable_df <- raster::extract(variable, regiones, fun=mean, na.rm=TRUE, df=TRUE)
regiones[[variable_name]] <- variable_df$layer
toc()


write.csv2(regiones,paste0(directorio_result, "regiones_pcp.csv"))
r <- read.csv2(paste0(directorio_result, "regiones_var.csv"))
r <- r[,-c(1,2,3,4)]

HFI <- terra::rast("C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/HFI/hfp2018.tif")
c <- terra::crs(HFI)
reg <- st_transform(regiones, c)
HFI <- terra::mask(terra::crop(HFI, reg), reg)
HFI <- project(HFI, "+proj=longlat +datum=WGS84 +no_defs")
HFI <- raster::extract(HFI, regiones, fun=mean, na.rm=TRUE)
regiones$HFI <- HFI$hfp2018


NPP <- terra::rast("C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/NPP/CHELSA_npp_1981-2010_AE.tif")
NPP <- terra::mask(terra::crop(NPP, regiones), regiones)
NPP <- raster::extract(NPP, regiones, fun=mean, na.rm=TRUE)
regiones$NPP <- NPP$`CHELSA_npp_1981-2010_AE`


CSI <- terra::rast("C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/ESTABILIDAD/Raster_layers_R_scripts/Layers/past/csi_past.tif")
c <- terra::crs(CSI)
reg <- st_transform(regiones, c)
CSI <- terra::mask(terra::crop(CSI, reg), reg)
CSI <- project(CSI, "+proj=longlat +datum=WGS84 +no_defs")
CSI <- raster::extract(CSI, regiones, fun=mean, na.rm=TRUE)
regiones$CSI <- CSI$csi_past

SDB1 <- terra::rast("C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/ESTABILIDAD/Raster_layers_R_scripts/Layers/past/sd_past_bio_1.tif")
c <- terra::crs(SDB1)
reg <- st_transform(regiones, c)
SDB1 <- terra::mask(terra::crop(SDB1, reg), reg)
SDB1 <- project(SDB1, "+proj=longlat +datum=WGS84 +no_defs")
SDB1 <- raster::extract(SDB1, regiones, fun=mean, na.rm=TRUE)
regiones$SDB1 <- SDB1$sd_past_bio_1

r2 <- regiones[,c(4,13,14)]
kk <- cbind(r,r2)

write.csv2(regiones,paste0(directorio_result, "regiones_var_all_PES.csv"))


TCD <- rast("C:/A_TRABAJO/A_ALBERTO_JIM/VARIABLES/TCD/TCD_2018_100m_es_25830_V1_0.tif")
TCD <- project(TCD, "+proj=longlat +datum=WGS84 +no_defs")
TCD <- terra::mask(terra::crop(TCD, regiones), regiones)
TCD <- raster::extract(TCD, regiones, fun=mean, na.rm=TRUE)
regiones$TCD <- TCD$Class_Name

TCD_df <- raster::extract(TCD, regiones, fun=mean, na.rm=TRUE, df=TRUE)

regiones[[TCD]] <- TCD_df$layer

plot(variable)
plot(regiones$geometry, add=T)

library(reshape2)
library(PupillometryR)
# Excluir la columna 'geometry'
regiones_no_geom <- st_drop_geometry(regiones)

# Aplicar melt() solo a las columnas necesarias, seleccionando la columna 6 y la columna 13 (tmin)
regiones_plot <- melt(regiones_no_geom[, c(6, 14)], id.vars = "UTM_kars_6")

# Convertir UTM_kars_6 a factor para que se trate como una variable categórica
regiones_plot$UTM_kars_6 <- as.factor(regiones_plot$UTM_kars_6)

# Realizar el violin plot agrupado por la columna 6 (UTM_kars_6)
ggplot(regiones_plot,
       aes(
         x = UTM_kars_6,     # Agrupar por la columna 6 (UTM_kars_6)
         y = value,          # Los valores de la columna "tmin" (o la columna que hayas seleccionado)
         fill = UTM_kars_6,
         colour = UTM_kars_6,
         group = UTM_kars_6  # Agrupar explícitamente por UTM_kars_6
       )) +
  geom_flat_violin(
    position = position_nudge(x = .2, y = 0),
    trim = FALSE,
    alpha = 0.5,
    colour = NA
  ) +
  geom_point(position = position_jitter(width = .1),
             size = 2,
             shape = 20,
             alpha =.3) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = .5,
    width = .1,
    colour = "black"
  ) +
  ggtitle("ppt") +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Para mostrar los valores de UTM_kars_6
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank()
  )


