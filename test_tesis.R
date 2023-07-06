library(raster)
library(sf)
library(tidyverse)
library(RStoolbox)
library(stringr)
library(corrplot)
library(caret)

gc(reset = T)

# Load data ----
# Climatic representativeness -----
present_climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Presente", ".tif", full.names = T))

# Climatic representativeness -----
future_climatic_variables <- raster::stack(list.files("T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Futuro/RCP85/2050", ".tif", full.names = T))

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "//.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/national_parks_2.shp")


# Reference system
reference_system <- projection(present_climatic_variables) # "+proj=longlat +datum=WGS84 +no_defs"

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
}


# Crop raster to study area
future_climatic_variables <- resample(future_climatic_variables, present_climatic_variables, method='bilinear')
present_climatic_variables <-  raster::mask(crop(present_climatic_variables, future_climatic_variables), future_climatic_variables)
future_climatic_variables  <-  raster::mask(crop(future_climatic_variables,  present_climatic_variables), present_climatic_variables)

# Extract raster data
data_present_climatic_variables <- raster::as.data.frame(present_climatic_variables, xy = TRUE)
data_future_climatic_variables <- raster::as.data.frame(future_climatic_variables, xy = TRUE)

# Delete NA
data_present_climatic_variables<-na.omit(data_present_climatic_variables)
data_future_climatic_variables <-na.omit(data_future_climatic_variables)

# Rename columns
colnames(data_present_climatic_variables) <- c("x", "y", "isotermal", "prec_meshum","rango_temp_anual","temp_max_mescal","temp_med_anual")

colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)

# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:7])

#Select variables less correlated 
drop  <-  findCorrelation(cor, cutoff = .8)
drop  <-  names(data_present_climatic_variables[,3:21])[drop]
data_present_climatic_variables <- data_present_climatic_variables[!names(data_present_climatic_variables) %in% drop]
data_future_climatic_variables <- data_future_climatic_variables[!names(data_future_climatic_variables) %in% drop]


present_climatic_variables <- dropLayer(present_climatic_variables, names(present_climatic_variables)[drop_1])
future_climatic_variables <- dropLayer(future_climatic_variables, names(future_climatic_variables)[drop_1])


corrplot(cor(data_present_climatic_variables[3:7]),
         method = "number",
         type = "upper")
corrplot(cor(data_future_climatic_variables[3:7]),
         method = "number",
         type = "upper")

# Add field period 
data_present_climatic_variables <- mutate(data_present_climatic_variables, Period = c("Present"),  .after = "y")
data_future_climatic_variables  <- mutate(data_future_climatic_variables, Period = c("Future"),  .after = "y")

# Join two dataset
data <- rbind(data_present_climatic_variables, data_future_climatic_variables)




# Create name object
names <- polygon$ORIG_NAME


# Mahalanobis distance ----
# Create empty data frame
mh_f <- data.frame(matrix(1,    
                          nrow = nrow(data),
                          ncol = length(names)))

names(mh_f) <- names

for (i in 1:length(names)){
  pol <- polygon[i,]
  raster_polygon <- raster::mask(raster::crop(present_climatic_variables, pol), pol)
  data_polygon <- raster::as.data.frame(raster_polygon, xy = TRUE)
  data_polygon <- na.omit(data_polygon)
  
  colnames(data_polygon) <- c("x","y","isotermal","prec_meshum","rango_temp_anual","temp_max_mescal","temp_med_anual")
  
  mh <- mahalanobis(data[,4:length(data)], 
                    colMeans(data_polygon[,3:7]), 
                    cov(data[,4:length(data)]), 
                    inverted = F)
  
  # Agregar información espacial al reusltado de mh
  mh_f[,i] <- mh
}

# Agregar informacion espacial al reusltado de mh
mh <- cbind(data[,c(1:3)], mh_f)


#Creamos raster
mh_present<- raster::brick()

for(j in 4:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Present")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_present <- raster::stack(mh_present, mh_f)
}
plot(mh_present)

mh_futuro <- raster::brick()

for(j in 4:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Future")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_futuro <- raster::stack(mh_futuro, mh_f)
}

for ( i in 1:nlayers(mh_present)){
  writeRaster(mh_present[[i]], paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/PN/kk/mh_present_", names[i], "_.tif"), overwrite=TRUE)
  writeRaster( mh_futuro[[i]], paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/PN/kk/mh_future_", names[i],   "_.tif"), overwrite=TRUE)
}


# Threshold selection
mh_present_umbral <- raster::brick()
mh_futuro_umbral <- raster::brick()

for (i in 1:nlayers(mh_present)){
  mh_poligono <- mask(crop(mh_present[[i]], polygon[i,]), polygon[i,])
  mh_poligono <- raster::as.data.frame(mh_poligono, xy = T)
  mh_poligono <- quantile(na.omit(mh_poligono[,3]), probs = c(.95))
  mh_present_bin <- reclassify(mh_present[[i]], c(mh_poligono,Inf,NA))
  mh_present_umbral <- raster::stack(mh_present_umbral, mh_present_bin)
  mh_futuro_bin <- reclassify(mh_futuro[[i]], c(mh_poligono,Inf,NA))
  mh_futuro_umbral <- raster::stack(mh_futuro_umbral, mh_futuro_bin)
}

# Export raster
for ( i in 1:nlayers(mh_present_umbral)){
  writeRaster(mh_present_umbral[[i]], paste0("T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Resultados/2/mh_present_50_85_", names[i], "_T.tif"), overwrite=TRUE)
  writeRaster( mh_futuro_umbral[[i]], paste0( "T:/MODCLIM_R_DATA/DATOS_ANTIGUOS/Resultados/2/mh_future_50_85_", names[i],   "_T.tif"), overwrite=TRUE)
}








#############################################################################################
############################################################################################
###########################################################################################


