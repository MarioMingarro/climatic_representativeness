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
present_climatic_variables <- raster::stack(list.files("D:/MODCLIM/CLIMA/PRESENTE/", "\\.tif$", full.names = T))

future_climatic_variables <- raster::stack(list.files("D:/MODCLIM/CLIMA/FUTURO/IPSL/", "\\.tif$", full.names = T))

names(present_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                       "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7","CHELSA_bio8","CHELSA_bio9")
names(future_climatic_variables) <- c("CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12","CHELSA_bio13","CHELSA_bio14",
                                       "CHELSA_bio15","CHELSA_bio16","CHELSA_bio17","CHELSA_bio18","CHELSA_bio19","CHELSA_bio2",
                                      "CHELSA_bio3","CHELSA_bio4","CHELSA_bio5","CHELSA_bio6","CHELSA_bio7","CHELSA_bio8","CHELSA_bio9")

study_area <- read_sf("D:/MODCLIM/CHILE/CHILE.shp")

polygon <- read_sf("D:/MODCLIM/CHILE/PN_Vicente_Perez.shp")


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
colnames(data_present_climatic_variables) <- c("x","y","CHELSA_bio1","CHELSA_bio10","CHELSA_bio11","CHELSA_bio12",
                                               "CHELSA_bio13","CHELSA_bio14","CHELSA_bio15","CHELSA_bio16","CHELSA_bio17",
                                               "CHELSA_bio18","CHELSA_bio19","CHELSA_bio2","CHELSA_bio3","CHELSA_bio4",
                                               "CHELSA_bio5","CHELSA_bio6","CHELSA_bio7","CHELSA_bio8","CHELSA_bio9")

colnames(data_future_climatic_variables) <- colnames(data_present_climatic_variables)



# Correlation between variables ----
cor <- cor(data_present_climatic_variables[,3:21])

# Select variables less correlated 
drop_1  <-  findCorrelation(cor, cutoff = .6)
drop  <-  names(data_present_climatic_variables[,3:21])[drop_1]
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
data <- data_present_climatic_variables

xxxxxxxx

# Create name object
names <- polygon$WDPAID


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
mh_present <- raster::brick()

for(j in 4:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Present")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_present <- raster::stack(mh_present, mh_f)
}



mh_future <- raster::brick()

for(j in 4:length(mh)){
  mh_f <- dplyr::filter(mh, Period == "Future")
  mh_f <- rasterFromXYZ(mh_f[, c(1:2,j)])
  names(mh_f) <- colnames(mh[j])
  mh_future <- raster::stack(mh_future, mh_f)
}

plot(mh_future)
# Export raster

#for ( i in 1:nlayers(mh_present)){
writeRaster(mh_present[[i]], paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/Slovenia/mh_present_IPSL_2040_2070_SSP85", names[i], ".tif"), overwrite=TRUE)
writeRaster( mh_future[[i]], paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/Slovenia/mh_future_IPSL_2040_2070_SSP85", names[i],   ".tif"), overwrite=TRUE)
#}
writeRaster(mh_present_umbral[[1]], "D:/MODCLIM/RESULT/mh_pre_madeira_umbral.tif")
writeRaster( mh_future_umbral[[1]], "D:/MODCLIM/RESULT/mh_fut_madeira_umbral.tif")

writeRaster(mh_present_umbral[[1]], "D:/MODCLIM/CHILE/PN_P_UMB.tif")
writeRaster(mh_future_umbral[[1]], "D:/MODCLIM/CHILE/PN_F_UMB.tif")
writeRaster(mh_present[[1]], "D:/MODCLIM/CHILE/PN_P.tif")
writeRaster(mh_future[[1]], "D:/MODCLIM/CHILE/PN_F.tif")

# Threshold selection
mh_present_umbral <- raster::brick()
mh_future_umbral <- raster::brick()

for (i in 1:nlayers(mh_present)){
  mh_polygon <- mask(crop(mh_present[[i]], polygon[i,]), polygon[i,])
  mh_polygon <- raster::as.data.frame(mh_polygon, xy = T)
  mh_polygon <- quantile(na.omit(mh_polygon[,3]), probs = c(.90))
  mh_present_bin <- reclassify(mh_present[[i]], c(mh_polygon,Inf,NA))
  mh_present_umbral <- raster::stack(mh_present_umbral, mh_present_bin)
  mh_future_bin <- reclassify(mh_future[[i]], c(mh_polygon,Inf,NA))
  mh_future_umbral <- raster::stack(mh_future_umbral, mh_future_bin)
}


# Export raster
#for ( i in 1:nlayers(mh_present_umbral)){
#  writeRaster(mh_present_umbral[[i]],  paste0("D:/REPRESENATTIVENESS/RAMSAR_RESULTS/mh_present_IPSL_2040_2070_SSP85", names[i], "_T.tif"), overwrite=TRUE)
#  writeRaster( mh_future_umbral[[i]],  paste0("D:/REPRESENATTIVENESS/RAMSAR_RESULTS/mh_future_IPSL_2040_2070_SSP85", names[i],   "_T.tif"), overwrite=TRUE)
#}




rm(list= ls()[!(ls() %in% c("mh_present_umbral", "mh_future_umbral", "polygon", "study_area", "cor", "drop", "drop_1"))])
rm(cor)
rm(drop)
rm(drop_1)

gc(reset=TRUE)


reference_system <- projection("+init=epsg:25828") # "+proj=longlat +datum=WGS84 +no_defs"   

projection(mh_present_umbral) = CRS("+proj=longlat +datum=WGS84 +no_defs")
projection(mh_future_umbral) = CRS("+proj=longlat +datum=WGS84 +no_defs")
mh_present_umbral <- projectRaster(mh_present_umbral, crs = reference_system)
mh_future_umbral <- projectRaster(mh_future_umbral, crs = reference_system)
study_area <- st_transform(study_area, crs(reference_system))
polygon <- st_transform(polygon, crs(reference_system))

isolation <- data.frame(name=character(), 
                        present_dist=numeric(), 
                        future_present_dist=numeric(), 
                        future_dist=numeric())



for (i in 1:nlayers(mh_present_umbral)){
  aa <- data.frame(name="a", 
                   present_dist=2,
                   future_present_dist=2,
                   future_dist=2)
  pol <- polygon[i,]
  pre_in <- raster::mask(mh_present_umbral[[i]], pol)
  n_pre_in <- length(na.omit(pre_in@data@values))
  pre_out <- raster::mask(mh_present_umbral[[i]], pol, inverse = T)
  fut_out <- raster::mask(mh_future_umbral[[i]], pol, inverse = T)
  pre_out <- rasterToPoints(pre_out, spatial = T)
  fut_out <- rasterToPoints(fut_out, spatial = T)
  
  pol <- as(pol, 'Spatial')
  distance_pre <- sort(rgeos::gDistance(pre_out, pol, byid=TRUE))
  distance_fut <- sort(rgeos::gDistance(fut_out, pol, byid=TRUE))
  
  aa$name <- pol$NAME
  aa$present_dist <-  mean(distance_pre[1:n_pre_in])
  aa$future_present_dist <-  mean(distance_fut[1:n_pre_in])
  aa$future_dist <-  mean(distance_fut)
  isolation <- rbind(isolation, aa)
}

toc()


writexl::write_xlsx(isolation, "D:/REPRESENATTIVENESS/isolation.xlsx")

mean(distance_fut)

proj4string(pol)
projection(pre_out)
i=3
plot(mh_present_umbral[[1]])
plot(mh_future_umbral[[1]])

#pre_out <- SpatialPoints(pre_out[,1:2])
#fut_out <- SpatialPoints(fut_out[,1:2])


#GEODIVERSE---- 

geomorphologic_variables <-  raster::stack(list.files("D:/MODCLIM/GEODIVERSIDAD/Macaronesia/reclass/", ".tif", full.names = T))
geomorphologic_variables <- projectRaster(geomorphologic_variables, crs = reference_system)
plot(geomorphologic_variables)


geomorphologic_variables <- raster::aggregate(geomorphologic_variables, 4)


# Extract raster data
data_geo <- raster::as.data.frame(geomorphologic_variables, xy = TRUE)

# Delete NA
data_geo<-na.omit(data_geo)


###################################################
# Correlation between variables ----
cor <- cor(data_geo[,3:7])

#Select variables less correlated 
drop_1  <-  findCorrelation(cor, cutoff = .8)
drop  <-  names(data_geo[,3:7])[drop_1]
data_geo <- data_geo[!names(data_geo) %in% drop]

geo <- dropLayer(geomorphologic_variables, names(geomorphologic_variables)[drop_1])



corrplot(cor(data_geo[4:length(data_geo)]),
         method = "number",
         type = "upper")


# Add field period 
data_geo <- mutate(data_geo, Period = c("PN"),  .after = "y")





mh <- mahalanobis(data[,4:length(data)], 
                  colMeans(data_geo[,4:length(data_geo)]), 
                  cov(data[,4:length(data)]), 
                  inverted = F)



# Agregar informacion espacial al reusltado de mh
mh2 <- cbind(data[,c(1:3)], mh)

mh_pre <- dplyr::filter(mh2, Period == "PN")
mh_pre <- rasterFromXYZ(mh_pre[, c(1,2,4)])

mh_fut <- dplyr::filter(mh2, Period == "OUT")
mh_fut <- rasterFromXYZ(mh_fut[, c(1,2,4)])

writeRaster(mh_pre, "T:/MODCLIM_R_DATA/ANALISIS/SLOVENIA/mh_geo_pre.tif",overwrite=TRUE)
writeRaster(mh_fut, "T:/MODCLIM_R_DATA/ANALISIS/SLOVENIA/mh_geo_fut.tif",overwrite=TRUE)






polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/PAS/national_parks_slovenia.shp")
if(compareCRS(polygon, geo) == TRUE) {
  print("Same Reference System")
} else {
  polygon <- st_transform(polygon, projection(geo))
}

  mh_polygon <- mask(crop(mh_pre, polygon), polygon)
  mh_polygon <- raster::as.data.frame(mh_polygon, xy = T)
  mh_polygon <- quantile(na.omit(mh_pre[,3]), probs = c(.90))
  mh_pre_bin <- reclassify(mh_pre, c(mh_polygon,Inf,NA))
  mh_fut_bin <- reclassify(mh_fut, c(mh_polygon,Inf,NA))

  writeRaster(mh_pre, "T:/MODCLIM_R_DATA/ANALISIS/SLOVENIA/mh_geo_pre_Threshold.tif",overwrite=TRUE)
  writeRaster(mh_fut, "T:/MODCLIM_R_DATA/ANALISIS/SLOVENIA/mh_geo_fut_Threshold.tif",overwrite=TRUE)
  
# Export raster
for ( i in 4){
  writeRaster(mh_present_umbral[[i]],  paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/Slovenia/RECIPIENT/mh_present_IPSL_2040_2070_SSP85", names[i], "_T.tif"), overwrite=TRUE)
  writeRaster( mh_future_umbral[[i]],  paste0("T:/MODCLIM_R_DATA/ANALISIS/RESULTADOS/Slovenia/RECIPIENT/mh_future_IPSL_2040_2070_SSP85", names[i],   "_T.tif"), overwrite=TRUE)
}






















# PCA ----
# Climatic representativeness
PCA_variables <- RStoolbox::rasterPCA(variables, spca = TRUE)


# PCA results
summary(PCA_variables$model)


round(PCA_variables$model$loadings[,1:6],3)


selected_PCA_variables <- PCA_variables$map[[1:6]]


plot(selected_PCA_variables)

data_PCA <- raster::as.data.frame(variables, xy = TRUE)
data_PCA <- na.omit(data_PCA)

presente <- data_PCA[,1:21]
futuro <- data_PCA[,c(1:2, 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)]

geo <- data_PCA[,c(1:2, 41,42,43,44,45,46)]


presente <- mutate(presente, Periodo = c("Presente"))
futuro <- mutate(futuro, Periodo = c("Futuro"))

presente <- left_join(presente, geo, by = c("x", "y"))
futuro <- left_join(futuro, geo, by = c("x", "y"))

colnames(futuro) <- colnames(presente)
data <- rbind(presente, futuro)




writeRaster(selected_PCA_present_climatic_variables, "T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")




###########################################################################################################
####################################################################################

PCA_present <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_present.tif ")
PCA_future <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_future.tif ")
PCA_geo <- raster::stack("T:/MODCLIM_R_DATA/ANALISIS/PCA/PCA_geo.tif ")


#PCA_geo <- resample(PCA_geo, PCA_present, method='bilinear')
#
#delete <- c("PCA_present", "PCA_future", "PCA_geo", "study_area", "polygon")
#rm(list=(ls()[ls()!= delete]))

# Load data
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/ANALISIS", "\\.shp$", full.names = T))

polygon <- read_sf("T:/MODCLIM_R_DATA/ANALISIS/NP/national_parks.shp")



data_PCA <- raster::as.data.frame(selected_PCA_variables, xy = TRUE)
data_PCA_present <- na.omit(data_PCA_present)
data_PCA_present <- mutate(data_PCA_present, Periodo = c("Presente"))
data_PCA_present <- mutate(data_PCA_present, Var = c("climatic"))

data_PCA_future <- raster::as.data.frame(PCA_future, xy = TRUE)
data_PCA_future <- na.omit(data_PCA_future)
data_PCA_future <- mutate(data_PCA_future, Periodo = c("Futuro"))
data_PCA_future <- mutate(data_PCA_future, Var = c("climatic"))

data_PCA_geo <- raster::as.data.frame(PCA_geo, xy = TRUE)
data_PCA_geo <- na.omit(data_PCA_geo)
data_PCA_geo <- mutate(data_PCA_geo, Periodo = c("Presente"))
data_PCA_geo <- mutate(data_PCA_geo, Var = c("geomor"))


colnames(data_PCA_present) <- c("x", "y", "PC_1", "PC_2", "PC_3", "PC_4", "Period", "Var")
colnames(data_PCA_geo) <- colnames(data_PCA_present)
colnames(data_PCA_future) <- colnames(data_PCA_present)


data <- rbind(data_PCA_present, data_PCA_geo)
data <- rbind(data, data_PCA_future)




