library(raster)
library(tidyverse)
library(sf)
#CRU
# Download: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/

monthly_tmin <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.tmn.dat.nc"), varname = "tmn")
monthly_tmax <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.tmx.dat.nc"), varname = "tmx")
monthly_pcp <- raster::brick(here::here("T:/MODCLIM_R_DATA/CRU/cru_ts4.07.1901.2022.pre.dat.nc"), varname = "pre")

res(monthly_tmin)

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))

reference_system <- "+proj=longlat +datum=WGS84 +no_defs"
study_area <- st_transform(study_area, crs(reference_system))

beginCluster()

monthly_tmin <- mask(crop(monthly_tmin, study_area), study_area)
monthly_tmax <- mask(crop(monthly_tmax, study_area), study_area)
monthly_pcp <- mask(crop(monthly_pcp, study_area), study_area)

endCluster()


### Monthly data to annual average ----

year <- seq(1901, 2022, 1)
### Maximum temperature ----
beginCluster()
annual_tmax <- raster::stack()

for (i in 1:length(year)){
  raster <- calc(raster::subset(monthly_tmax, grep(paste0(year[i]), names(monthly_tmax), value = T)), mean) # MEAN
  annual_tmax <- raster::stack(annual_tmax, raster)
}
names(annual_tmax) <- paste0("Y_", seq(1901, 2022, 1))
endCluster()

### Minimum temperature ----
beginCluster()
annual_tmin <- raster::stack()

for (i in 1:length(year)){
  raster <- calc(raster::subset(monthly_tmin, grep(paste0(year[i]), names(monthly_tmin), value = T)), mean) # MEAN
  annual_tmin <- raster::stack(annual_tmin, raster)
}
names(annual_tmin) <- paste0("Y_", seq(1901, 2022, 1))
endCluster()

### Precipitation ----
beginCluster()
annual_pcp <- raster::stack()

for (i in 1:length(year)){
  raster <- calc(raster::subset(monthly_pcp, grep(paste0(year[i]), names(monthly_pcp), value = T)), sum) # MEAN
  annual_pcp <- raster::stack(annual_pcp, raster)
}
names(annual_pcp) <- paste0("Y_", seq(1901, 2022, 1))
endCluster()

# Bioclimatic creation ----

# BIO1 = Annual Mean Temperature
beginCluster()
BIO1 <- calc(raster::stack(annual_tmax, annual_tmin), mean)
endCluster()
# BIO2 = Mean Diurnal Range (Mean of monthly * (max temp - min temp))
beginCluster()
BIO2 <- raster::brick()

for (i in 1:length(year)){
  max_monthly_tmax <- calc(monthly_tmax, max)
  min_monthly_tmin <- calc(monthly_tmin, min)
  monthly_tmed <- stack(monthly_tmax, monthly_tmin)
  monthly_tmed <- calc(monthly_tmed, mean)
  raster <- monthly_tmed*(max_monthly_tmax - min_monthly_tmin)
  BIO2 <- raster::stack(BIO2, raster)
}

BIO2 <- calc(BIO2, mean)
endCluster()
plot(BIO2)

# BIO3 = Isothermality (BIO2/BIO7) (×100)
BIO3 <- BIO2/BIO7
# BIO4 = Temperature Seasonality (standard deviation ×100)
# 
# BIO5 = Max Temperature of Warmest Month
BIO5 <- calc(monthly_tmax, mean)

# BIO6 = Min Temperature of Coldest Month
BIO6 <- calc(monthly_tmin, mean)

# BIO7 = Temperature Annual Range (BIO5-BIO6)
BIO7 <- BIO5-BIO6

# BIO8 = Mean Temperature of Wettest Quarter
# 
# BIO9 = Mean Temperature of Driest Quarter
# 
# BIO10 = Mean Temperature of Warmest Quarter
# 
# BIO11 = Mean Temperature of Coldest Quarter
# 
# BIO12 = Annual Precipitation
BIO12 <- calc(PCP_annual, mean)

# BIO13 = Precipitation of Wettest Month
# 
# BIO14 = Precipitation of Driest Month
# 
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# 
# BIO16 = Precipitation of Wettest Quarter
# 
# BIO17 = Precipitation of Driest Quarter
# 
# BIO18 = Precipitation of Warmest Quarter
# 
# BIO19 = Precipitation of Coldest Quarter


