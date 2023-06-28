library(sf)
library(raster)
# Load study area shp
study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/CUENCAS", "\\.shp$", full.names = T))

# Reference system
study_area <- st_transform(study_area, 
                           crs("+proj=longlat +datum=WGS84 +no_defs"))


# Period selection
  # 1: 2011-2040
  # 2: 2041-2070
  # 3: 2071-2100
p <- 1

# Model selection
  # 1: GFDL-ESM4 (National Oceanic and Atmospheric Administration, Geophysical Fluid Dynamics Laboratory, Princeton, NJ 08540, USA)
  # 2: IPSL-CM6A-LR (Institut Pierre Simon Laplace, Paris 75252, France)
  # 3: MPI-ESM1-2-HR (Max Planck Institute for Meteorology, Hamburg 20146, Germany)
  # 4: MRI-ESM2-0 (Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan)
  # 5: UKESM1-0-LL (Met Office Hadley Centre, Fitzroy Road, Exeter, Devon, EX1 3PB, UK)
m <- 2

# Scenario selection
  # 1: ssp126 (SSP1-RCP 2.6 climate as simulated by the GCMs)
  # 2: ssp370 (SSP3-RCP 7 climate as simulated by the GCMs)
  # 3: ssp585 (SSP5-RCP 8.5 climate as simulated by the GCMs)
s <- 1

period <- c("2011-2040", "2041-2070", "2071-2100")
model <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
scenario <- c("ssp126", "ssp370", "ssp585")


period   <- period[p]
model    <- model[m]
scenario <- scenario[s]


download_path <- "T:/CHELSA_FUTURE/"

chelsa_future_tmax(period, model, scenario, download_path)

# Monthly precipitation amount

chelsa_future_pcp <- function(period, model, scenario, download_path){
  data <- readLines("https://raw.githubusercontent.com/MarioMingarro/climatic_representativeness/main/pcp_future_download_path.txt")
  data_pediod <- str_subset(data, pattern = paste0(period))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario))
  for(i in 1:length(data_pediod_model_scenario)){
    dir.create(paste0(download_path, "pcp"))
    dir.create(paste0(download_path, "/pcp/", period, "_", model, "_", scenario))
    directory <- paste0(download_path, "pcp/", period, "_", model, "_", scenario, "/")
    url <- data_pediod_model_scenario[i]
    name <- paste0("pr_", str_sub(sub(".*pr_", "", data_pediod_model_scenario[i]),1,12), ".tif")
    download.file(url, destfile = paste0(directory,"raster.tif"), mode="wb")
    raster <- raster(paste0(directory,"raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    raster <- raster/100
    writeRaster(raster, paste0(directory,name))
  }
}

# Mean daily maximum air temperature 

chelsa_future_tmax <- function(period, model, scenario, download_path){
  data <- readLines("https://raw.githubusercontent.com/MarioMingarro/climatic_representativeness/main/tmax_future_download_path.txt")
  data_pediod <- str_subset(data, pattern = paste0(period))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario))
  for(i in 1:length(data_pediod_model_scenario)){
    dir.create(paste0(download_path, "tmax"))
    dir.create(paste0(download_path, "/tmax/", period, "_", model, "_", scenario))
    directory <- paste0(download_path, "tmax/", period, "_", model, "_", scenario, "/")
    url <- data_pediod_model_scenario[i]
    name <- paste0("tmax_", str_sub(sub(".*tasmax_", "", data_pediod_model_scenario[i]),1,12), ".tif")
    download.file(url, destfile = paste0(directory,"raster.tif"), mode="wb")
    raster <- raster(paste0(directory,"raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    writeRaster(raster, paste0(directory,name))
  }
}

# Mean daily minimum air temperature 

chelsa_future_tmax <- function(period, model, scenario, download_path){
  data <- readLines("https://raw.githubusercontent.com/MarioMingarro/climatic_representativeness/main/tmin_future_download_path.txt")
  data_pediod <- str_subset(data, pattern = paste0(period))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario))
  for(i in 1:length(data_pediod_model_scenario)){
    dir.create(paste0(download_path, "tmin"))
    dir.create(paste0(download_path, "/tmin/", period, "_", model, "_", scenario))
    directory <- paste0(download_path, "tmin/", period, "_", model, "_", scenario, "/")
    url <- data_pediod_model_scenario[i]
    name <- paste0("tmin_", str_sub(sub(".*tasmin_", "", data_pediod_model_scenario[i]),1,12), ".tif")
    download.file(url, destfile = paste0(directory,"raster.tif"), mode="wb")
    raster <- raster(paste0(directory,"raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    writeRaster(raster, paste0(directory,name))
  }
}
