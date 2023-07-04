# Download current chelsa predictions
# https://chelsa-climate.org/timeseries/

study_area <- read_sf(list.files("T:/MODCLIM_R_DATA/analysis/", "\\.shp$", full.names = T))
download_path <- "T:/CHELSA_PRESENTE/"


# Monthly precipitation amount

chelsa_present_pcp <- function(download_path){
  data <- readLines("pcp_present_1980_2018_download_path.txt")
  dir.create(paste0(download_path, "pcp"))
  directory <- paste0(download_path, "pcp/")
  for(i in 1:length(data_pediod_model_scenario)){
    url <- data[i]
    name <- paste0("pr_", str_sub(sub(".*pr_", "", data[i]),1,7), ".tif")
    download.file(url, destfile = paste0(directory,"raster.tif"), mode="wb")
    raster <- raster(paste0(directory,"raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    raster <- raster/100
    writeRaster(raster, paste0(directory,name))
  }
}


library(doParallel)
library(foreach)

chelsa_present_pcp <- function(download_path, study_area) {
  data <- readLines("pcp_present_1980_2018_download_path.txt")
  dir.create(paste0(download_path, "pcp"))
  directory <- paste0(download_path, "pcp/")
  
  # Set up parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(data), .packages = c("raster")) %dopar% {
    url <- data[i]
    name <- paste0("pr_", str_sub(sub(".*pr_", "", data[i]), 1, 7), ".tif")
    download.file(url, destfile = paste0(directory, "raster.tif"), mode = "wb")
    raster <- raster(paste0(directory, "raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    raster <- raster / 100
    writeRaster(raster, paste0(directory, name))
  }
  
  stopCluster(cl)
}


chelsa_present_pcp(download_path, study_area)







# Download future chelsa predictions
# https://chelsa-climate.org/cmip6/

# Monthly precipitation amount

library(doParallel)
library(foreach)

chelsa_present_pcp <- function(download_path) {
  data <- readLines("pcp_present_1980_2018_download_path.txt")
  dir.create(paste0(download_path, "pcp"))
  directory <- paste0(download_path, "pcp/")
  
  # Set up parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(data), .packages = c("raster")) %dopar% {
    url <- data[i]
    name <- paste0("pr_", str_sub(sub(".*pr_", "", data[i]), 1, 7), ".tif")
    download.file(url, destfile = paste0(directory, "raster.tif"), mode = "wb")
    raster <- raster(paste0(directory, "raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    raster <- raster / 100
    writeRaster(raster, paste0(directory, name))
  }
  
  stopCluster(cl)
}




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



library(doParallel)
library(foreach)

chelsa_future_pcp <- function(period, model, scenario, download_path) {
  data <- readLines("https://github.com/MarioMingarro/climatic_representativeness/blob/main/Download_paths/pcp_future_download_path.txt")
  data_pediod <- str_subset(data, pattern = paste0(period))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario))
  
  # Set up parallel processing
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(data_pediod_model_scenario), .packages = c("raster", "stringr")) %dopar% {
    dir.create(paste0(download_path, "pcp"), showWarnings = FALSE)
    dir.create(paste0(download_path, "/pcp/", period, "_", model, "_", scenario), showWarnings = FALSE)
    directory <- paste0(download_path, "pcp/", period, "_", model, "_", scenario, "/")
    url <- data_pediod_model_scenario[i]
    name <- paste0("pr_", str_sub(sub(".*pr_", "", data_pediod_model_scenario[i]), 1, 12), ".tif")
    download.file(url, destfile = paste0(directory, "raster.tif"), mode = "wb")
    raster <- raster(paste0(directory, "raster.tif"))
    raster <- raster::mask(crop(raster, study_area), study_area)
    raster <- raster / 100
    writeRaster(raster, paste0(directory, name))
  }
  
  stopCluster(cl)
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
