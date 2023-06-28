# Period selection
p <- 
# Model selection
m <- 
# scenario selection
s <- 

period <- c("2011-2040", "2041-2070", "2071-2100")
model <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
scenario <- c("ssp126", "ssp370", "ssp585")


period   <- period[p]
model    <- model[m] 
scenario <- scenario[s]


download_path <- "T:/CHELSA_FUTURE/"
data_path <- "T:/GITHUB_REP/climatic_representativeness/"


chelsa_future_precipitation <- function(variable, period, model, scenario, download_path){
  data <- readLines(paste0(data_path, "precipitation_download_path.txt"))
  data_pediod <- str_subset(data, pattern = paste0(period[1]))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model[1]))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario[1]))
  dir.create(paste0(download_path, "precipitation"))
  dir.create(paste0(download_path, "/precipitation/", period[1], "_", model[1], "_", scenario[1]))
  download.file(data_pediod_model_scenario, download_path)
}

chelsa_future_maximum_temperature <- function(variable, period, model, scenario, download_path){
  data <- readLines(paste0(data_path, "precipitation_download_path.txt"))
  data_pediod <- str_subset(data, pattern = paste0(period[1]))
  data_pediod_model <- str_subset(data_pediod, pattern = paste0(model[1]))
  data_pediod_model_scenario <- str_subset(data_pediod_model, pattern = paste0(scenario[1]))
  dir.create(paste0(download_path, "precipitation"))
  dir.create(paste0(download_path, "/precipitation/", period[1], "_", model[1], "_", scenario[1]))
  download.file(data_pediod_model_scenario, download_path)
}
