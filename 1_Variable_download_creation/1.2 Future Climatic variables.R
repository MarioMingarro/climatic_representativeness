library(sf)
library(raster)

source("T:/GITHUB_REP/climatic_representativeness/download_climatic_functions.R")
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
download_path <- "T:/MODCLIM_R_DATA/DESCARGA"

chelsa_future_pcp(period, model, scenario, download_path)

