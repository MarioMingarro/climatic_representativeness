packages.to.use <- c("raster", ,"doParallel","foreach","ggplot2", 
                     "sf", "rnaturalearth", "rnaturalearthdata", 
                     "viridis", "dismo")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages(package) }
  #if( ! package %in% rownames(installed.packages()) & package == "VoCC" ) { devtools::install_github("JorGarMol/VoCC", dependencies = TRUE) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  suppressWarnings( library(package, character.only = TRUE) )
}



