# --------------- Task --------------- #
# Prepare rasters to run deforestation model in matlab. As of now, we have to
# prepare 3 data sets. One for productivity, one for initial deforestation and
# agribusiness and the last one for protected areas.
#
# Productivity and deforestation are rasters and protected areas in vectorial.
# This script prepare rasters with resolution 0.01 0.01 which is the resolution
# of the deforestation data (prepared by Jialing in previous work). We also use 
# PRODES deforestation data, which is vectorial as a small extension.
# ------------------------------------ #

# use R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
rm(list = ls())
gc()

Sys.getlocale()
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# PREAMBLE ---------------------------------------------------------------------

# Loading libraries
library(sf)
library(dplyr)
library(ggplot2)
library(raster)
library(tiff)
library(terra)

# Worling directory
setwd("C:/Users/vquadros/Dropbox/tradeDeforestation/theory/simulation_exercises")

# STEP 0 - Crop the state of Mato Grosso ---------------------------------------
# Note for the future: We  will do this to the whole legal amazon territory at 
# some point)

file <- read_sf("02_inputs/dashboard_states-static-layer.shp")
# Check CRS
if (is.na(st_crs(file)$epsg) || st_crs(file)$epsg != 4326) {
  file <- st_transform(file, 4326)
}

# Filtering Mato Grosso's geometry
MT <- file %>% filter(name == "MATO GROSSO")

grid <- rast(ext(-61.2, -50.2, -18, -7.3), resolution = 0.01, crs = "epsg:4326") %>% 
  setValues(1)

# crop the 0.01deg cell to obtain those only in Brazil
template <- mask(grid, vect(MT))
#plot(template)

# STEP 1 - Mapbiomas deforestastion data ---------------------------------------

deforestation <- read_sf("02_inputs/sf_file_2002.gpkg") %>% filter(land_use %in% c("pasture","agri"))
deforestation <- st_intersection(deforestation, MT)

deforestation_raster <- rasterize(deforestation, template, field = 1, fun = "sum", cover = TRUE, background = 0)

writeRaster(deforestation_raster, "03_outputs/exercise_mt/deforestation001.tif", overwrite=TRUE)

# STEP 2 - Soy productivity data -----------------------------------------------

# Load
productivity <- rast("02_inputs/High_soyb200b_yld.tif")

# Crop
productivity <- crop(productivity, st_bbox(MT))

# Mask
productivity <- mask(productivity, MT)

# Resample
productivity_raster <- resample(productivity, deforestation_raster, method = "near")

writeRaster(productivity_raster, "03_outputs/exercise_mt/productivity001.tif", overwrite=TRUE)
# STEP 3 - Legal Status data ---------------------------------------------------

# Loading indiginous sites and conservation units and intersecting with MT 
conservation_units <- st_read("02_inputs/conservation_units_legal_amazon.shp")
indiginous_sites <- st_read("02_inputs/indigenous_area_legal_amazon.shp")

# CRS sanity check
if (is.na(st_crs(conservation_units)$epsg) || st_crs(conservation_units)$epsg != 4326) {
  conservation_units <- st_transform(conservation_units, st_crs(MT))
}
if (is.na(st_crs(indiginous_sites)$epsg) || st_crs(indiginous_sites)$epsg != 4326) {
  indiginous_sites <- st_transform(indiginous_sites, st_crs(MT))
}

# Some polygon here has a duplicated vertex, so we use this function
conservation_units <- st_make_valid(conservation_units)

conservation_units <- st_intersection(conservation_units,MT)
indiginous_sites <- st_intersection(indiginous_sites,MT)

protected <- st_union(conservation_units,indiginous_sites)

save(protected, file = "03_outputs/protected_areas.RData")
protected_raster <- rasterize(protected, template, field = 1, fun = "sum", cover = TRUE, background = 0)

writeRaster(protected_raster, "03_outputs/exercise_mt/protected001.tif", overwrite=TRUE)

