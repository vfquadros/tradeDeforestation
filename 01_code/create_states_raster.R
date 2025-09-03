# use R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
rm(list = ls())
gc()

Sys.getlocale()
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# PREAMBLE ---------------------------------------------------------------------

# Loading libraries
library(sf)

crs <- "epsg:4326"

# Worling directory
setwd("C:/Users/vquadros/Dropbox/tradeDeforestation/theory/tradeDeforestation")

state_boundaries <- st_read("02_inputs/dashboard_states-static-layer.shp")
legal_amazon <- st_read("02_inputs/brazilian_legal_amazon.shp")

final_shp <- st_intersection(legal_amazon, state_boundaries)

st_write(final_shp, "02_inputs/states_legal_amazon.shp")

grid <- rast(ext(-73.99097, -43.95183, -18.0417, 5.272225), resolution = 0.01, crs = "epsg:4326") %>% 
  setValues(1)

# crop the 0.01deg cell to obtain those only in Brazil
template <- mask(grid, vect(final_shp))

r <- rasterizeGeom(vect(final_shp),template,"count")
