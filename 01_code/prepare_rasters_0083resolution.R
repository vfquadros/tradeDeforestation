# --------------- Task --------------- #
# Replicate Jialings work to MT
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
setwd("C:/Users/vquadros/Desktop/amzlgl/data")

# STEP 2 - Obtain the productivity raster --------------------------------------
# Note: this will be the most agregated raster

# read the Brazil admin1 level geometry shp file
file <- read_sf(file.path(wd,"dashboard_states-static-layer/dashboard_states-static-layer.shp"))

# If the EPSG code is missing or not 4326, transform the CRS to WGS 84 (EPSG: 4326)
if (is.na(st_crs(file)$epsg) || st_crs(file)$epsg != 4326) {
  file <- st_transform(file, 4326)
}

# Filtering one particular state, MT. Afterwards we can do the whole country.
MT <- file %>% filter(name == "MATO GROSSO")
#plot(MT)

# Gets soy productivity data. This will become a more complicated process when
# we introduce more crops. Also crops and masks the raster to MT
raster_data <- raster(file.path(wd,"High_soyb200b_yld.tif"))
productivity <- crop(raster_data, st_bbox(MT))
productivity <- mask(raster_data)

writeRaster(productivity,file.path(wd,"outputs/productivity.tif"), overwrite = TRUE)

template <- rast(productivity)

# STEP 3 - Legal status raster -------------------------------------------------

protected <- read_sf(file.path(wd,"protected_areas/conservation_units_legal_amazon.shp"))

# Template = same grid as 'productivity'
template <- rast(productivity)

# Ensure polygons are a SpatVector in the same CRS
vprot <- vect(protected)
if (!same.crs(template, vprot)) vprot <- project(vprot, template)

# 2) Fractional cover of each cell by protected polygons (0..1)
#    'field=1' with fun='sum' + cover=TRUE returns the covered fraction (after dissolve)
cover <- rasterize(vprot, template, field = 1, fun = "sum", cover = TRUE, background = 0)

# 3) Majority rule: 1 if >= 50% covered, else 0
prot_major <- ifel(cover >= 0.5, 1, 0)

# 4) Keep NA pattern identical to 'productivity' (optional but common)

# 5) Save + plot
#writeRaster(prot_major, "protected_majority.tif", overwrite = TRUE, datatype = "INT1U")
#plot(prot_major, main = "Protected (majority rule, ≥50%)")

writeRaster(prot_major, file.path(wd,"outputs/protected.tif"), overwrite=TRUE)

# STEP 3 - MapBiomes Deforestation ---------------------------------------------

deforestation <- read_sf(file.path(wd,"sf_file_2002.gpkg"))

deforestation <- st_intersection(deforestation, MT)
deforestation <- deforestation %>% filter(land_use %in% c("pasture","agri"))
# 2) Fractional cover of each cell by protected polygons (0..1)
#    'field=1' with fun='sum' + cover=TRUE returns the covered fraction (after dissolve)
cover <- rasterize(deforestation, template, field = 1, fun = "sum", cover = TRUE, background = 0)

# 3) Majority rule: 1 if >= 50% covered, else 0
deforestation_major <- ifel(cover >= 0.5, 1, 0)

# 4) Keep NA pattern identical to 'productivity' (optional but common)

# 5) Save + plot
#writeRaster(prot_major, "protected_majority.tif", overwrite = TRUE, datatype = "INT1U")
plot(deforestation_major, main = "Protected (majority rule, ≥50%)")

writeRaster(deforestation_major, file.path(wd,"outputs/deforestation.tif"), overwrite=TRUE)

