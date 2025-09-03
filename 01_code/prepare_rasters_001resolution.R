# ------------------------------- Task --------------------------------------- #
# Prepare rasters to run deforestation model in matlab. As of now, we have to
# prepare 3 data sets. One for productivity, one for initial deforestation and
# agribusiness and the last one for protected areas.
#
# Productivity and deforestation are rasters and protected areas in vectorial.
# This script prepare rasters with resolution 0.01 0.01 which is the resolution
# of the deforestation data (prepared by Jialing in previous work).
# ---------------------------------------------------------------------------- #

# use R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
rm(list = ls())
gc()

Sys.getlocale()
Sys.setlocale("LC_ALL", "en_US.UTF-8")

# PREAMBLE ---------------------------------------------------------------------

# Loading libraries
library(sf)
library(parallel)
library(stars)
library(dplyr)
library(ggplot2)
library(raster)
library(tiff)
library(terra)
library(exactextractr)
library(tidyverse)
library(ggmap)

# Worling directory
setwd("C:/Users/vquadros/Dropbox/tradeDeforestation/theory/tradeDeforestation")

# STEP 0 - CREATE MODEL GRID ----------- ---------------------------------------
# Note for the future: We  will do this to the whole legal amazon territory at 
# some point)

file <- read_sf("02_inputs/brazilian_legal_amazon.shp")

# Check CRS
if (is.na(st_crs(file)$epsg) || st_crs(file)$epsg != 4326) {
  file <- st_transform(file, 4326)
}

# these are the limits of the borders of the brazilian legal amazon.
grid <- rast(ext(-73.99097, -43.95183, -18.0417, 5.272225), resolution = 0.01, crs = "epsg:4326") %>% 
  setValues(1)

# crop the 0.01deg cell to obtain those only in Brazil
template <- mask(grid, vect(file))
#plot(template)

writeRaster(template, "03_outputs/exercise_amabiome/mask001.tif", overwrite=TRUE)

# STEP 1 - MapBiomes Agribusiness ----------------------------------------------

# creating a vector version of template

final_grid_0_01 <- st_read("02_inputs/final_grid_0_01.gpkg")

dfs <- list()
for (i in 1:8){
  load(paste0("02_inputs/landuse_allyr_extracted_0_01_",i,".RData"))
  df_extract_fin <- df_extract_fin %>% filter(year==2002)
  dfs[[i]] <- df_extract_fin
}

mb <- bind_rows(dfs)

aux <- mb %>% # remember the name for the RData file is "df_extract_fin"
  as.data.frame() %>%
  mutate(totl_area = rowSums(across(-c(cell_id_0_01, year)))) %>%
  mutate(across(-c(cell_id_0_01, year, totl_area), ~ .x / totl_area)) %>% # calculate the share
  rename_with(~ str_replace(., "^(.*)$", "\\1_share"), -c(cell_id_0_01, year, totl_area))

# I need to obtain those cells that are outside BRA, it appears because of geometry boundary differences from different sources
cell_out <- aux %>% 
  filter(nodata_share == 1) %>% 
  pull(cell_id_0_01)

sf <- aux %>% 
  filter(!cell_id_0_01 %in% cell_out) %>% # filter those totally outside BRA, it appears because of geometry boundary differences from different sources
  mutate(nodata_share = ifelse(nodata_share > 0, 0, nodata_share),
         sum_shares = rowSums(dplyr::select(., -c(cell_id_0_01, year, totl_area, nodata_share)))) %>%   # those cells with nodata_share > 0 are all boundary cells
  mutate(across(-c(cell_id_0_01, year, totl_area, nodata_share, sum_shares), ~ . / sum_shares)) %>%  
  mutate(forest_share = fores_form_share + savan_form_share + Mangrove_share + flood_fores_share + wood_sand_veg_share,
         noforest_natur_form_share = wetland_share + grassland_share + hypersa_tid_fl_share + rocky_outcrop_share + herba_sand_veg_share + oth_nonforst_form_share,
         novege_share = bea_du_sand_share + urban_infra_share + mining_share + oth_non_veg_area_share,
         water_share = riv_lak_oce_share + Aquaculture_share,
         pasture_share = pasture_share,
         other_agri_share = temp_crp_share + sugarcane_share + rice_share + cotton_share + oth_temp_crp_share + peren_crp_share + coffee_share + citrus_share + palm_oil_share + oth_peren_crp_share,
         fores_plt_share = fores_plant_share,
         mosaic_share = mosa_uses_share,
         soy_share = soybean_share,
         max_share = pmax(forest_share, noforest_natur_form_share, novege_share, water_share, pasture_share, other_agri_share, fores_plt_share, mosaic_share,soy_share)) %>% 
  mutate(land_use = case_when(
    max_share == forest_share ~ "forest",
    max_share == noforest_natur_form_share ~ "noforest_natur_form",
    max_share == novege_share ~ "novege_area",
    max_share == water_share ~ "water",
    max_share == pasture_share ~ "pasture",
    max_share == other_agri_share ~ "other_agri",
    max_share == fores_plt_share ~ "fores_plt",
    max_share == mosaic_share ~ "mosaic",
    max_share == soy_share ~ "soy",
    TRUE ~ NA_character_
  )) %>% 
  dplyr::select(c(cell_id_0_01, year, totl_area, max_share, land_use))

sf_file <- final_grid_0_01 %>% 
  left_join(sf %>% dplyr::select(c(cell_id_0_01,land_use)), by = "cell_id_0_01") 

agribusiness <- sf_file %>% filter(land_use == "soy") %>% 
  st_intersection(file) %>% st_union()
  
agribusiness_raster <- rasterize(vect(agribusiness), template, field = 1,
                                    fun = "sum", cover = TRUE, background = 0)

agribusiness_df <- as.data.frame(agribusiness_raster, cells = TRUE, na.rm = FALSE)
names(agribusiness_df) <- c("cell", "coverage")
save(agribusiness_df, file = "03_outputs/agribusiness_df.RData")

agribusiness_raster <- ifel(agribusiness_raster >= 0.1, 1, 0)

writeRaster(agribusiness_raster, "03_outputs/exercise_amabiome/agribusiness001.tif", overwrite=TRUE)

# STEP 2 - MapBiomes deforestation ---------------------------------------------

# load, filter to agribusiness
deforestation <- sf_file %>% filter(land_use %in% c("soy", "other_agri", "pasture")) %>% 
  st_intersection(file) %>% st_union()

deforestation_raster <- rasterize(vect(deforestation), template, field = 1,
                                 fun = "sum", cover = TRUE, background = 0)

deforestation_raster <- ifel(deforestation_raster >= 0.5, 1, 0)

writeRaster(deforestation_raster, "03_outputs/exercise_amabiome/deforestation001.tif", overwrite=TRUE)

# STEP 3 - Soy productivity data -----------------------------------------------

# Load
productivity <- rast("02_inputs/High_soyb200b_yld.tif")
# Crop
productivity <- crop(productivity, st_bbox(file))

#productivity <- mask(productivity, template)

# Resample
productivity_raster <- resample(productivity, template, method = "near")

productivity_df <- as.data.frame(productivity_raster, cells = TRUE, na.rm = FALSE)
names(productivity_df) <- c("cell", "soy_productivity")
save(productivity_df, file = "03_outputs/productivity_df.RData")

writeRaster(productivity_raster, "03_outputs/exercise_amabiome/productivity001.tif", overwrite=TRUE)
# STEP 4 - Legal Status data ---------------------------------------------------

# Loading indiginous sites and conservation units and intersecting with MT 
conservation_units <- st_read("02_inputs/conservation_units_legal_amazon.shp")
indiginous_sites <- st_read("02_inputs/indigenous_area_legal_amazon.shp")

# CRS sanity check
if (is.na(st_crs(conservation_units)$epsg) || st_crs(conservation_units)$epsg != 4326) {
  conservation_units <- st_transform(conservation_units, st_crs(file))
}
if (is.na(st_crs(indiginous_sites)$epsg) || st_crs(indiginous_sites)$epsg != 4326) {
  indiginous_sites <- st_transform(indiginous_sites, st_crs(file))
}

# Some polygon here has a duplicated vertex, so we use this function
conservation_units_raster <- rasterize(conservation_units, template, field = 1, 
                                       fun = "sum", cover = TRUE, background = 0)
indiginous_sites_raster <- rasterize(indiginous_sites, template, field = 1, 
                                     fun = "sum", cover = TRUE, background = 0)

conservation_units_raster <- ifel(conservation_units_raster >= 0.5, 1, 0)
indiginous_sites_raster <- ifel(indiginous_sites_raster >= 0.5, 1, 0)


writeRaster(conservation_units_raster, "03_outputs/exercise_amabiome/conservation001.tif", overwrite=TRUE)
writeRaster(indiginous_sites_raster, "03_outputs/exercise_amabiome/indiginous001.tif", overwrite=TRUE)


# SOME SANITY CHECKS -----------------------------------------------------------

# Regression 1: whole amazon

df <- agribusiness_df %>% left_join(productivity_df) %>% 
  mutate(soy_productivity = soy_productivity/1000,
         coverage = 100*coverage)
reg <- feols(coverage ~ soy_productivity, data = df)


# Regression 2: only buffer
buff <- st_read("02_inputs/boundary_3deg_final.gpkg")

template2 <- rast(ext(-63.4675523, 40.4281765 , -19.6397217, -0.8849127), resolution = 0.01, crs = "epsg:4326") %>% 
  setValues(1)

agri_temp <- sf_file %>% filter(land_use == "soy") %>% 
  st_intersection(buff) %>% st_union()

agri_temp <- rasterize(vect(agri_temp), template2, field = 1,
                                 fun = "sum", cover = TRUE, background = 0)

agri_df <- as.data.frame(agri_temp, cells = TRUE, na.rm = FALSE)
names(agri_df) <- c("cell", "coverage")

prod <- resample(productivity, template2, method = "near")

prod_df <- as.data.frame(prod, cells = TRUE, na.rm = FALSE)
names(prod_df) <- c("cell", "soy_productivity")

df <- agri_df %>% left_join(prod_df) %>%    
  mutate(soy_productivity = soy_productivity/1000,
  coverage = 100*coverage)
reg2 <- feols(coverage ~ soy_productivity, data = df)

# Regression 3: subset of non-forest

deforestation_df <- as.data.frame(deforestation_raster, cells = TRUE, na.rm = FALSE)
names(deforestation_df) <- c("cell","non_forest")

df <- agribusiness_df %>% left_join(productivity_df) %>% left_join(deforestation_df) %>% 
  mutate(soy_productivity = soy_productivity/1000,
         coverage = 100*coverage)

reg3 <- feols(coverage ~ soy_productivity, data = df %>%  filter(non_forest==1))
  
etable(list(reg,reg2,reg3),
       tex = TRUE,
       digits = "r4",
       digits.stats = 2,
       #   se.row = T,
       coef.just = "c",
       se.below = T,
       signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.1), 
       depvar = F,
       placement = "h!", 
       #label = auxLabel,
       #postprocess.tex = newAdjustbox_etable, # adjustbox only for table
       adjustbox = TRUE, # argument for newAdjustbox_etabl
       file = "04_results/sanity_check.tex",
       replace = T
       
)
# Productivity map
productivity_masked <- mask(productivity_raster, file)

prod_df <- as.data.frame(productivity_masked, xy = TRUE)

names(prod_df)[3] <- "value"                  # ensure a clean column name

plt <- ggplot(prod_df, aes(x = x, y = y, fill = value/1000)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c() +
  labs(fill = "High Input Soy Productivity") +
  theme_void() +  # removes axes, ticks, grid
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
ggsave(plot = plt, filename = "04_results/prod.png", dpi = 300, height = 6, width = 8)

# Land use map 

sf_file_amz <- sf_file %>% st_intersection(file)

plt <- ggplot() +
  geom_sf(data = sf_file_amz, aes(fill = land_use), color = NA, alpha = 0.6) +
  scale_fill_manual(
    name = "",
    values = c(
      "forest"              = "#1b7837",
      "noforest_natur_form" = "#a6dba0",
      "pasture"             = "#d9f0d3",
      "water"               = "#2c7fb8",
      "novege_area"         = "#f6e8c3",
      "other_agri"          = "#dfc27d",
      "soy"                 = "#fdae61",
      "fores_plt"           = "#80cdc1",
      "mosaic"              = "#9970ab"
    ),
    labels = c("Forest", "Nat-Non-Forest", "Pasture", "Water",
               "Urban/Mining", "Other-Agri", "Soy", "Forest Plt", "Mosaic"), 
    na.value = "grey80"
  ) +
  geom_sf(data = file, fill = NA, color ="black", size = 5) +  # Thicker boundary and legend inclusion
  geom_sf(data = conservation_units, fill = NA, aes(color = "Conservation Units"), size = 5) +  # Thicker boundary and legend inclusion
  geom_sf(data = indiginous_sites, fill = NA, aes(color = "Indiginous sites"), size = 5) +           # Thicker boundary and legend inclusion
  scale_color_manual(values = c("Conservation Units" = "darkgreen", "Indiginous sites" = "darkblue"),
                     name = "",
                     guide = guide_legend(override.aes = list(size = 1.5))) +  # Adjust legend line thickness
  theme_void() +
  theme(legend.position = "bottom") 
ggsave(plot = plt, filename = "04_results/land_use.png")


# Land use map 

sf_file_amz2 <- sf_file_amz %>% filter(land_use %in% c("soy", "other_agri", "pasture"))

plt <- ggplot() +
  geom_sf(data = file, fill = "lightgreen", color ="black", size = 5) +  # Thicker boundary and legend inclusion
  geom_sf(data = conservation_units, fill = NA, aes(color = "Conservation Units"), size = 5) +  # Thicker boundary and legend inclusion
  geom_sf(data = indiginous_sites, fill = NA, aes(color = "Indiginous sites"), size = 5) +           # Thicker boundary and legend inclusion
  scale_color_manual(values = c("Conservation Units" = "darkgreen", "Indiginous sites" = "darkblue"),
                     name = "",
                     guide = guide_legend(override.aes = list(size = 1.5))) +  # Adjust legend line thickness
  geom_sf(data = sf_file_amz2, aes(fill = land_use), color = NA, alpha = 0.6) +
  scale_fill_manual(
    name = "",
    values = c("soy" = "#fdae61", "other_agri" = "red", "pasture" = "purple"), 
    breaks = c("soy", "other_agri", "pasture"),          # <-- match your data values
    labels = c("Soy", "Other Agri" , "Pasture"),
    na.value = "grey80"
  ) +
  theme_void()+
  theme(legend.position = "bottom") 
ggsave(plot = plt, filename = "04_results/land_use_simplified.png", dpi = 300, height = 6, width = 8)
