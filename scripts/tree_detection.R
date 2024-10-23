#------------------------------------------------------------------------
# Name:         tree_detection.R
# Description:  In this script, different individual tree detection (ITD)
#               methods are compared to a subset of a larger forest area.
#               The packages lidR and ForestTools are used for ITD.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)




# 01 - set file paths
#-------------------------------------

# input path to normalized canopy height model (CHM)
ndsm_path <- file.path(raw_data_dir, 'nDSM_mosaic')

# input path to forestry data
forst_path <- file.path(raw_data_dir, 'forst')

# name of the input CHM
ndsm_name <- 'ndsm_np_kellerwald_edersee.tif'

# name of the forestry data
forst_name <- 'hf_fe_geometrien.shp'



# 02 - data reading
#-------------------------------------

# read full CHM
ndsm <- terra::rast(file.path(ndsm_path, ndsm_name))
ndsm

# read forestry data
fwk <- sf::st_read(file.path(forst_path, forst_name))
fwk



# 03 - pre-processing
#-------------------------------------

# select forest subregion based on the 'Abteilung'
subregion <- fwk %>% 
  dplyr::filter(FO_HRW4ABT_BEZ == '84' | FO_HRW4ABT_BEZ == '83')
subregion

# crop CHM to a smaller subregion
ndsm_subregion <- terra::crop(ndsm, subregion, mask = T)
ndsm_subregion

# quick overview
terra::plot(ndsm_subregion,
            col = lidR::height.colors(50))



# 04 - ITD
# local maximum filter (LMF) with
# variable window size (WS) depending
# on tree heights (TH)
#-------------------------------------

# source functions for calculating ws
source('src/calc_ws.R', local = TRUE)

# using package ForestTools
# https://github.com/andrew-plowright/ForestTools
minHgt <- 2
calc_ws_ForestTools <- function(x){
  ws <- ifelse(x < minHgt, NA, (x * 0.06 + 0.5) + 2.38)
  return(ws)
}

# plot showing the calculation of WS
# depending on the choosen function
heights <- seq(-5,35,0.5)
ws_linear <- calc_ws_linear(heights)
ws_exp_decay <- calc_ws_exp_decay(heights)
ws_spline_int <- calc_ws_spline_int(heights)
ws_ForestTools <- calc_ws_ForestTools(heights)

colors <- RColorBrewer::brewer.pal(4, 'Set1')

plot(heights, ws_linear, type = 'l', col = colors[1], lwd = 2,
     ylim = c(0,6), xlab = 'height (m)', ylab = 'window size')
lines(heights, ws_exp_decay, col = colors[2], lwd = 2)
lines(heights, ws_spline_int, col = colors[3], lwd = 2)
lines(heights, ws_ForestTools, col = colors[4], lwd = 2)
legend('bottomright', legend = c('linear', 'exponential decay', 'spline interpolation', 'ForestTools'), 
       col = c(colors[1], colors[2], colors[3], colors[4]), lty = 1, lwd = 2)

# detect tree tops using LMF with 
# the functions defined before
ws_methods <- list(
  linear = calc_ws_linear,
  exp_decay = calc_ws_exp_decay,
  spline_int = calc_ws_spline_int
)

ttops <- list()

for (method in names(ws_methods)) {
  
  calc_ws_func <- ws_methods[[method]]
  
  ttops[[method]] <- lidR::locate_trees(
    ndsm_subregion,
    lidR::lmf(ws = function(x) calc_ws_func(x), shape = 'circular')
  )
  
  output_file_name <- file.path(output_dir, paste0('ttops_', method, '.gpkg'))
  
  sf::st_write(ttops[[method]], output_file_name)
  
  cat('tree top detection using method', method, 'completed.\n')
  
}

ttops_forest_tools <- ForestTools::vwf(ndsm_subregion, calc_ws_ForestTools, minHgt)
sf::st_write(ttops_forest_tools, file.path(output_dir, 'ttops_ForestTools.gpkg'))

# quick view
ttops_2d <- lapply(ttops, function(t) {
  # remove Z (height) coordinate
  t_2d <- sf::st_zm(t)
  return(t_2d)
})

mapview::mapview(ndsm_subregion, col = colorRampPalette(c("black", "white"))(50),
                 map.types = c('OpenStreetMap.DE', 'Esri.WorldImagery'),
                 na.color = NA) +
  
  mapview::mapview(ttops_2d$exp_decay, cex = 2 ,label = 'Z', col.regions = 'blue', color = 'blue', layer.name = 'ttops_exp_decay') + 
  mapview::mapview(ttops_2d$linear, cex = 2 ,label = 'Z', col.regions = 'red', color = 'red', layer.name = 'ttops_linear') +
  mapview::mapview(ttops_2d$spline_int, cex = 2 ,label = 'Z', col.regions = 'green', color = 'green', layer.name = 'ttops_spline_int')



# 04 - comparison - statistical analysis
#----------------------------------------

# Create a data frame for plotting
df <- data.frame(
  height = c(ttops$linear$Z, ttops$exp_decay$Z, ttops$spline_int$Z),
  method = c(rep('linear', length(ttops$linear$Z)), rep('exponential decay', length(ttops$exp_decay$Z)), rep('spline interpolation', length(ttops$spline_int$Z)))
)

# density plot of detected heights
ggplot(df, aes(x = height, fill = method)) +
  geom_density(color = NA) +
  labs(title = 'height distribution of tree tops',
       x = 'height (Z values)', 
       y = 'density') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'Set1')) +
  theme_minimal()



# 05 - application to whole area
#----------------------------------------

# plot whole area with forest stands
terra::plot(ndsm, col = lidR::height.colors(50))
terra::plot(fwk$geom, alpha = 1, border = 'black', add = T)

# mask CHM to forest stands
ndsm_aoi <- terra::crop(ndsm, fwk, mask = T)
terra::plot(ndsm_aoi, col = lidR::height.colors(50))

# define ws calculation
winFunction <- function(x){x * 0.06 + 0.5} + 2.38

# set minimum tree height
minHgt <- 2

# ITD with ForestTools
ttops_ndsm_aoi <- ForestTools::vwf(ndsm_aoi, winFunction, minHgt)
head(ttops_ndsm_aoi)
sf::st_write(ttops_ndsm_aoi, file.path(output_dir, 'ttops_np_kellerwald_edersee.gpkg'))



# 06 - add tree species information,
# stand number ("Abteilungsnummer"),
# and stand area size
#--------------------------------------

### tree species ###

# clean FB column (tree species) to keep only the species number
# (part before the comma) and exclude the age
fwk$FB_clean <- ifelse(is.na(fwk$FB), NA, substr(fwk$FB, 1, 1))

# convert FB_clean column to numeric
fwk$FB_clean <- as.numeric(fwk$FB_clean)

# define corresponding names of the tree species
tree_species_names <- c("Eiche", "Buche", "Edellaubbäume", "Weichlaubbäume",
                        "Fichte, Tanne", "Douglasie", "Kiefer", "Lärche")

# add new column with species names based on the FB_clean column
fwk$tree_species_name <- tree_species_names[fwk$FB_clean]

# assign tree species to each detected tree top
ttops_ndsm_aoi_ts <- sf::st_join(ttops_ndsm_aoi, 
                                 fwk[, c('FB_clean', 'tree_species_name')])
head(ttops_ndsm_aoi_ts)

# rename columns
ttops_ndsm_aoi_ts <- ttops_ndsm_aoi_ts %>%
  dplyr::rename(
    Baumartcode = FB_clean,
    Baumartname = tree_species_name
  )

### stand number ###

# convert columns FO_HRW4ABT and FO_HRW2_BE to numeric
fwk$FO_HRW4ABT <- as.numeric(fwk$FO_HRW4ABT)
fwk$FO_HRW2_BE <- as.numeric(fwk$FO_HRW2_BE)

# assign three forest classification levels (Abteilung, Unterabteilung, Bestand)
# to each detected tree top
ttops_ndsm_aoi_ts_stands <- sf::st_join(ttops_ndsm_aoi_ts, 
                                        fwk[, c('FO_HRW4ABT', 'FO_HRW3_BE', 'FO_HRW2_BE')])
head(ttops_ndsm_aoi_ts_stands)

# rename columns
ttops_ndsm_aoi_ts_stands <- ttops_ndsm_aoi_ts_stands %>%
  dplyr::rename(
    Abteilung = FO_HRW4ABT,
    Unterabteilung = FO_HRW3_BE ,
    Bestand = FO_HRW2_BE
    )

### stand area size ###

# calculate area of each polygon (Bestand) in fwk
fwk$area <- sf::st_area(fwk)

# convert to numeric
fwk$area <- as.numeric(fwk$area)

# convert from sqm to ha
fwk$area <- fwk$area / 10000

# assign stand area sizes to detected tree tops
ttops_ndsm_aoi_ts_stands_area <- sf::st_join(ttops_ndsm_aoi_ts_stands, 
                                             fwk[, 'area'])

head(ttops_ndsm_aoi_ts_stands_area)

# rename columns
ttops_ndsm_aoi_ts_stands_area <- ttops_ndsm_aoi_ts_stands_area %>%
  dplyr::rename(
    Bestandsflaeche = area
  )

# write to disk
sf::st_write(ttops_ndsm_aoi_ts_stands_area, 
             file.path(output_dir, 'ttops_ts_stnr_area_np_kellerwald_edersee.gpkg'))








