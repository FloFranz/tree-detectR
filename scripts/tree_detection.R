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

# select forest subregions (stands)
subregion_1 <- fwk %>% 
  dplyr::filter(FO_IDFBBCH == '0718,1,2212,6504,B,2,2013,2005')
subregion_1

subregion_2 <- fwk %>% 
  dplyr::filter(FO_IDFBBCH == '0718,1,2212,519,B,1,2016,2005')
subregion_2

# crop CHM to a smaller subregions
ndsm_subregion_1 <- terra::crop(ndsm, subregion_1, mask = T)
ndsm_subregion_2 <- terra::crop(ndsm, subregion_2, mask = T)
ndsm_subregion_1
ndsm_subregion_2

# quick overview
par_org <- par()
par(mfrow = c(1,2))
terra::plot(ndsm_subregion_1, col = lidR::height.colors(50))
terra::plot(ndsm_subregion_2, col = lidR::height.colors(50))
par(par_org)

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


#######################################################
### comparison of different window size calculation ###
#######################################################

# set minimum tree height
minHgt <- 2

calc_ws <- function(x, minHgt){
  ws_1 <- ifelse(x < minHgt, NA, (x * 0.06 + 0.5) + 2.38)
  ws_2 <- ifelse(x < minHgt, NA, (x * 0.04 + 0.75) + 2.38)
  ws_3 <- ifelse(x < minHgt, NA, (x * 0.03 + 0.95) + 2.38)
  ws_4 <- ifelse(x < minHgt, NA, (x * 0.02 + 1.2) + 2.38)
  return(data.frame(ws_1 = ws_1, ws_2 = ws_2, ws_3 = ws_3, ws_4 = ws_4))
}

# plot showing the calculation of WS
# depending on the choosen ws function
heights <- seq(-5,40,0.5)
ws_values <- calc_ws(heights, minHgt)

colors <- RColorBrewer::brewer.pal(4, 'Set1')

plot(heights, ws_values$ws_1, type = 'l', col = colors[1], lwd = 2,
     ylim = c(0, 10), xlab = 'Height (m)', ylab = 'Window Size (m)')
lines(heights, ws_values$ws_2, col = colors[2], lwd = 2)
lines(heights, ws_values$ws_3, col = colors[3], lwd = 2)
lines(heights, ws_values$ws_4, col = colors[4], lwd = 2)
legend('bottomright', 
       legend = c('original', 'flat 1', 'flat 2', 'flat 3'), 
       col = colors, lty = 1, lwd = 2)

# define ws calculation
winFunction_1 <- function(x){x * 0.06 + 0.5} + 2.38
winFunction_2 <- function(x){x * 0.04 + 0.75} + 2.38
winFunction_3 <- function(x){x * 0.03 + 0.95} + 2.38
winFunction_4 <- function(x){x * 0.02 + 1.2} + 2.38

# ITD with ForestTools
ttops1_subregion_1 <- ForestTools::vwf(ndsm_subregion_1, winFunction_1, minHgt)
ttops2_subregion_1 <- ForestTools::vwf(ndsm_subregion_1, winFunction_2, minHgt)
ttops3_subregion_1 <- ForestTools::vwf(ndsm_subregion_1, winFunction_3, minHgt)
ttops4_subregion_1 <- ForestTools::vwf(ndsm_subregion_1, winFunction_4, minHgt)

ttops1_subregion_2 <- ForestTools::vwf(ndsm_subregion_2, winFunction_1, minHgt)
ttops2_subregion_2 <- ForestTools::vwf(ndsm_subregion_2, winFunction_2, minHgt)
ttops3_subregion_2 <- ForestTools::vwf(ndsm_subregion_2, winFunction_3, minHgt)
ttops4_subregion_2 <- ForestTools::vwf(ndsm_subregion_2, winFunction_4, minHgt)

# new trees detected in flat 1 but not in original
new_trees_2_vs_1 <- st_difference(ttops2_subregion_1, st_union(ttops1_subregion_1))

# missing trees detected in original but not in flat 1
missing_trees_2_vs_1 <- st_difference(ttops1_subregion_1, st_union(ttops2_subregion_1))

# count number of new and missing trees
num_new_trees_2_vs_1 <- nrow(new_trees_2_vs_1)
num_missing_trees_2_vs_1 <- nrow(missing_trees_2_vs_1)

# calculate net difference
net_difference_2_vs_1 <- num_new_trees_2_vs_1 - num_missing_trees_2_vs_1

cat('Flat 1 vs Original:\n')
cat('New Trees Detected:', num_new_trees_2_vs_1, '\n')
cat('Missing Trees:', num_missing_trees_2_vs_1, '\n')
cat('Net Difference:', net_difference_2_vs_1, '\n')

# write to disk
sf::st_write(ttops1_subregion_1, file.path(output_dir, 'ttops_subregion_1_org.gpkg'))
sf::st_write(ttops2_subregion_1, file.path(output_dir, 'ttops_subregion_1_flat1.gpkg'))
sf::st_write(ttops3_subregion_1, file.path(output_dir, 'ttops_subregion_1_flat2.gpkg'))
sf::st_write(ttops4_subregion_1, file.path(output_dir, 'ttops_subregion_1_flat3.gpkg'))



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
winFunction_new <- function(x){x * 0.04 + 0.75} + 2.38

# set minimum tree height
minHgt <- 2

# ITD with ForestTools
ttops_ndsm_aoi <- ForestTools::vwf(ndsm_aoi, winFunction_new, minHgt)
head(ttops_ndsm_aoi)
sf::st_write(ttops_ndsm_aoi, file.path(output_dir, 'ttops_np_kellerwald_edersee_new.gpkg'))



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
# if the following message is printed:
# "Fehler in st_geos_binop("intersects", x, y, sparse = sparse, prepared = prepared,  : 
# st_crs(x) == st_crs(y) ist nicht TRUE
# add this command:
# sf::st_transform(ttops_ndsm_aoi, sf::st_crs(fwk))
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

# assign three forest classification levels (Abteilung, Unterabteilung, Bestand),
# and all three in one key to each detected tree top
ttops_ndsm_aoi_ts_stands <- sf::st_join(
  ttops_ndsm_aoi_ts, 
  fwk[, c('FO_IDFBBCH', 'FO_HRW4ABT', 'FO_HRW3_BE', 'FO_HRW2_BE')])

head(ttops_ndsm_aoi_ts_stands)

# rename columns
ttops_ndsm_aoi_ts_stands <- ttops_ndsm_aoi_ts_stands %>%
  dplyr::rename(
    Schluessel = FO_IDFBBCH,
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
             file.path(output_dir, 'ttops_ts_stnr_area_np_kellerwald_edersee_new.gpkg'))








