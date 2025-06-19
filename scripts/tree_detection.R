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

# input path to raw forestry data
forst_path <- file.path(raw_data_dir, 'forst')

# input path to processed forestry data (PSI)
forst_psi_path <- processed_data_dir

# name of the input CHM
ndsm_name <- 'ndsm_np_kellerwald_edersee.tif'

# name of the forestry data
forst_name <- 'hf_fe_geometrien.shp'

# name of the processed forestry data (PSI)
psi_name <- 'psi_kellerwald_crown_diameter.csv'



# 02 - data reading
#-------------------------------------

# read full CHM
ndsm <- terra::rast(file.path(ndsm_path, ndsm_name))
ndsm

# read raw forestry data
fwk <- sf::st_read(file.path(forst_path, forst_name))
fwk

# read processed forestry data (PSI)
psi <- read.csv(file.path(forst_psi_path, psi_name))
psi



# 03 - pre-processing
#-------------------------------------

# select forest subregions (stands)
stand_1 <- fwk %>% 
  dplyr::filter(FO_IDFBBCH == '0718,1,2212,51,A,1,2016,2005')

stand_2 <- fwk %>% 
  dplyr::filter(FO_IDFBBCH == '0718,1,2212,93,,1,2016,2005')

stand_3 <- fwk %>% 
  dplyr::filter(FO_IDFBBCH == '0718,1,2212,308,,1,2016,2005')

# crop CHM to a smaller subregions
ndsm_stand_1 <- terra::crop(ndsm, stand_1, mask = T)
ndsm_stand_2 <- terra::crop(ndsm, stand_2, mask = T)
ndsm_stand_3 <- terra::crop(ndsm, stand_3, mask = T)
ndsm_stand_1
ndsm_stand_2
ndsm_stand_3

# quick overview
par_org <- par()
par(mfrow = c(1,3))
terra::plot(ndsm_stand_1, col = lidR::height.colors(50))
terra::plot(ndsm_stand_2, col = lidR::height.colors(50))
terra::plot(ndsm_stand_3, col = lidR::height.colors(50))
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

# fit models to crown width-height relationship
lm_fit <- stats::lm(KB ~ hoe_mod, data = psi)
loess_fit <- stats::loess(KB ~ hoe_mod, data = psi)

# print model summaries
print(summary(lm_fit))
print(summary(loess_fit))

# set minimum tree height
minHgt <- 2

# create window functions
current_fun <- function(x) ifelse(x < minHgt, NA, (x * 0.1 + 1.2) * 2)  # multiply by 2 to show diameter
lm_fun <- function(x) ifelse(x < minHgt, NA, stats::predict(lm_fit, newdata = data.frame(hoe_mod = x)))
loess_fun <- function(x) ifelse(x < minHgt, NA, stats::predict(loess_fit, newdata = data.frame(hoe_mod = x)))

# create data for plotting
heights_seq <- seq(0, 40, 0.5)
comparison_df <- data.frame(
  height = heights_seq,
  current = current_fun(heights_seq),
  linear_fit = lm_fun(heights_seq),
  loess_fit = loess_fun(heights_seq)
)

# plot showing the calculation of WS
# depending on the chosen ws function
ggplot(psi, aes(x = hoe_mod, y = KB)) +
  geom_point(color = 'black', alpha = 0.5) +
  geom_line(data = comparison_df, aes(x = height, y = current, color = "current linear"), size = 2) +
  geom_line(data = comparison_df, aes(x = height, y = linear_fit, color = "fitted linear"), size = 2) +
  geom_line(data = comparison_df, aes(x = height, y = loess_fit, color = "fitted LOESS"), size = 2) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(name = "models:",
                     values = c("current linear" = "#7A76C2",     
                                "fitted linear" = "#E69F00",         
                                "fitted LOESS" = "#56B4E9")) +       
  labs(
    x = 'tree height [m]',
    y = 'crown width [m]'
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  annotate("text", x = 2.8, y = 11.2, label = "2 m height threshold", hjust = 0, color = "black") +
  scale_y_continuous(breaks = 1:ceiling(max(c(psi$KB, comparison_df$current, 
                                              comparison_df$linear_fit, comparison_df$loess_fit), 
                                            na.rm = T)))

# define ws calculation functions
winFunctions <- list(
  current = function(x) ifelse(x < minHgt, NA, x * 0.1 + 1.2),
  lm = function(x) ifelse(x < minHgt, NA, (coef(lm_fit)[1] + coef(lm_fit)[2] * x) / 2),
  loess = function(x) ifelse(x < minHgt, NA, stats::predict(loess_fit, newdata = data.frame(hoe_mod = x)) / 2)
)

# list of subregions (CHMs of the 3 stands)
ndsm_stands <- list(
  stand_1 = ndsm_stand_1,
  stand_2 = ndsm_stand_2,
  stand_3 = ndsm_stand_3
)

# loop through subregions and window functions
ttops_stands <- list()

for (stand in names(ndsm_stands)) {
  
  for (win_fun in names(winFunctions)) {
    
    # ITD with ForestTools
    ttops <- ForestTools::vwf(
      ndsm_stands[[stand]], winFunctions[[win_fun]], minHgt
    )
    
    # store detected trees in the list
    ttops_stands[[paste0('ttops_', win_fun, '_', stand)]] <- ttops
    
    # write ttops to disk
    file_name <- paste0('ttops_', win_fun, '_', stand, '.gpkg')
    sf::st_write(ttops, file.path(output_dir, file_name), delete_layer = T)
  }
}



# 04 - comparison - statistical analysis
#----------------------------------------

# initialize empty list to collect data
df_list <- list()

# loop through all ttops calculated
# with different ws in all stands
for (name in names(ttops_stands)) {
  
  # extract ttops
  ttops <- ttops_stands[[name]]
  
  # ensure object is not empty
  if (!is.null(ttops) && nrow(ttops) > 0) {
    
    # extract tree IDs and heights
    df_list[[name]] <- data.frame(
      treeID = ttops$treeID,
      height = ttops$height,
      method = sub("ttops_([^_]+)_.*", "\\1", name),
      subregion = sub(".*_(stand_[0-9])", "\\1", name)
    )
  }
}

# combine all extracted trees into one dataframe
df_trees <- bind_rows(df_list)

# summarize number of trees per method and subregion
tree_counts <- df_trees %>%
  group_by(subregion, method) %>%
  summarise(count = n())

writexl::write_xlsx(tree_counts, path = file.path(output_dir, 'n_trees_detected_method_comp.xlsx'))

# define colors to match the earlier plot
plot_colors <- c("#7A76C2", "#E69F00", "#56B4E9") 
names(plot_colors) <- c("current", "lm", "loess")

# boxplot of tree heights per method and subregion
ggplot(df_trees, aes(x = subregion, y = height, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_jitter(
    aes(color = method),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.2, alpha = 0.5, show.legend = F
  ) +
  labs(
    title = 'Tree Height Distribution by Method',
    x = '',
    y = 'Tree Height (m)'
  ) +
  scale_fill_manual(
    values = plot_colors,
    labels = c("Current", "Linear Model", "LOESS")
  ) +
  scale_color_manual(
    values = plot_colors
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )

# cumulative distribution of detected tree heights
ggplot(df_trees, aes(x = height, color = method)) +
  stat_ecdf(geom = "step", size = 1) +
  facet_wrap(~subregion) +
  labs(title = 'Cumulative Distribution of Tree Heights by Method',
       x = 'Height (m)', 
       y = 'Cumulative Proportion') +
  scale_color_manual(values = plot_colors,
                     name = "Method",
                     labels = c("Current", "Linear Model", "LOESS")) +
  theme_minimal()

library(ggrepel)
library(dplyr)

# density plot with lines and direct labels (all subregions together)
dens_data <- df_trees %>%
  group_by(method) %>%
  do({
    dens <- density(.$height)
    data.frame(height = dens$x, density = dens$y)
  })

label_data <- dens_data %>%
  group_by(method) %>%
  filter(height == max(height))

ggplot(dens_data, aes(x = height, y = density, color = method)) +
  geom_line(size = 1.2) +
  geom_text_repel(
    data = label_data,
    aes(label = method),
    nudge_x = 1,
    direction = "y",
    hjust = 0,
    segment.color = NA,
    fontface = "bold",
    show.legend = FALSE
  ) +
  labs(
    title = 'Density of Tree Heights by Method (with Labels)',
    x = 'Tree Height (m)',
    y = 'Density'
  ) +
  scale_color_manual(
    values = plot_colors,
    labels = c("Current", "Linear Model", "LOESS")
  ) +
  theme_minimal() +
  guides(color = 'none')

# density plot with lines and direct labels (faceted by subregion)
dens_data_facet <- df_trees %>%
  group_by(subregion, method) %>%
  do({
    dens <- density(.$height)
    data.frame(height = dens$x, density = dens$y)
  })

label_data_facet <- dens_data_facet %>%
  group_by(subregion, method) %>%
  filter(height == max(height))

ggplot(dens_data_facet, aes(x = height, y = density, color = method)) +
  geom_line(size = 1.2) +
  geom_text_repel(
    data = label_data_facet,
    aes(label = method),
    nudge_x = 1,
    direction = "y",
    hjust = 0,
    segment.color = NA,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~subregion) +
  labs(
    title = 'Density of Tree Heights by Method and Subregion (with Labels)',
    x = 'Tree Height (m)',
    y = 'Density'
  ) +
  scale_color_manual(
    values = plot_colors,
    labels = c("Current", "Linear Model", "LOESS")
  ) +
  theme_minimal() +
  guides(color = 'none')



# 05 - application to whole area
#----------------------------------------

# plot whole area with forest stands
terra::plot(ndsm, col = lidR::height.colors(50))
terra::plot(fwk$geom, alpha = 1, border = 'black', add = T)

# mask CHM to forest stands
ndsm_aoi <- terra::crop(ndsm, fwk, mask = T)
terra::plot(ndsm_aoi, col = lidR::height.colors(50))

# define ws calculation
winFunction <- function(x){x * 0.1 + 1.2}

# set minimum tree height
minHgt <- 2

# ITD with ForestTools
ttops_ndsm_aoi <- ForestTools::vwf(ndsm_aoi, winFunction, minHgt)
head(ttops_ndsm_aoi)
sf::st_write(ttops_ndsm_aoi, file.path(output_dir, 'ttops_np_kellerwald_edersee.gpkg'))

### new detection ###

# define linear fit for crown width-height relationship
lm_fit <- stats::lm(KB ~ hoe_mod, data = psi)

# set minimum tree height
minHgt <- 2

# define window size function using the linear fit
winFunction <- function(x) ifelse(x < minHgt, NA, (coef(lm_fit)[1] + coef(lm_fit)[2] * x) / 2)

# ITD with ForestTools using linear fit
ttops_ndsm_aoi_lm <- ForestTools::vwf(ndsm_aoi, winFunction, minHgt)
head(ttops_ndsm_aoi_lm)
sf::st_write(ttops_ndsm_aoi_lm, file.path(output_dir, 'ttops_np_kellerwald_edersee_lm.gpkg'))



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
sf::st_transform(ttops_ndsm_aoi_lm, sf::st_crs(fwk))
ttops_ndsm_aoi_ts <- sf::st_join(ttops_ndsm_aoi_lm, 
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
             file.path(output_dir, 'ttops_ts_stnr_area_np_kellerwald_edersee_lm.gpkg'))








