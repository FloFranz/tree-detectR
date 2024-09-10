#------------------------------------------------------------------------
# Name:         tree_detection.R
# Description:  In this script, different individual tree detection (ITD)
#               methods are compared to a subset of a larger forest area.
#               The lidR package is used for ITD.
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
forst_name <- 'fwk_np_kellerwald_edersee.gpkg'



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
#-------------------------------------

# local maximum filter (LMF) with
# variable window size (WS) depending on tree heights (TH)

# function implementing a linear relationship
# between TH and WS
calc_ws_linear = function(x) { 
  
  # define end points for the linear relationship
  # TH < 2m --> WS = 3
  # TH > 20m --> WS = 5
  x1 <- 2; y1 <- 3
  x2 <- 20; y2 <- 5
  
  # calculate slope (m)
  m <- (y2 - y1) / (x2 - x1)
  
  # calculate intercept (b)
  b <- y1 - m * x1
  
  # calculate WS using linear relationship
  y <- m * x + b
  y[x < 2] <- 3
  y[x > 20] <- 5
  
  return(y)
}

# function implementing a non-linar relationship
# (exponential decay) between TH and WS
# see: https://r-lidar.github.io/lidRbook/itd-its.html
calc_ws_exp_decay <- function(x) {
  
  # calculate exponential decay
  y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3
  
  # define end points for the non-linear relationship
  # TH < 2m --> WS = 3
  # TH > 20m --> WS = 5
  y[x < 2] <- 3
  y[x > 20] <- 5
  
  return(y)
}

# function implementing a non-linear relationship
# (spline interpolation) between TH and WS
calc_ws_spline_int <- function(x) {
  
  # control points for the spline
  control_heights <- c(0, 2, 20, 25)
  control_ws <- c(3, 3, 5, 5)
  
  # create spline function
  ws_func <- stats::splinefun(control_heights, control_ws, method = 'natural')
  
  # calculate WS
  y <- ws_func(x)
  
  # ensure end points
  y[x < 2] <- 3
  y[x > 20] <- 5
  
  return(y)
}

# plot showing the calculation of WS
# depending on the choosen function
heights <- seq(-5,30,0.5)
ws_linear <- calc_ws_linear(heights)
ws_exp_decay <- calc_ws_exp_decay(heights)
ws_spline_int <- calc_ws_spline_int(heights)

plot(heights, ws_linear, type = 'l', col = 'red', ylim = c(0,5), 
     xlab = 'height (m)', ylab = 'window size')
lines(heights, ws_exp_decay, col = 'blue')
lines(heights, ws_spline_int, col = 'green')
legend('bottomright', legend = c('linear', 'exponential decay', 'spline interpolation'), 
       col = c('blue', 'red', 'green'), lty = 1)

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

# quick view
terra::plot(ndsm_subregion,
            col = lidR::height.colors(50))
terra::plot(sf::st_geometry(ttops), add = T, pch = 3)

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







