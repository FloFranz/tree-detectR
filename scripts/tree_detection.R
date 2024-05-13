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
calc_ws = function(x) { 
  
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
calc_ws <- function(x) {
  
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
calc_ws <- function(x) {
  
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

# example plot showing the calculation of WS
# depending on the choosen function
heights <- seq(-5,30,0.5)
ws <- calc_ws(heights)
plot(heights, ws, type = 'l',  ylim = c(0,5))


# detect tree tops using LMF with 
# the functions defined before
ttops <- lidR::locate_trees(ndsm_subregion,
                            lidR::lmf(ws = calc_ws,
                                      shape = 'circular'))

# write to disk
sf::st_write(ttops, file.path(output_dir, 'ttops_nonlinear_spline.gpkg'))

# quick view
terra::plot(ndsm_subregion,
            col = lidR::height.colors(50))
terra::plot(sf::st_geometry(ttops), add = T, pch = 3)







