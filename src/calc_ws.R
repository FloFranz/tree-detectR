#------------------------------------------------------------------------------
# Name:         calc_ws.R
# Description:  In this script, different methods to calculate window size (WS)
#               depending on tree height (TH) are defined. They can be used for
#               individual tree detection (ITD) with variable WS.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#------------------------------------------------------------------------------

# function implementing a linear relationship
# between TH and WS
calc_ws_linear = function(x) { 
  
  # define end points for the linear relationship
  # TH < 2m --> WS = 3
  # TH > 20m --> WS = 5
  x1 <- 2; y1 <- 3
  x2 <- 30; y2 <- 6
  
  # calculate slope (m)
  m <- (y2 - y1) / (x2 - x1)
  
  # calculate intercept (b)
  b <- y1 - m * x1
  
  # calculate WS using linear relationship
  y <- m * x + b
  y[x < 2] <- 3
  y[x > 30] <- 6
  
  return(y)
}

# function implementing a non-linar relationship
# (exponential decay) between TH and WS
# see: https://r-lidar.github.io/lidRbook/itd-its.html
calc_ws_exp_decay <- function(x) {
  
  # calculate exponential decay
  y <- 3.47 * (-(exp(-0.07*(x-2)) - 1)) + 3
  
  # define end points for the non-linear relationship
  # TH < 2m --> WS = 3
  # TH > 20m --> WS = 5
  y[x < 2] <- 3
  y[x > 30] <- 6
  
  return(y)
}

# function implementing a non-linear relationship
# (spline interpolation) between TH and WS
calc_ws_spline_int <- function(x) {
  
  # control points for the spline
  control_heights <- c(0, 2, 30, 35)
  control_ws <- c(3, 3, 6, 6)
  
  # create spline function
  ws_func <- stats::splinefun(control_heights, control_ws, method = 'natural')
  
  # calculate WS
  y <- ws_func(x)
  
  # ensure end points
  y[x < 2] <- 3
  y[x > 30] <- 6
  
  return(y)
}