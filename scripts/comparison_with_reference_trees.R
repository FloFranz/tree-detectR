#------------------------------------------------------------------------
# Name:         comparison_ref_trees.R
# Description:  In this script, tree heights obtained with different
#               individual tree detection (ITD) methods from an ALS-based 
#               CHM are compared with terrestrial tree heights.
#               This is done for three different parcels.
#               The packages lidR and ForestTools are used for ITD.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - set file paths
#-------------------------------------

# input path to normalized canopy height model (CHM)
chm_path <- file.path(raw_data_dir, 'nDSM_mosaic')

# input path to forestry data
forst_path <- file.path(raw_data_dir, 'forst')

# name of the input CHM
chm_name <- 'chm_solling_2024_als.tif'

# name of the reference tree data
ref_trees_name <- 'Stammv.shp'



# 02 - data reading
#-------------------------------------

# read full CHM
chm <- terra::rast(file.path(chm_path, chm_name))
chm

# read reference tree data
reference_trees <- sf::st_read(file.path(forst_path, ref_trees_name))
head(reference_trees)



# 03 - pre-processing
#-------------------------------------

# transform the reference trees into the desired CRS
# (ETRS89 / UTM zone 32N (EPSG:25832))
# and remove z values from the geometry
reference_trees <- sf::st_transform(reference_trees, 'EPSG:25832')
reference_trees <- sf::st_zm(reference_trees, drop = T)
head(reference_trees)

# extract the coordinates from the three different parcels
colnames(reference_trees)
parz_x <- reference_trees$x_UTM[reference_trees$art != 611]
parz_x <- as.numeric(substr(parz_x, 3, 20))
parz_y <- reference_trees$y_UTM[reference_trees$art != 611]

# create matrices of the coordinates for each parcel
parcel_1_coords <- matrix(c(parz_x[1:4], parz_x[1], 
                            parz_y[1:4], parz_y[1]),
                          ncol = 2)

parcel_2_coords <- matrix(c(parz_x[5:9], parz_x[5],
                            parz_y[5:9], parz_y[5]),
                          ncol = 2)

parcel_3_coords <- matrix(c(parz_x[12:16], parz_x[12],
                            parz_y[12:16], parz_y[12]),
                          ncol = 2)

# create individual sf polygons for each parcel
parcel_1 <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(parcel_1_coords))), crs = 'EPSG:25832')
parcel_2 <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(parcel_2_coords))), crs = 'EPSG:25832')
parcel_3 <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(parcel_3_coords))), crs = 'EPSG:25832')

# filter the reference trees and adjust column naming
reference_trees <- reference_trees %>%
  dplyr::filter(art == 611) %>%
  dplyr::mutate(x_UTM = as.numeric(substr(x_UTM, 3, 20))) %>%
  dplyr::rename(baumhoehe = 4) %>%
  dplyr::mutate(baumhoehe = as.numeric(baumhoehe) / 10)

# crop CHM to the extent of the reference trees
# and project remove height information from the CRS
reference_trees_extent <- sf::st_bbox(reference_trees)
chm_exp_site <- terra::crop(chm, reference_trees_extent)
chm_exp_site <- terra::project(chm_exp_site, 'EPSG:25832')

# reduce resolution to 1 m
chm_exp_site <- terra::aggregate(chm_exp_site, fact = 2,fun = 'mean')
chm_exp_site

# filter reference trees within each parcel
reference_trees_parcel_1 <- reference_trees[sf::st_intersects(reference_trees, parcel_1, sparse = FALSE), ]
reference_trees_parcel_2 <- reference_trees[sf::st_intersects(reference_trees, parcel_2, sparse = FALSE), ]
reference_trees_parcel_3 <- reference_trees[sf::st_intersects(reference_trees, parcel_3, sparse = FALSE), ]

# plot CHM with parcels and corresponding reference trees
terra::plot(chm_exp_site)
terra::plot(parcel_1$geometry, border = 'white', add = T)
terra::plot(parcel_2$geometry, border = 'white', add = T)
terra::plot(parcel_3$geometry, border = 'white', add = T)
terra::plot(reference_trees_parcel_1$geometry,  pch = 4, add = T)
terra::plot(reference_trees_parcel_2$geometry,  pch = 4, add = T)
terra::plot(reference_trees_parcel_3$geometry,  pch = 4, add = T)

# crop CHM to each parcel
chm_parcel_1 <- terra::crop(chm_exp_site, parcel_1, mask = T)
chm_parcel_2 <- terra::crop(chm_exp_site, parcel_2, mask = T)
chm_parcel_3 <- terra::crop(chm_exp_site, parcel_3, mask = T)

# plot individual CHM parcels with reference trees
terra::plot(chm_parcel_1)
terra::plot(reference_trees_parcel_1$geometry, pch = 4, add = T)
terra::plot(chm_parcel_2)
terra::plot(reference_trees_parcel_2$geometry, pch = 4, add = T)
terra::plot(chm_parcel_3)
terra::plot(reference_trees_parcel_3$geometry, pch = 4, add = T)



# 04 - ITD using different methods
#-------------------------------------

# source functions for calculating ws
source('src/calc_ws.R', local = TRUE)

# define the different tree detection methods
ws_methods <- list(
  linear = calc_ws_linear,
  exp_decay = calc_ws_exp_decay,
  spline_int = calc_ws_spline_int,
  forest_tools = function(chm, minHgt) {
    winFunction <- function(x){x * 0.06 + 0.5} + 2.38
    ForestTools::vwf(chm, winFunction, minHgt)
  }
)

# define the datasets for the different parcels
parcels <- list(
  list(reference = reference_trees_parcel_1, chm = chm_parcel_1),
  list(reference = reference_trees_parcel_2, chm = chm_parcel_2),
  list(reference = reference_trees_parcel_3, chm = chm_parcel_3)
)

# initialize empty list to store results and plots for each method and parcel
results_all_methods <- list()
plots_all_methods <- list()

# loop over each method (linear, exp_decay, spline_int, forest_tools)
for (method in names(ws_methods)) {
  
  ###########
  ### ITD ###
  ###########
  
  cat('\nProcessing method:', method, '\n')
  
  # initialize empty list to store results and plots for each parcel
  method_results <- list()
  method_plots <- list()
  
  # get function for the current method
  calc_ws_func <- ws_methods[[method]]
  
  # loop over each parcel
  for (parcel_idx in 1:length(parcels)) {
    
    # extract reference trees and CHM for the current parcel
    reference_trees <- parcels[[parcel_idx]]$reference
    chm <- parcels[[parcel_idx]]$chm
    
    # ensure reference trees are not empty
    if (nrow(reference_trees) == 0) {
      cat('No reference trees found for parcel', parcel_idx, 'with method', method, '\n')
      next
    }
    
    # detect tree tops using the current method
    if (method == 'forest_tools') {
      detected_trees <- calc_ws_func(chm, minHgt = 2)
    } else {
      detected_trees <- lidR::locate_trees(
        chm,
        lidR::lmf(ws = function(x) calc_ws_func(x), shape = 'circular')
      )
      
      # rename "Z" to "height" for consistency
      detected_trees <- detected_trees %>%
        dplyr::rename(height = Z)
    }
    
    # check if detected_trees has valid data
    if (nrow(detected_trees) == 0) {
      cat('No detected trees found for parcel', parcel_idx, 'with method', method, '\n')
      next
    }
    
    ################
    ### Matching ###
    ################
    
    # extract coordinates for both detected and reference trees
    coords_detected <- sf::st_coordinates(detected_trees)
    coords_reference <- sf::st_coordinates(reference_trees)
    
    # initialize empty vectors to store results
    ID_vec_NWFVA <- numeric()       # store matched reference tree IDs
    dist_vec <- numeric()           # store minimum distance
    x_NWFVA <- numeric()            # store x coordinate of nearest reference tree
    y_NWFVA <- numeric()            # store y coordinate of nearest reference tree
    reference_heights <- numeric()  # store the height of the matched reference tree
    
    # loop over each detected tree in the remote sensing dataset
    for (i in 1:nrow(coords_detected)) {
      temp_x <- coords_detected[i, 1] 
      temp_y <- coords_detected[i, 2]
      
      # calculate euclidean distances to all reference trees
      temp_dist_vec <- sqrt((temp_x - coords_reference[, 1])^2 + (temp_y - coords_reference[, 2])^2)
      
      # find the index of the nearest reference tree
      if (length(temp_dist_vec) > 0) {
        min_index <- which.min(temp_dist_vec)
        
        # extract corresponding reference tree information
        possible_treeID <- reference_trees$baumid[min_index]
        ID_vec_NWFVA <- c(ID_vec_NWFVA, possible_treeID)
        dist_vec <- c(dist_vec, min(temp_dist_vec))
        x_NWFVA <- c(x_NWFVA, coords_reference[min_index, 1])
        y_NWFVA <- c(y_NWFVA, coords_reference[min_index, 2])
        reference_heights <- c(reference_heights, as.numeric(reference_trees$baumhoehe[min_index]))
      }
    }
    
    # ensure the number of detected tree IDs matches the length of results
    detected_tree_ids <- detected_trees$treeID
    detected_heights <- detected_trees$height
    detected_geometries <- sf::st_geometry(detected_trees)
    
    # create data frame with the matching results for this parcel
    df <- data.frame(
      NWFVA_Baum_ID = ID_vec_NWFVA,
      NWFVA_height = reference_heights,
      x_NWFVA = x_NWFVA,
      y_NWFVA = y_NWFVA,
      ttop_ID = detected_tree_ids,
      ttop_height = detected_heights,
      ttop_geometry = detected_geometries,
      dist_min = dist_vec
    )
    
    # store the results for the current parcel and method
    method_results[[paste0('Parcel_', parcel_idx)]] <- df
    
    # filter the data frame to exclude rows where NWFVA_height is NA
    parcel_filtered <- df %>%
      dplyr::filter(!is.na(NWFVA_height))
    
    # write data frame to disk
    output_file <- file.path(output_dir, paste0('results_', method, '_parcel_', parcel_idx, '.txt'))
    write.table(parcel_filtered, file = output_file, sep = ';', row.names = F)
    cat('Data written to:', output_file, '\n')
    
    ###################
    ### Correlation ###
    ###################
    
    # calculate pearson correlation and R-squared
    if (nrow(parcel_filtered) > 1) {
      corr <- stats::cor(parcel_filtered$NWFVA_height, parcel_filtered$ttop_height)
      r_squared <- corr^2
      
      # print correlation result
      cat('Correlation for method', method, 'in parcel', parcel_idx, ':', corr, '\n')
      cat('R² for method', method, 'in parcel', parcel_idx, ':', r_squared, '\n')
    } else {
      cat('Not enough valid data for correlation for method', method, 'in parcel', parcel_idx, '\n')
      r_squared <- NA
    }
    
    ################
    ### Plotting ###
    ################
    
    # find the range to ensure both axes have the same scaling
    axis_limits <- range(c(parcel_filtered$NWFVA_height, parcel_filtered$ttop_height), na.rm = TRUE)
    
    # create scatterplot
    plot <- ggplot(parcel_filtered, aes(x = NWFVA_height, y = ttop_height)) +
      geom_point(color = 'black') +
      geom_abline(slope = 1, intercept = 0, color = 'red', lwd = 1) +
      labs(title = paste(method, 'parcel', parcel_idx),
           x = 'Reference tree height (m)',
           y = 'Remote sensing detected height (m)') +
      coord_equal(xlim = axis_limits, ylim = axis_limits) +
      annotate('text', x = min(axis_limits), y = max(axis_limits),
               label = paste0('R² = ', round(r_squared, 2)), 
               hjust = 0, vjust = 1, size = 5)
    
    # store plot for this parcel and method
    method_plots[[paste0('Parcel_', parcel_idx)]] <- plot
  }
  
  # store results and plots for the current method
  results_all_methods[[method]] <- method_results
  plots_all_methods[[method]] <- method_plots
}


# arrange plots in a 3-column by 4-row grid
combined_plots <- cowplot::plot_grid(
  
  # linear method
  plots_all_methods$linear$Parcel_1, plots_all_methods$linear$Parcel_2, plots_all_methods$linear$Parcel_3,
  
  # exponential decay method
  plots_all_methods$exp_decay$Parcel_1, plots_all_methods$exp_decay$Parcel_2, plots_all_methods$exp_decay$Parcel_3,
  
  # spline interpolation method
  plots_all_methods$spline_int$Parcel_1, plots_all_methods$spline_int$Parcel_2, plots_all_methods$spline_int$Parcel_3,
  
  # ForestTools method
  plots_all_methods$forest_tools$Parcel_1, plots_all_methods$forest_tools$Parcel_2, plots_all_methods$forest_tools$Parcel_3,
  
  ncol = 3,
  nrow = 4,
  align = 'hv'
)

# add title for the entire plot
title_plot <- cowplot::ggdraw() + 
  cowplot::draw_label('Comparison of remote sensing and reference tree top heights',
                      fontface = 'bold', size = 16)

# combine title and grid of plots
final_plot <- cowplot::plot_grid(title_plot, combined_plots, ncol = 1, rel_heights = c(0.1, 1))

# save to disk
ggsave(file.path(output_dir, 'comparison_scatterplots.pdf'), final_plot, width = 16, height = 20)

# display final plot
print(final_plot)

# example plot parcel 1
pdf(file.path(output_dir, 'comparison_parcel_1.pdf'), width = 8, height = 12)
par_org <- par()
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, ncol = 2, byrow = T))

terra::plot(chm_parcel_1, main = 'ForestTools')
terra::plot(reference_trees_parcel_1$geometry, pch = 4, add = T)
terra::plot(results_all_methods$forest_tools$Parcel_1$geometry, pch = 4, col = 'white', add = T)
terra::plot(sf::st_geometry(sf::st_as_sf(results_all_methods$forest_tools$Parcel_1, coords = c('x_NWFVA', 'y_NWFVA'), 
                                         crs = sf::st_crs(reference_trees))), pch = 4, col = 'red', add = T)
terra::plot(chm_parcel_1, main = 'linear')
terra::plot(reference_trees_parcel_1$geometry, pch = 4, add = T)
terra::plot(results_all_methods$linear$Parcel_1$geometry, pch = 4, col = 'white', add = T)
terra::plot(sf::st_geometry(sf::st_as_sf(results_all_methods$linear$Parcel_1, coords = c('x_NWFVA', 'y_NWFVA'), 
                                         crs = sf::st_crs(reference_trees))), pch = 4, col = 'red', add = T)
terra::plot(chm_parcel_1, main = 'exp_decay')
terra::plot(reference_trees_parcel_1$geometry, pch = 4, add = T)
terra::plot(results_all_methods$exp_decay$Parcel_1$geometry, pch = 4, col = 'white', add = T)
terra::plot(sf::st_geometry(sf::st_as_sf(results_all_methods$exp_decay$Parcel_1, coords = c('x_NWFVA', 'y_NWFVA'), 
                                         crs = sf::st_crs(reference_trees))), pch = 4, col = 'red', add = T)
terra::plot(chm_parcel_1, main = 'spline_int')
terra::plot(reference_trees_parcel_1$geometry, pch = 4, add = T)
terra::plot(results_all_methods$spline_int$Parcel_1$geometry, pch = 4, col = 'white', add = T)
terra::plot(sf::st_geometry(sf::st_as_sf(results_all_methods$spline_int$Parcel_1, coords = c('x_NWFVA', 'y_NWFVA'), 
                                         crs = sf::st_crs(reference_trees))), pch = 4, col = 'red', add = T)

plot.new()
legend('bottom', 
       legend = c('Reference trees', 'Remote sensing detected trees', 'Matched trees'), 
       pch = c(4, 4, 4), 
       col = c('black', 'white', 'red'), 
       pt.bg = c('black', 'white', 'red'),
       cex = 2,
       bg = 'gray')
dev.off()
par(par_org)



