#-----------------------------------------------------------------------
# Name:         n_trees_comp.R
# Description:  Script compares number of trees per ha detected with
#               remote sensing-based individual tree detection (ITD) vs. 
#               terrestrial counted trees in sample plots.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-----------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------

# read stand data from the Kellerwald
stand_data_kw <- read.csv(file.path(processed_data_dir, 'stand_info_new_DBH_too.csv'))

# read PSI data
psi_kw <- read.csv(
  file.path(processed_data_dir, 'psi_kellerwald_summary_statistics.csv'), 
  header = F, 
  skip = 1,
  col.names = c("id1", "Baeume_pro_Parzelle", "nha", "Mittel_BHD_Parzelle", "hoe_mod")
  )


# 02 - data preparation
#-------------------------------------

# aggregate PSI data by plot, 
# calculate mean height and mean number of trees per plot
psi_agg <- psi_kw %>%
  dplyr::group_by(id1) %>%
  dplyr::summarise(
    mean_height = mean(hoe_mod, na.rm = T),
    mean_ntrees = mean(nha, na.rm = T)
  )



# 03 - comparison of detected trees
#-------------------------------------

# combine both datasets
rs_stands <- stand_data_kw %>%
  select(Mean_DBH_Bestand, mean_height, n_trees_per_ha) %>%
  rename(dbh = Mean_DBH_Bestand, height = mean_height, trees_per_ha = n_trees_per_ha) %>%
  mutate(dataset = "stand_data_kw")

psi_plots <- psi_kw %>%
  select(Mittel_BHD_Parzelle, hoe_mod, nha) %>%
  rename(dbh = Mittel_BHD_Parzelle, height = hoe_mod, trees_per_ha = nha) %>%
  mutate(dataset = "psi_kw")

combined_data <- bind_rows(rs_stands, psi_plots)

# define colors for boxplot boxes
plot_colors <- c("psi_kw" = "black", "stand_data_kw" = "grey70")

# boxplot n trees per ha
p <- ggplot(combined_data, aes(x = dataset, y = trees_per_ha, fill = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, 
               fatten = 2) +
  geom_jitter(
    aes(color = dataset),
    position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1),
    size = 1, alpha = 0.5, show.legend = F
  ) +
  labs(
    x = '',
    y = 'Stammzahl pro ha'
  ) +
  scale_fill_manual(
    values = plot_colors,
    labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")
  ) +
  scale_color_manual(
    values = plot_colors
  ) +
  scale_x_discrete(labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")) +
  ylim(0, 2000) +
  theme_minimal() +
  theme(legend.position = "none")

# extract boxplot data and add colored median lines
dat <- ggplot_build(p)$data[[1]]
p + geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), 
                 colour = "#E69F00", linewidth = 1, inherit.aes = F)

ggsave(
  filename = file.path(output_dir, 'boxplot_nha_psi_rs.png'),
  width = 9,
  height = 5
)

# scatterplot mean height vs n trees per ha
ggplot(combined_data, aes(x = height, y = trees_per_ha)) +
  geom_point(aes(color = dataset)) +
  geom_smooth(data = subset(combined_data, dataset == "psi_kw"), 
              aes(x = height, y = trees_per_ha), 
              method = "loess", se = F, color = "red") +
  geom_smooth(data = subset(combined_data, dataset == "stand_data_kw"), 
              aes(x = height, y = trees_per_ha), 
              method = "loess", se = F, color = "blue") +
  annotate("text", x = 8.2, y = 600, label = "PSI", color = "red", size = 4) +
  annotate("text", x = 6, y = 300, label = "Fernerkundung", color = "blue", size = 4) +
  ylim(0, 3000) +
  xlim(5, 40) +
  labs(
    x = 'Plot-/Bestandesmittelhöhe [m]',
    y = 'Stammzahl pro ha',
    color = ""
  ) +
  scale_color_manual(values = c("psi_kw" = "black", "stand_data_kw" = "grey70"),
                     labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(output_dir, 'scatterplot_height_vs_nha_psi_rs.emf'),
  width = 9,
  height = 5
)

ggplot(combined_data, aes(x = height, y = trees_per_ha)) +
  geom_point(aes(color = dataset)) +
  geom_smooth(data = subset(combined_data, dataset == "psi_kw"), 
              aes(x = height, y = trees_per_ha), 
              method = "loess", se = F, color = "red") +
  geom_smooth(data = subset(combined_data, dataset == "stand_data_kw"), 
              aes(x = height, y = trees_per_ha), 
              method = "loess", se = F, color = "blue") +
  annotate("text", x = 7, y = 180, label = "PSI", color = "red", size = 4) +
  annotate("text", x = 7.5, y = 265, label = "Fernerkundung", color = "blue", size = 4) +
  ylim(0, 500) +
  xlim(5, 40) +
  labs(
    x = 'Plot-/Bestandesmittelhöhe [m]',
    y = 'Stammzahl pro ha',
    color = ""
  ) +
  scale_color_manual(values = c("psi_kw" = "black", "stand_data_kw" = "grey70"),
                     labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(output_dir, 'scatterplot_height_vs_nha_psi_rs_ylim.emf'),
  width = 9,
  height = 5
)

# scatterplot mean dbh vs n trees per ha
ggplot(combined_data, aes(x = dbh, y = trees_per_ha)) +
  geom_point(aes(color = dataset)) +
  geom_smooth(data = subset(combined_data, dataset == "psi_kw"), 
              aes(x = dbh, y = trees_per_ha), 
              method = "loess", se = F, color = "red") +
  geom_smooth(data = subset(combined_data, dataset == "stand_data_kw"), 
              aes(x = dbh, y = trees_per_ha), 
              method = "loess", se = F, color = "blue") +
  annotate("text", x = 8.2, y = 1000, label = "PSI", color = "red", size = 4) +
  annotate("text", x = 6, y = 300, label = "Fernerkundung", color = "blue", size = 4) +
  ylim(0, 3000) +
  xlim(5, 40) +
  labs(
    x = 'Plot-/Bestandesmittel-BHD [cm]',
    y = 'Stammzahl pro ha',
    color = ""
  ) +
  scale_color_manual(values = c("psi_kw" = "black", "stand_data_kw" = "grey70"),
                     labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(output_dir, 'scatterplot_dbh_vs_nha_psi_rs.emf'),
  width = 9,
  height = 5
)

ggplot(combined_data, aes(x = dbh, y = trees_per_ha)) +
  geom_point(aes(color = dataset)) +
  geom_smooth(data = subset(combined_data, dataset == "psi_kw"), 
              aes(x = dbh, y = trees_per_ha), 
              method = "loess", se = F, color = "red") +
  geom_smooth(data = subset(combined_data, dataset == "stand_data_kw"), 
              aes(x = dbh, y = trees_per_ha), 
              method = "loess", se = F, color = "blue") +
  annotate("text", x = 7, y = 180, label = "PSI", color = "red", size = 4) +
  annotate("text", x = 7.5, y = 265, label = "Fernerkundung", color = "blue", size = 4) +
  ylim(0, 500) +
  xlim(5, 40) +
  labs(
    x = 'Plot-/Bestandesmittelhöhe [m]',
    y = 'Stammzahl pro ha',
    color = ""
  ) +
  scale_color_manual(values = c("psi_kw" = "black", "stand_data_kw" = "grey70"),
                     labels = c("psi_kw" = "PSI", "stand_data_kw" = "Fernerkundung")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(output_dir, 'scatterplot_height_vs_nha_psi_rs_ylim.emf'),
  width = 9,
  height = 5
)







