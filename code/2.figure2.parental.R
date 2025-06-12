# Clear console and remove all variables to reset the environment
cat("\014")
rm(list = ls())

# Load required packages
# readr: Reading CSV files
# dplyr: Data manipulation
# ggplot2: Data visualization
# patchwork: Combining multiple plots
# ggpubr: Enhanced ggplot2 visualizations and statistical annotations
# skimr: Data overview
library(readr) 
library(dplyr) 
library(ggplot2) 
library(patchwork) 
library(ggpubr) 
library(skimr)     

# Set working directory
setwd("C:/Users/qianh/Desktop/R/gene.N")
getwd()

# Read data
data <- read_csv("./result/effectsize.parental.csv")
data <- data %>%
  mutate(gene = paste0("G", as.character(gene)))  # Add "G" prefix to genotypes

# Data overview
skim(data)
str(data)

# Factor design for residual effect
data$gene <- factor(data$gene, levels = c("G12", "G11", "G10", "G9", "G8", "G7", "G6", "G5", "G4", "G3", "G2", "G1"))

# Define dependent variables for analysis
response_variables <- c("leaf_area","node_number", 
          "leafbiomass", "stembiomass", "rootbiomass", "totalbiomass")

# Define the color mapping for each response variable
color_mapping <- list(
  point_color_leaf_area = c("G4", "G5", "G7", "G8", "G9", "G10", "G11", "G12"),
  point_color_node_number = c("G3", "G6", "G8"),
  point_color_leafbiomass = c("G3", "G4", "G6", "G7", "G8", "G10", "G11", "G12"),
  point_color_stembiomass = c("G3", "G4", "G6", "G7", "G8", "G10", "G11"),
  point_color_rootbiomass = c("G1", "G3", "G7", "G8"),
  point_color_totalbiomass = c("G3", "G6", "G7", "G8", "G10")
)

# Apply the color mapping to each response variable
for (col_name in names(color_mapping)) {
  data[[col_name]] <- ifelse(data$gene %in% color_mapping[[col_name]], "#414986", "#c2bddb")
}

# Define significance annotations for each response variable
signif_leaf_area <- data.frame(
  x = rep(-1, 8),
  annotations = c("***", "*", "***", "*", "*", "*", "**", "***"),
  y_position = c(8.8, 7.8, 5.8, 4.8, 3.8, 2.8, 1.8, 0.8))
signif_node_number <- data.frame(
  x = rep(-1, 3),
  annotations = c("**", "*", "**"),
  y_position = c(9.8, 6.8, 4.8))
signif_leafbiomass <- data.frame(
  x = rep(-1, 8),
  annotations = c("***", "***", "***", "***", "***", "***", "**", "**"),
  y_position = c(9.8, 8.8, 6.8, 5.8, 4.8, 2.8, 1.8, 0.8))
signif_stembiomass <- data.frame(
  x = rep(-1, 7),
  annotations = c("**", "**", "**", "***", "*", "**", "*"),
  y_position = c(9.8, 8.8, 6.8, 5.8, 4.8, 2.8, 1.8))
signif_rootbiomass <- data.frame(
  x = rep(-1, 4),
  annotations = c("*", "*", "**", "*"),
  y_position = c(11.8, 9.8, 5.8, 4.8))
signif_totalbiomass <- data.frame(
  x = rep(-1, 5),
  annotations = c("**", "*", "***", "**", "**"),
  y_position = c(9.8, 6.8, 5.8, 4.8, 2.8))

# Define gene levels and positions for plot alignment
gene_levels <- unique(data$gene)
gene_positions <- seq_along(gene_levels)


# Create reusable plot template
plot_common <- function(data, x_var, se_var, color_var, title, tag, y_label = "Genotype", x_limits = c(-1.2, 2.2), y_limits = c(1.1, 11.9), signif_data = NULL) {
  ggplot(data, aes(x = !!sym(x_var), y = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "black", alpha = 0.8) +
    geom_errorbar(aes(xmin = !!sym(x_var) - !!sym(se_var),
                      xmax = !!sym(x_var) + !!sym(se_var)),
                  width = 0.2, color = "#c2bddb", size = 0.8) +
    geom_point(size = 2.5, aes(color = !!sym(color_var))) +
    geom_segment(data = data.frame(gene = gene_levels),
                 aes(x = -1.2, xend = 2.2, y = gene_positions - 0.5, yend = gene_positions - 0.5),
                 color = "black", alpha = 1) +
    coord_cartesian(ylim = y_limits) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_color_identity() +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.title.x = element_blank(),
      plot.margin = margin(0, 0, 8, 0),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.2, 0.975)
    ) +
    labs(title = title, y = y_label, tag = tag) +
    geom_text(
      data = signif_data,
      aes(x = x, y = y_position, label = annotations),
      inherit.aes = FALSE,
      size = 6,
      color = "black"
    )
}

# Generate individual plots 
plot_leaf_area <- plot_common(data, "eff_value_funcleaf_area_mean", "eff_value_funcleaf_area_se", "point_color_leaf_area", "Functional leaf area", "A", signif_data = signif_leaf_area)
plot_node_number <- plot_common(data, "eff_value_node_number_mean", "eff_value_node_number_se", "point_color_node_number", "Node numbers", "B", y_label = NULL, signif_data = signif_node_number) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))  # Adjust panel label position
plot_leafbiomass <- plot_common(data, "eff_value_leafbiomass_mean", "eff_value_leafbiomass_se", "point_color_leafbiomass", "Leaf biomass", "C", signif_data = signif_leafbiomass)
plot_stembiomass <- plot_common(data, "eff_value_stembiomass_mean", "eff_value_stembiomass_se", "point_color_stembiomass", "Stem biomass", "D", y_label = NULL, signif_data = signif_stembiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))  # Adjust panel label position
plot_rootbiomass <- plot_common(data, "eff_value_rootbiomass_mean", "eff_value_rootbiomass_se", "point_color_rootbiomass", "Root biomass", "E", signif_data = signif_rootbiomass)
plot_totalbiomass <- plot_common(data, "eff_value_totalbiomass_mean", "eff_value_totalbiomass_se", "point_color_totalbiomass", "Total biomass", "F", y_label = NULL, signif_data = signif_totalbiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))  # Adjust panel label position

# Combine plots using patchwork 
plot_parental <- (plot_leaf_area + plot_node_number + 
                  plot_leafbiomass + plot_stembiomass + 
                  plot_rootbiomass + plot_totalbiomass) +
                plot_layout(ncol = 2) +
                plot_annotation(
                     caption = "Effect size", 
                     theme = theme(
                       plot.caption = element_text(
                               size = 12, 
                               color = "black", 
                               hjust = 0.5,        
                               margin = margin(t = 5, b = 5) 
                               ))) +
                theme(
                     plot.tag = element_text(size = 12, face = "bold"),
                     plot.margin = margin(b = 20) 
                )
plot_parental

# Export high-resolution image 
ggexport(plot_parental, filename = "./figure/figure2.parental.png",
         width = 1800,
         height = 2600,
         pointsize = 12,
         res = 300)
