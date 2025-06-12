# Clear memory
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
data <- read_csv("./result/effectsize.cumulative.effect.csv")
data <- data %>%
  mutate(gene = paste0("G", as.character(gene)))

# Data overview
skim(data)
str(data)

# Factor design for residual effect
data$gene <- factor(data$gene, levels = c("G12", "G11", "G10", "G9", "G8", "G7", "G6", "G5", "G4", "G3","G2", "G1"))

# Define dependent variables for analysis
response_variables <- c("node_number", "underbiomass", "abovebiomass", "totalbiomass")

# Define the color mapping for each response variable
color_mapping <- list(
  point_color_node_number = c("G8", "G9"),
  point_color_underbiomass = c("G8"),
  point_color_abovebiomass = c("G12"),
  point_color_totalbiomass = NULL
)

# Apply the color mapping to each response variable
for (col_name in names(color_mapping)) {
  data[[col_name]] <- ifelse(data$gene %in% color_mapping[[col_name]], "#414986", "#c2bddb")
}

# Define significance annotations for each response variable
signif_node_number <- data.frame(
  x = rep(-1, 2), 
  annotations = c("**", "*"),
  y_position = c(4.8, 3.8))
signif_underbiomass <- data.frame(
  x = -1, 
  annotations = "**",
  y_position = 4.8)
signif_abovebiomass <- data.frame(
  x = -1, 
  annotations = "*",
  y_position = 0.8)

# Define gene levels and positions for plot alignment
gene_levels <- unique(data$gene)
gene_positions <- seq_along(gene_levels)


# Create reusable plot template
plot_common <- function(data, x_var, se_var, color_var, title, tag, y_label = "Genotype", x_limits = c(-1.2, 3.2), y_limits = c(1.1, 11.9), signif_data = NULL) {
  p <- ggplot(data, aes(x = !!sym(x_var), y = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 1, color = "black", alpha = 0.8) +
    geom_errorbar(aes(xmin = !!sym(x_var) - !!sym(se_var),
                      xmax = !!sym(x_var) + !!sym(se_var)),
                  width = 0.2, color = "#c2bddb", size = 0.8) +
    geom_point(size = 2.5, aes(color = !!sym(color_var))) +
    geom_segment(data = data.frame(gene = gene_levels),
                 aes(x = -1.2, xend = 3.2, y = gene_positions - 0.5, yend = gene_positions - 0.5),
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
    labs(title = title, y = y_label, tag = tag) # 只在有显著性数据时添加标注
if (!is.null(signif_data)) {
  p <- p + geom_text(
    data = signif_data,
    aes(x = x, y = y_position, label = annotations),
    inherit.aes = FALSE,
    size = 6,
    color = "black"
  )
}

return(p)
}


# Generate individual plots 
plot_node_number <- plot_common(data, "eff_value_node_number_mean", "eff_value_node_number_se", "point_color_node_number", "Node numbers", "A", signif_data = signif_node_number)
plot_underbiomass <- plot_common(data, "eff_value_underbiomass_mean", "eff_value_underbiomass_se", "point_color_underbiomass", "Below-ground biomass", "B", y_label = NULL, signif_data = signif_underbiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))  # Adjust panel label position
plot_abovebiomass <- plot_common(data, "eff_value_abovebiomass_mean", "eff_value_abovebiomass_se", "point_color_abovebiomass", "Aboveground biomass", "C", signif_data = signif_abovebiomass)
plot_totalbiomass <- plot_common(data, "eff_value_totalbiomass_mean", "eff_value_totalbiomass_se", "point_color_totalbiomass", "Total biomass", "D", y_label = NULL, signif_data = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))  # Adjust panel label position

# Combine plots using patchwork 
plot_cumulative.effect <- (plot_node_number + plot_underbiomass + 
                  plot_abovebiomass + plot_totalbiomass +
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
                     plot.margin = margin(b = 10) 
                ))
plot_cumulative.effect

# Export high-resolution image 
ggexport(plot_cumulative.effect, filename = "./figure/figure4.cumulative.effect.png",
         width = 1800,
         height = 1766,
         pointsize = 12,
         res = 300)
