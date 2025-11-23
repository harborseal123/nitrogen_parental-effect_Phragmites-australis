# Clear console and remove all variables to reset the environment
cat("\014")
rm(list = ls())

# Load required packages
library(dplyr)
library(readr)
library(car)
library(emmeans)
library(MASS)
library(skimr)
library(lme4)

# Set working directory
setwd("C:/Users/qianh/Desktop/R/gene.N")
getwd()

# Read data
data <- read_csv("./data/parentaldata.csv")

# Data overview
skim(data)
str(data)

# Factorize categorical variables
# Nir: Nitrogen treatment with levels "C" (Control) and "N" (Nitrogen)
# gene: Genotype with 12 levels (1-12)
data$Nir <- factor(data$Nir, levels = c("C", "N"))
data$gene <- factor(data$gene, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

# Apply Box-Cox transformation to improve normality for specific variables (totol biomass and number of stolonsegments)
# Box-Cox identifies optimal lambda via maximum likelihood estimation
boxcox_out <- boxcox(lm(totalbiomass ~ Nir * gene, data = data))
boxcox_out <- boxcox(lm(node_number ~ Nir * gene, data = data))

# Transform variables using Box-Cox lambdas
data <- data %>%
  mutate(
    totalbiomass = (totalbiomass^0.8 - 1) / 0.8,
    node_number = (node_number^0.43 - 1) / 0.43
  )

# Define dependent variables for analysis
vars <- c("funcleaf_area", "node_number", "leafbiomass", "stembiomass", "rootbiomass", "totalbiomass")

# Homogeneity of variance test (Levene's Test)
levene_tests <- lapply(vars, function(var) {
  formula <- as.formula(paste(var, "~ Nir * gene"))
  leveneTest(formula, data = data %>% filter(!is.na(get(var))))
})
names(levene_tests) <- vars

# Extract Levene's test p-values
levene_p_values <- sapply(levene_tests, function(x) x$`Pr(>F)`[1])

# Create Levene's test results table
output_leveneTest <- data.frame(
  variable = vars,
  p_value = levene_p_values
)
output_leveneTest

# Build mixed-effects models for each dependent variable
lm_models <- lapply(vars, function(var) {
  formula <- as.formula(paste(var, "~ Nir * gene + (1|block)"))
  lmer(formula, data = data %>% filter(!is.na(get(var))))
})

# Residual normality check (Shapiro-Wilk Test)
shapiro_results <- lapply(lm_models, function(model) {
  shapiro.test(residuals(model))
})

# Organize and export Shapiro-Wilk test results
shapiro_p <- sapply(shapiro_results, function(x) x$p.value)
shapiro_summary <- data.frame(
  variable = vars,
  shapiro_p = shapiro_p
)
print(shapiro_summary)
write.csv(shapiro_summary, file = "./result/shapiro_summary.parental.csv", row.names = FALSE)

# Type III Wald chi-square tests for mixed models
aov_models <- lapply(lm_models, function(model) {
  Anova(model, type = "III", singular.ok = TRUE)
})
names(aov_models) <- vars
aov_models

# Custom function to extract ANOVA results into a formatted dataframe
anova_fp <- function(anova, trait) {
  tab <- as.data.frame(anova)
  rownames(tab) <- gsub(" ", "", rownames(tab))
  effs <- c("Nir", "gene", "Nir:gene")
  df <- data.frame(
    Trait = trait,
    Effect = effs,
    df = as.numeric(tab[effs, "Df"]),
    Chisq = as.numeric(tab[effs, "Chisq"]),
    P_value = as.numeric(tab[effs, "Pr(>Chisq)"])
  )
  return(df)
}

# Batch extraction of ANOVA results
res_list <- mapply(anova_fp, aov_models, names(aov_models), SIMPLIFY = FALSE)
anova_summary <- do.call(rbind, res_list)

# Round numeric values
anova_summary$Chisq <- round(anova_summary$Chisq, 3)
anova_summary$P_value <- round(anova_summary$P_value, 3)

# Export ANOVA results
print(anova_summary)
write.csv(anova_summary, file = "./result/aov_models.parental.csv", row.names = FALSE)

# Extract Sample sizes (n)
lmer_models <- lapply(vars, function(var) {
  fm <- as.formula(paste(var, "~ Nir * gene + (1|block)"))
  lmer(fm, data = filter(data, !is.na(.data[[var]])))
})
names(lmer_models) <- vars
print(lmer_models)

# Post-hoc pairwise comparisons (emmeans)
posthoc_results <- lapply(lm_models, function(model) {
  em <- emmeans(model, pairwise ~ Nir | gene,
                lmer.df = "satterthwaite", 
                adjust  = "holm")
  contrasts <- as.data.frame(summary(em$contrasts))
  result <- contrasts %>% 
    dplyr::select(contrast, t.ratio, p.value) %>% 
    rename(T_value = t.ratio, P_value = p.value) %>%
    mutate(
      Signif = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE            ~ ""
      ),
      T_value = round(as.numeric(T_value), 3),
      P_value = round(as.numeric(P_value), 3))
  return(result)
})
names(posthoc_results) <- vars
posthoc_results

# Export post-hoc results
write.csv(posthoc_results, file = "./result/posthoc_results.parental.csv", row.names = FALSE)

# Calculate control (C) means for each genotype
mean_C <- data %>%
  filter(Nir == "C") %>%
  group_by(gene) %>%
  summarise(across(all_of(vars), mean, na.rm = TRUE, .names = "mean_{.col}"), .groups = "drop")

# Merge control means with original data
data_with_mean <- data %>%
  left_join(mean_C, by = "gene")

# Compute effect sizes (N response relative to C)
effect_values <- data_with_mean %>%
  filter(Nir == "N") %>%
  mutate(across(all_of(vars), 
                ~ (. - get(paste0("mean_", cur_column()))) / 
                  get(paste0("mean_", cur_column())), 
                .names = "eff_value_{.col}"))

# Summarize effect sizes (mean ± SE)
effect_summary <- effect_values %>%
  group_by(gene) %>%
  summarise(across(starts_with("eff_value_"), 
                   list(mean = ~ mean(., na.rm = TRUE), 
                        se = ~ sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))), 
                   .names = "{.col}_{.fn}"), 
            .groups = "drop")
print(effect_summary)

# Export effect size results
write.csv(effect_summary, file = "./result/effectsize.parental.csv", row.names = FALSE)

# ================================================================
# ============================Graphing============================
# ================================================================
# Clear console and remove all variables to reset the environment
cat("\014")
rm(list = ls())

# Load required packages
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
  x = rep(-1, 8),
  annotations = c("**", "**", "**", "***", "*", "**", "*", "*"),
  y_position = c(9.8, 8.8, 6.8, 5.8, 4.8, 2.8, 1.8, 0.8))
signif_rootbiomass <- data.frame(
  x = rep(-1, 4),
  annotations = c("*", "*", "**", "*"),
  y_position = c(11.8, 9.8, 5.8, 4.8))
signif_totalbiomass <- data.frame(
  x = rep(-1, 6),
  annotations = c("***", "**", "***", "**", "**", "*"),
  y_position = c(9.8, 6.8, 5.8, 4.8, 2.8, 0.8))

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
plot_node_number <- plot_common(data, "eff_value_node_number_mean", "eff_value_node_number_se", "point_color_node_number", "Number of rhizome internodes", "B", y_label = NULL, signif_data = signif_node_number) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))
plot_leafbiomass <- plot_common(data, "eff_value_leafbiomass_mean", "eff_value_leafbiomass_se", "point_color_leafbiomass", "Leaf biomass", "C", signif_data = signif_leafbiomass)
plot_stembiomass <- plot_common(data, "eff_value_stembiomass_mean", "eff_value_stembiomass_se", "point_color_stembiomass", "Stem biomass", "D", y_label = NULL, signif_data = signif_stembiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))
plot_rootbiomass <- plot_common(data, "eff_value_rootbiomass_mean", "eff_value_rootbiomass_se", "point_color_rootbiomass", "Belowground biomass", "E", signif_data = signif_rootbiomass)
plot_totalbiomass <- plot_common(data, "eff_value_totalbiomass_mean", "eff_value_totalbiomass_se", "point_color_totalbiomass", "Total biomass", "F", y_label = NULL, signif_data = signif_totalbiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))

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
