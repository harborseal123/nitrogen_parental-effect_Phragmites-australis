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
data <- read_csv("./data/cumulative.effect.data.csv")

# Data overview
skim(data)
str(data)

# Factorize categorical variables
# N: Nitrogen treatment with levels "CN" (parental control and offspring nitrogen) 
#                               and "NN" (parental nitrogen and offspring nitrogen)
# gene: Genotype with 12 levels (1-12)
data$N <- factor(data$N, levels = c("CN", "NN"))
data$gene <- factor(data$gene, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

# Remove NA values in N column
# Create log-transformed underbiomass variable to improve normality for total 
data <- data %>% 
  filter(!is.na(N)) %>%
  mutate(
    underbiomass = log(data$underbiomass)
  )

# Define dependent variables for analysis
vars <- c("node_number", "underbiomass", "abovebiomass", "totalbiomass")

# Homogeneity of variance test (Levene's Test)
levene_tests <- lapply(vars, function(var) {
  formula <- as.formula(paste(var, "~ N * gene"))
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
  formula <- as.formula(paste(var, "~ N * gene + (1|block)"))
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
write.csv(shapiro_summary, file = "./result/shapiro_summary.cumulative.effect.csv", row.names = FALSE)

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
  effs <- c("N", "gene", "N:gene")
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
write.csv(anova_summary, file = "./result/aov_models.cumulative.effect.csv", row.names = FALSE)

# Extract Sample sizes (n)
lmer_models <- lapply(vars, function(var) {
  fm <- as.formula(paste(var, "~ N * gene + (1|block)"))
  lmer(fm, data = filter(data, !is.na(.data[[var]])))
})
names(lmer_models) <- vars
print(lmer_models)

# Post-hoc pairwise comparisons (emmeans)
posthoc_results <- lapply(lm_models, function(model) {
  em <- emmeans(model, pairwise ~ N | gene,
                lmer.df = "satterthwaite", 
                adjust  = "holm")
  contrasts <- as.data.frame(summary(em$contrasts))
  result <- contrasts %>% 
    dplyr::select(contrast, t.ratio, p.value) %>% 
    rename(T_value = t.ratio, P_value = p.value)%>%
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
write.csv(posthoc_results, file = "./result/posthoc_results.cumulative.effect.csv", row.names = FALSE)

# Calculate control (CN) means for each genotype
mean_C <- data %>%
  filter(N == "CN") %>%
  group_by(gene) %>%
  summarise(across(all_of(vars), mean, na.rm = TRUE, .names = "mean_{.col}"), .groups = "drop")

# Merge control means with original data
data_with_mean <- data %>%
  left_join(mean_C, by = "gene")

# Compute effect sizes (NN response relative to CN)
effect_values <- data_with_mean %>%
  filter(N == "NN") %>%
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
write.csv(effect_summary, file = "./result/effectsize.cumulative.effect.csv", row.names = FALSE)

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
data <- read_csv("./result/effectsize.cumulative.effect.csv")
data <- data %>%
  mutate(gene = paste0("G", as.character(gene)))

# Data overview
skim(data)
str(data)

# Factor design for residual effect
data$gene <- factor(data$gene, levels = c("G12", "G11", "G10", "G9", "G8", "G7", "G6", "G5", "G4", "G3","G2", "G1"))
# For visualization purposes, 
# the effect size of belowground biomass for genotype G9 (10.04±1.97) was truncated and displayed as 2.5±0.5 in the figure.
data <- data %>%
  mutate(
    eff_value_underbiomass_mean = ifelse(gene == "G9", 2.5, eff_value_underbiomass_mean),
    eff_value_underbiomass_se   = ifelse(gene == "G9", 0.5, eff_value_underbiomass_se)
  )

# Define dependent variables for analysis
response_variables <- c("node_number", "underbiomass", "abovebiomass", "totalbiomass")

# Define the color mapping for each response variable
color_mapping <- list(
  point_color_node_number = c("G8", "G9"),
  point_color_underbiomass = c("G8", "G9"),
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
  x = rep(-1, 2), 
  annotations = c("*", "*"),
  y_position = c(4.8, 3.8))
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
    labs(title = title, y = y_label, tag = tag)
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
plot_node_number <- plot_common(data, "eff_value_node_number_mean", "eff_value_node_number_se", "point_color_node_number", "Number of rhizome internodes", "A", signif_data = signif_node_number)
plot_underbiomass <- plot_common(data, "eff_value_underbiomass_mean", "eff_value_underbiomass_se", "point_color_underbiomass", "Belowground biomass", "B", y_label = NULL, signif_data = signif_underbiomass) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))
plot_abovebiomass <- plot_common(data, "eff_value_abovebiomass_mean", "eff_value_abovebiomass_se", "point_color_abovebiomass", "Aboveground biomass", "C", signif_data = signif_abovebiomass)
plot_totalbiomass <- plot_common(data, "eff_value_totalbiomass_mean", "eff_value_totalbiomass_se", "point_color_totalbiomass", "Total biomass", "D", y_label = NULL, signif_data = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.tag.position = c(0.05, 0.975))

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
