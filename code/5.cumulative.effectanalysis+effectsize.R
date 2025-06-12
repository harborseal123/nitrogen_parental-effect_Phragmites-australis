# Clear console and remove all variables to reset the environment
cat("\014")
rm(list = ls())

# Load required packages
# dplyr: Data manipulation
# readr: Reading CSV files
# car: ANOVA and statistical tests (e.g., Levene's test)
# emmeans: Post-hoc multiple comparisons
# MASS: Box-Cox transformation
# skimr: Data overview
library(dplyr)
library(readr)
library(car)
library(emmeans)
library(MASS)
library(skimr)

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

# Build linear models for each dependent variable
lm_models <- lapply(vars, function(var) {
  formula <- as.formula(paste(var, "~ N * gene"))
  lm(formula, data = data %>% filter(!is.na(get(var))))
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

# Two-way ANOVA (Type III SS)
aov_models <- lapply(lm_models, function(model) {
  Anova(model, type = "III", singular.ok = TRUE)
})

names(aov_models) <- vars
aov_models

# Custom function to extract ANOVA results into a formatted dataframe
anova_fp <- function(anova, trait) {
  tab <- as.data.frame(anova)
  # 整理rownames，防止空格
  rownames(tab) <- gsub(" ", "", rownames(tab))
  effs <- c("N", "gene", "N:gene")
  df <- data.frame(
    Trait = trait,
    Effect = effs,
    F_value = as.numeric(tab[effs, "F value"]),
    P_value = as.numeric(tab[effs, "Pr(>F)"])
  )
  return(df)
}

# Batch extraction of ANOVA results
res_list <- mapply(anova_fp, aov_models, names(aov_models), SIMPLIFY = FALSE)
anova_summary <- do.call(rbind, res_list)

# Round numeric values
anova_summary$F_value <- round(anova_summary$F_value, 3)
anova_summary$P_value <- round(anova_summary$P_value, 3)

# Export ANOVA results
print(anova_summary)
write.csv(anova_summary, file = "./result/aov_models.cumulative.effect.csv", row.names = FALSE)

# Post-hoc pairwise comparisons (emmeans)
posthoc_results <- lapply(lm_models, function(model) {
  em <- emmeans(model, pairwise ~ N | gene)
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