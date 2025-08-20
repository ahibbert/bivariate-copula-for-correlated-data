#load("Data/CoefSimData_GA_2_0.5_NA_1_2_1_0.01_1000_2025-08-20.RData");dist="GA";a=2;b=0.5;c=NA;mu1=1;mu2=2;x1=1;x2=.1;n=1000
#load("Data/CoefSimData_GA_0.5_2_NA_1_2_1_0.01_1000_2025-08-20.RData");dist="GA";a=.5;b=2;c=NA;mu1=1;mu2=2;x1=1;x2=.1;n=1000
#load("Data/CoefSimData_GA_1_1_NA_1_2_1_0.01_1000_2025-08-20.RData");dist="GA";a=1;b=1;c=NA;mu1=1;mu2=2;x1=1;x2=.1;n=1000

load("Data/CoefSimData_NO_1_1_0.75_1_2_1_0.01_1000_2025-08-20.RData"); dist="NO";a=1;b=1;c=.75;mu1=1;mu2=2;x1=1;x2=0.1;n=1000
#load("Data/CoefSimData_NO_1_1_0.25_1_2_1_0.01_1000_2025-08-20.RData"); dist="NO";a=1;b=1;c=.25;mu1=1;mu2=2;x1=1;x2=0.1;n=1000

#######VARIOGRAM AND LOGLIK#####
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# Exclude logliks models not in other matrices
ref_models <- colnames(score_items[[which(names(score_items) != "logliks")[1]]])
score_items$logliks <- score_items$logliks[, ref_models, drop = FALSE]

# Convert each matrix in the list to a long-format data frame, combine into one data frame
score_long <- map2_dfr(score_items, names(score_items), ~ {
  df <- as.data.frame(.x)
  df$row <- seq_len(nrow(df))
  df_long <- pivot_longer(df, -row, names_to = "method", values_to = "value")
  df_long$score_type <- .y
  df_long
})

#score_long <- score_long %>% select(-row)

# Set the correct order and labels for method
method_order <- c(
  "glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
  "cop", "cop_n", "cop_j", "cop_g", "cop_f",
  "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t"
)

rename_model <- function(x) {
  main_map <- c(
    glm = "GLM",
    gee = "GEE",
    re_nosig = "GAMLSS",
    re_np = "GAMLSS NP",
    lme4 = "LME4",
    gamm = "GAMM"
  )
  cop_map <- c(
    cop = "GJRM (C)",
    cop_n = "GJRM (N)",
    cop_j = "GJRM (J)",
    cop_g = "GJRM (G)",
    cop_f = "GJRM (F)",
    cop_amh = "GJRM (AMH)",
    cop_fgm = "GJRM (FGM)",
    cop_pl = "GJRM (PL)",
    cop_h = "GJRM (H)",
    cop_t = "GJRM (T)"
  )
  if(x %in% names(main_map)) return(main_map[x])
  if(x %in% names(cop_map)) return(cop_map[x])
  return(x)
}

method_label_order <- sapply(method_order, rename_model)

score_long <- score_long %>%
  mutate(method_label = sapply(as.character(method), rename_model),
         method_label = factor(method_label, levels = method_label_order))

# Set order and labels for score_type
score_type_map <- c(
  es = "Energy Score",
  vs1 = "Variogram Score (p=1)",
  vs2 = "Variogram Score (p=2)",
  vs2_wt = "Variogram Score (p=2, Weighted)",
  vs2_wt_coronly = "Variogram Score (p=2, Correlated Obs Only)",
  logliks = "Log Likelihood"
)
score_type_order <- score_type_map[c("es", "vs1", "vs2", "vs2_wt", "vs2_wt_coronly", "logliks")]

score_long <- score_long %>%
  mutate(score_type_label = score_type_map[score_type],
         score_type_label = factor(score_type_label, levels = score_type_order))

# Calculate lowest median per score_type
meds <- score_long %>%
  group_by(score_type, method_label) %>%
  summarize(med = median(value), .groups = "drop") %>%
  group_by(score_type) %>%
  summarize(min_median = min(med), .groups = "drop") %>%
  mutate(score_type_label = score_type_map[score_type],
         score_type_label = factor(score_type_label, levels = score_type_order))

# Plot (legend removed)
p=ggplot(score_long, aes(x = method_label, y = value, fill = score_type_label)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(data = meds, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
  facet_wrap(~score_type_label, scales = "free", ncol = 2, labeller = labeller(score_type_label = label_value)) +
  theme_bw() +
  labs(title = "Model evaluation",
       y = "Score Value", x = "Method") +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(paste("Charts/Variogram_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), plot = p, width = 9, height = 15, dpi = 900)



#### All four bottom panels in one plot, outliers removed, in one chart ####
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# Exclude logliks models not in other matrices
ref_models <- colnames(score_items[[which(names(score_items) != "logliks")[1]]])
score_items$logliks <- score_items$logliks[, ref_models, drop = FALSE]

# Convert each matrix in the list to a long-format data frame, combine into one data frame
score_long <- map2_dfr(score_items, names(score_items), ~ {
  df <- as.data.frame(.x)
  df$row <- seq_len(nrow(df))
  df_long <- pivot_longer(df, -row, names_to = "method", values_to = "value")
  df_long$score_type <- .y
  df_long
})

method_order <- c(
  "glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
  "cop", "cop_n", "cop_j", "cop_g", "cop_f",
  "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t"
)
rename_model <- function(x) {
  main_map <- c(
    glm = "GLM",
    gee = "GEE",
    re_nosig = "GAMLSS",
    re_np = "GAMLSS NP",
    lme4 = "LME4",
    gamm = "GAMM"
  )
  cop_map <- c(
    cop = "GJRM (C)",
    cop_n = "GJRM (N)",
    cop_j = "GJRM (J)",
    cop_g = "GJRM (G)",
    cop_f = "GJRM (F)",
    cop_amh = "GJRM (AMH)",
    cop_fgm = "GJRM (FGM)",
    cop_pl = "GJRM (PL)",
    cop_h = "GJRM (H)",
    cop_t = "GJRM (T)"
  )
  if(x %in% names(main_map)) return(main_map[x])
  if(x %in% names(cop_map)) return(cop_map[x])
  return(x)
}
method_label_order <- sapply(method_order, rename_model)
score_long <- score_long %>%
  mutate(method_label = sapply(as.character(method), rename_model),
         method_label = factor(method_label, levels = method_label_order))

score_type_map <- c(
  es = "Energy Score",
  vs1 = "Variogram Score (p=1)",
  vs2 = "Variogram Score (p=2)",
  vs2_wt = "Variogram Score (p=2, Weighted)",
  vs2_wt_coronly = "Variogram Score (p=2, Correlated Obs Only)",
  logliks = "Log Likelihood"
)
score_type_order <- score_type_map[c("es", "vs1", "vs2", "vs2_wt", "vs2_wt_coronly", "logliks")]
score_long <- score_long %>%
  mutate(score_type_label = score_type_map[score_type],
         score_type_label = factor(score_type_label, levels = score_type_order))

meds <- score_long %>%
  group_by(score_type, method_label) %>%
  summarize(med = median(value), .groups = "drop") %>%
  group_by(score_type) %>%
  summarize(min_median = min(med), .groups = "drop") %>%
  mutate(score_type_label = score_type_map[score_type],
         score_type_label = factor(score_type_label, levels = score_type_order))

bottom_four_labels <- tail(levels(score_long$score_type_label), 4)
score_long_bottom4 <- score_long %>% filter(score_type_label %in% bottom_four_labels)
meds_bottom4 <- meds %>% filter(score_type_label %in% bottom_four_labels)

# Calculate whisker limits (exclude outliers) for each facet/panel and method
whisker_limits <- score_long_bottom4 %>%
  group_by(score_type_label, method_label) %>%
  summarize(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 1.5 * IQR,
    upper = Q3 + 1.5 * IQR
  )

# Join limits back and filter out outliers
score_no_outliers <- score_long_bottom4 %>%
  left_join(whisker_limits, by = c("score_type_label", "method_label")) %>%
  filter(value >= lower & value <= upper)

# Plot as usual, but outliers are gone
p = ggplot(score_no_outliers, aes(x = method_label, y = value, fill = score_type_label)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) + # outlier.shape=NA disables drawing dots just in case
  geom_hline(data = meds_bottom4, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
  facet_wrap(~score_type_label, scales = "free", ncol = 2, labeller = labeller(score_type_label = label_value)) +
  theme_bw() +
  labs(
    title = "Model evaluation (Bottom Four Score Types, outliers removed)",
    y = "Score Value", x = "Method"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p
ggsave(paste("Charts/Variogram_p2only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),"_bottom4_no_outliers.png",sep=""), plot = p, width = 9, height = 12, dpi = 900)

####Coefficients WITH ALL COEFFICIENTS INCLUDED####
results_list=par_estimates
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Desired model order
model_order <- c(
  "glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
  "cop_n", "cop", "cop_j", "cop_g", "cop_f"
)

# Model label maps
main_map <- c(
  glm = "GLM",
  gee = "GEE",
  re_nosig = "GAMLSS",
  re_np = "GAMLSS NP",
  lme4 = "LME4",
  gamm = "GAMM"
)
cop_map <- c(
  cop = "GJRM (C)",
  cop_n = "GJRM (N)",
  cop_j = "GJRM (J)",
  cop_g = "GJRM (G)",
  cop_f = "GJRM (F)",
  cop_amh = "GJRM (AMH)",
  cop_fgm = "GJRM (FGM)",
  cop_pl = "GJRM (PL)",
  cop_h = "GJRM (H)",
  cop_t = "GJRM (T)"
)

# Combine maps for relabeling
model_labels <- c(main_map, cop_map)

# Tidy data creation, restrict to desired models
results_long <- imap_dfr(results_list, function(mat, varname) {
  # Keep only desired models, in the given order
  mat <- mat[, intersect(model_order, colnames(mat)), drop = FALSE]
  mat_df <- as.data.frame(mat)
  mat_df$run <- 1:nrow(mat_df)
  mat_long <- pivot_longer(mat_df, -run, names_to = "model", values_to = "value")
  mat_long$variable <- varname
  return(mat_long)
})

results_long$model <- factor(results_long$model, levels = model_order,
                             labels = model_labels[model_order])

q=ggplot(results_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Comparison of Models Across Variables", y = "Value", x = "Model")

q
ggsave(plot=q,file=paste("Charts/AllCoefficients_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 12, height = 12, dpi = 900)
#####Coefficients WITH ONLY x1, x2####

get_true_ses <- function(n, a, b, c, mu1, mu2, dist, x1, x2) {
  # sims=100 or 200 recommended for stable SE estimate
  sim_out <- simCovariateMLEs(sims = 200, n = n, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2,
                              dist = dist, x1 = x1, x2 = x2, trace = FALSE)
  # The returned ses vector is: c("mu1","mu2","x1","x2","s1","s2") -- only first 4 are for parameters of interest
  ses <- sim_out$ses[1:4]
  names(ses) <- c("t1", "t2", "x1", "x2")
  return(ses)
}

ses_out=get_true_ses(n, a, b, c, mu1, mu2, dist, x1, x2)

true = c(x1,x2,ses_out[c("x1","x2")])
names(true)=c("x1","x2","x1_se","x2_se")

results_list=par_estimates
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Desired model order and label maps (same as before)
model_order <- c(
  "glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
  "cop_n", "cop", "cop_j", "cop_g", "cop_f"
)
main_map <- c(
  glm = "GLM",
  gee = "GEE",
  re_nosig = "GAMLSS",
  re_np = "GAMLSS NP",
  lme4 = "LME4",
  gamm = "GAMM"
)
cop_map <- c(
  cop = "GJRM (C)",
  cop_n = "GJRM (N)",
  cop_j = "GJRM (J)",
  cop_g = "GJRM (G)",
  cop_f = "GJRM (F)",
  cop_amh = "GJRM (AMH)",
  cop_fgm = "GJRM (FGM)",
  cop_pl = "GJRM (PL)",
  cop_h = "GJRM (H)",
  cop_t = "GJRM (T)"
)
model_labels <- c(main_map, cop_map)

vars_to_plot <- c("x1", "x1_se", "x2", "x2_se")

results_long <- imap_dfr(results_list[vars_to_plot], function(mat, varname) {
  mat <- mat[, intersect(model_order, colnames(mat)), drop = FALSE]
  mat_df <- as.data.frame(mat)
  mat_df$run <- 1:nrow(mat_df)
  mat_long <- pivot_longer(mat_df, -run, names_to = "model", values_to = "value")
  mat_long$variable <- varname
  return(mat_long)
})

results_long$model <- factor(results_long$model, levels = model_order,
                             labels = model_labels[model_order])

# true is assumed to be a named numeric vector like:
# true <- c(x1 = 1, x2 = 1, x1_se = 0.19795725, x2_se = 0.03757125)
true_vals_df <- data.frame(variable = names(true), true_value = as.numeric(true))

# Facet labels
facet_labels <- c(
  x1 = "X1 Coefficient",
  x1_se = "X1 Standard Error",
  x2 = "X2 Coefficient",
  x2_se = "X2 Standard Error"
)

coefplot=ggplot(results_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free", ncol = 2, 
             labeller = as_labeller(facet_labels)) +
  geom_hline(
    data = true_vals_df,
    aes(yintercept = true_value),
    colour = "red", linetype = "dashed", linewidth = 1
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Coefficient Estimates and Standard Error for Each Model (10 Runs)",
    y = "Value",
    x = "Model"
  ) +
  theme(legend.position = "none")

coefplot
ggsave(plot=coefplot,file=paste("Charts/Coef_x1x2_only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 9, height = 8, dpi = 900)


