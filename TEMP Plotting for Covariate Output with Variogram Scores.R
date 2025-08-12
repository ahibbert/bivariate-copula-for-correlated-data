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
  facet_wrap(~score_type_label, scales = "free", ncol = 1, labeller = labeller(score_type_label = label_value)) +
  theme_bw() +
  labs(title = "Model evaluation",
       y = "Score Value", x = "Method") +
  theme(legend.position = "none")

ggsave(paste("Charts/Variogram_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), plot = p, width = 9, height = 15, dpi = 900)


#####Coefficients WITH ALL COEFFICIENTS INCLUDED####
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

ggplot(results_long, aes(x = model, y = value, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free", ncol = 2) +
  geom_hline(
    data = true_vals_df,
    aes(yintercept = true_value),
    colour = "red", linetype = "dashed", linewidth = 1
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Comparison of Models for x1, x1_se, x2, x2_se",
    y = "Value",
    x = "Model"
  )