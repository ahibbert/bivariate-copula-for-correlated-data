# This script loops through a list of input datasets, loads them, 
# and generates all the plots as in original script.
# Assumes working directory is set so paths work.
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

# ---- Dataset list ----
input_datasets <- list(
   # list(
   #   rdata = "Data/CoefSimData_GA_2_0.5_NA_1_2_1_0.1_1000_100_2025-08-21.RData",
   #   meta = list(dist="GA",a=2,b=0.5,c=NA,mu1=1,mu2=2,x1=1,x2=.1,n=1000)
   # ),
   # list(
   #   rdata = "Data/CoefSimData_GA_0.5_2_NA_1_2_1_0.1_1000_100_2025-08-20.RData",
   #   meta = list(dist="GA",a=.5,b=2,c=NA,mu1=1,mu2=2,x1=1,x2=.1,n=1000)
   # ),
   # list(
   #   rdata = "Data/CoefSimData_GA_1_1_NA_1_2_1_0.1_1000_100_2025-08-21.RData",
   #   meta = list(dist="GA",a=1,b=1,c=NA,mu1=1,mu2=2,x1=1,x2=.1,n=1000)
   # ),
   # list(
   #   rdata = "Data/CoefSimData_NO_1_1_0.9_1_2_1_0.1_1000_100_2025-08-20.RData",
   #   meta = list(dist="NO",a=1,b=1,c=.9,mu1=1,mu2=2,x1=1,x2=0.1,n=1000)
   # ),
   # list(
   #   rdata = "Data/CoefSimData_NO_1_1_0.75_1_2_1_0.1_1000_100_2025-08-20.RData",
   #   meta = list(dist="NO",a=1,b=1,c=.75,mu1=1,mu2=2,x1=1,x2=0.1,n=1000)
   # ),
   # list(
   #   rdata = "Data/CoefSimData_NO_1_1_0.25_1_2_1_0.1_1000_100_2025-08-20.RData",
   #   meta = list(dist="NO",a=1,b=1,c=.25,mu1=1,mu2=2,x1=1,x2=0.1,n=1000)
   # ),
  # list(
  #   rdata= "Data/CoefSimData_PO_NA_0.2_5_0.5_5_1_0.1_1000_2025-08-20.RData",
  #   meta=list(dist="PO", a=NA, b=.2, c=5, mu1=.5, mu2=5, x1=1, x2=0.1, n=1000)
  # )
  #  ,
  #  list(
  #    rdata= "Data/CoefSimData_PO_NA_5_0.2_0.5_0.5_1_0.1_1000_2025-08-20.RData",
  #    meta=list(dist="PO", a=NA, b=5, c=.2, mu1=.5, mu2=.5, x1=1, x2=0.1, n=1000)
  #  ),
  #  list(
  #    rdata= "Data/CoefSimData_PO_NA_0.2_5_0.5_5_1_0.1_1000_2025-08-20.RData",
  #    meta=list(dist="PO", a=NA, b=.2, c=5, mu1=.5, mu2=5, x1=1, x2=0.1, n=1000)
  #  ),
  #  list(
  #    rdata= "Data/CoefSimData_PO_NA_5_0.2_0.5_0.5_1_0.1_1000_2025-08-20.RData",
  #    meta=list(dist="PO", a=NA, b=5, c=.2, mu1=.5, mu2=.5, x1=1, x2=0.1, n=1000)
  #  ),
  # list(
  #   rdata= "Data/CoefSimData_LO_NA_NA_0.25_0.25_0.5_1_0.1_1000_100_2025-08-21.RData",
  #   meta=list(dist="LO", a=NA, b=NA, c=.25, mu1=.25, mu2=.5, x1=1, x2=0.1, n=1000)
  # ),
  # list(
  #   rdata= "Data/CoefSimData_LO_NA_NA_0.75_0.25_0.5_1_0.1_1000_100_2025-08-21.RData",
  #   meta=list(dist="LO", a=NA, b=NA, c=.75, mu1=.25, mu2=.5, x1=1, x2=0.1, n=1000)
  # )
  
  #list(dist="PO", a=NA, b=.2, c=5, mu1=5, mu2=10, x1=1, x2=0.1, n=1000),
  
  # list(
  #   rdata = "Data/CoefSimData_PO_NA_0.2_5_5_10_1_0.1_1000_10_2025-08-21.RData",
  #   meta = list(dist="PO", a=NA, b=.2, c=5, mu1=5, mu2=10, x1=1, x2=0.1, n=1000)
  # ),
  # #list(dist="PO", a=NA, b=5, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000),
  # list(
  #   rdata = "Data/CoefSimData_PO_NA_5_0.2_5_10_1_0.1_1000_10_2025-08-21.RData",
  #   meta = list(dist="PO", a=NA, b=5, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000)
  # ),
  # #list(dist="PO", a=NA, b=.2, c=5, mu1=5, mu2=10, x1=1, x2=0.1, n=1000),
  # list(
  #   rdata = "Data/CoefSimData_PO_NA_5_0.2_5_10_1_0.1_1000_10_2025-08-21.RData",
  #   meta = list(dist="PO", a=NA, b=5, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000)
  # ),
  # #list(dist="PO", a=NA, b=5, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000),
  # list(
  #   rdata = "Data/CoefSimData_PO_NA_0.2_0.2_5_10_1_0.1_1000_10_2025-08-21.RData",
  #   meta = list(dist="PO", a=NA, b=.2, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000)
  # ),
  # #list(dist="PO", a=NA, b=5, c=5, mu1=5, mu2=10, x1=1, x2=0.1, n=1000),
  # list(
  #  rdata = "Data/CoefSimData_PO_NA_0.2_0.2_5_10_1_0.1_1000_10_2025-08-21.RData",
  #  meta = list(dist="PO", a=NA, b=.2, c=.2, mu1=5, mu2=10, x1=1, x2=0.1, n=1000)
  # ),
  # list(dist="PO", a=NA, b=2, c=1, mu1=5, mu2=10, x1=1, x2=0.01, n=1000),
  #list(
  #  rdata= "Data/CoefSimData_PO_NA_2_1_5_10_1_0.01_1000_10_2025-08-21.RData",
  #  meta=list(dist="PO", a=NA, b=2, c=1, mu1=5, mu2=10, x1=1, x2=0.01, n=1000)
  #),
  #list(dist="PO", a=NA, b=1, c=5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
  # list(
  #   rdata= "Data/CoefSimData_PO_NA_1_5_1_2_1_0.01_1000_10_2025-08-21.RData",
  #   meta=list(dist="PO", a=NA, b=1, c=5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
  # ),
  #list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
  list(
    rdata= "Data/CoefSimData_NO_1_1_0.25_1_2_1_0.01_1000_10_2025-08-21.RData",
    meta=list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
  )
)

# ---- Helper functions ----
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
score_type_map <- c(
  es = "Energy Score",
  vs1 = "Variogram Score (p=1)",
  vs2 = "Variogram Score (p=2)",
  vs2_wt = "Variogram Score (p=2, Weighted)",
  vs2_wt_coronly = "Variogram Score (p=2, Correlated Obs Only)",
  logliks = "Log Likelihood"
)
get_true_ses <- function(n, a, b, c, mu1, mu2, dist, x1, x2) {
  sim_out <- simCovariateMLEs(sims = 100, n = n, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2,
                              dist = dist, x1 = x1, x2 = x2, trace = FALSE)
  ses <- sim_out$ses[1:4]
  names(ses) <- c("t1", "t2", "x1", "x2")
  return(ses)
}

for (ds in input_datasets) {
  # ---- Load dataset ----
  load(ds$rdata)
  # attach meta variables
  dist <- ds$meta$dist; a <- ds$meta$a; b <- ds$meta$b; c <- ds$meta$c
  mu1 <- ds$meta$mu1; mu2 <- ds$meta$mu2; x1 <- ds$meta$x1; x2 <- ds$meta$x2; n <- ds$meta$n
  
  # --- Score boxplots ---
  # Exclude logliks models not in other matrices
  ref_models <- colnames(score_items[[which(names(score_items) != "logliks")[1]]])
  score_items$logliks <- score_items$logliks[, ref_models, drop = FALSE]
  
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
  method_label_order <- sapply(method_order, rename_model)
  
  score_long <- score_long %>%
    mutate(method_label = sapply(as.character(method), rename_model),
           method_label = factor(method_label, levels = method_label_order))
  
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
  
  p=ggplot(score_long, aes(x = method_label, y = value, fill = score_type_label)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(data = meds, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~score_type_label, scales = "free", ncol = 2, labeller = labeller(score_type_label = label_value)) +
    theme_bw() +
    labs(title = "Model evaluation",
         y = "Score Value", x = "Method") +
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(paste("Charts/Variogram_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), plot = p, width = 9, height = 12, dpi = 900)
  
  # ---- Bottom 4 panels, outliers removed ----
  bottom_four_labels <- tail(levels(score_long$score_type_label), 4)
  score_long_bottom4 <- score_long %>% filter(score_type_label %in% bottom_four_labels)
  meds_bottom4 <- meds %>% filter(score_type_label %in% bottom_four_labels)
  
  whisker_limits <- score_long_bottom4 %>%
    group_by(score_type_label, method_label) %>%
    summarize(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower = Q1 - 1.5 * IQR,
      upper = Q3 + 1.5 * IQR
    )
  
  score_no_outliers <- score_long_bottom4 %>%
    left_join(whisker_limits, by = c("score_type_label", "method_label")) %>%
    filter(value >= lower & value <= upper)
  
  p2 = ggplot(score_no_outliers, aes(x = method_label, y = value, fill = score_type_label)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
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
  
  print(p2)
  ggsave(paste("Charts/Variogram_p2only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),"_bottom4_no_outliers.png",sep=""), plot = p2, width = 9, height = 9, dpi = 900)
  
  # ---- Coefficient plots: all coefficients ----
  results_list=par_estimates
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
  results_long <- imap_dfr(results_list, function(mat, varname) {
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
    labs(title = "Comparison of Models for x1, x2", y = "Value", x = "Model")
  
  print(q)
  ggsave(plot=q,file=paste("Charts/AllCoefficients_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 12, height = 12, dpi = 900)
  
  # ---- Coefficient plots: only x1, x2, with true lines ----
  ses_out=get_true_ses(n, a, b, c, mu1, mu2, dist, x1, x2)
  true = c(x1,x2,ses_out[c("x1","x2")])
  names(true)=c("x1","x2","x1_se","x2_se")
  
  vars_to_plot <- c("x1", "x1_se", "x2", "x2_se")
  results_long2 <- imap_dfr(results_list[vars_to_plot], function(mat, varname) {
    mat <- mat[, intersect(model_order, colnames(mat)), drop = FALSE]
    mat_df <- as.data.frame(mat)
    mat_df$run <- 1:nrow(mat_df)
    mat_long <- pivot_longer(mat_df, -run, names_to = "model", values_to = "value")
    mat_long$variable <- varname
    return(mat_long)
  })
  results_long2$model <- factor(results_long2$model, levels = model_order,
                                labels = model_labels[model_order])
  true_vals_df <- data.frame(variable = names(true), true_value = as.numeric(true))
  facet_labels <- c(
    x1 = "X1 Coefficient",
    x1_se = "X1 Standard Error",
    x2 = "X2 Coefficient",
    x2_se = "X2 Standard Error"
  )
  coefplot=ggplot(results_long2, aes(x = model, y = value, fill = model)) +
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
  print(coefplot)
  ggsave(plot=coefplot,file=paste("Charts/Coef_x1x2_only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 9, height = 8, dpi = 900)
}