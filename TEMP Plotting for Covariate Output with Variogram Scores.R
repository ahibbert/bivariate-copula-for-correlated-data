# This script loops through a list of input datasets, loads them,
# and generates all the plots as in original script.
# Assumes working directory is set so paths work.
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
source("common_functions.R")
source("data_file_selector.R")  # Load the file selector
#source("start_httpgd.R")

# ---- Dataset Selection ----
cat("=== Data File Selection ===\n")
cat("Choose your data selection method:\n")
cat("1. Interactive selection (recommended)\n")
cat("2. Auto-select all files\n")
cat("3. Manual specification (original method)\n\n")

selection_method <- readline("Enter choice (1, 2, or 3): ")

if (selection_method == "1") {
  # Interactive selection
  input_datasets <- select_data_files()
  if (is.null(input_datasets)) {
    stop("No files selected. Exiting.")
  }
} else if (selection_method == "2") {
  # Auto-select all files
  input_datasets <- auto_generate_input_datasets()
  if (length(input_datasets) == 0) {
    stop("No data files found. Exiting.")
  }
} else {
  # Manual specification (original method)
  cat("Using manual specification...\n")
  input_datasets <- list(
    list(
      rdata= "Data/CoefSimData_NO_1_1_0.25_1_2_1_0.01_1000_10_2025-08-21.RData",
      meta=list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
    )
  )
}

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
  # Create a unique cache key from all parameters
  cache_key <- paste(n, a, b, c, mu1, mu2, dist, x1, x2, sep="_")
  cache_file <- paste0("Cache/true_ses_", cache_key, ".rds")

  # Create cache directory if it doesn't exist
  if (!dir.exists("Cache")) {
    dir.create("Cache", recursive = TRUE)
  }

  # Check if cached result exists
  if (file.exists(cache_file)) {
    cat("Loading cached true SEs for parameters:", cache_key, "\n")
    return(readRDS(cache_file))
  }

  # If not cached, compute and save
  cat("Computing true SEs for parameters:", cache_key, "(this may take a while...)\n")
  sim_out <- simCovariateMLEs(sims = 1000, n = n, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2, dist = dist, x1 = x1, x2 = x2, trace = FALSE)
  ses <- sim_out$ses[1:4]
  names(ses) <- c("t1", "t2", "x1", "x2")

  # Save to cache
  saveRDS(ses, cache_file)
  cat("Cached true SEs to:", cache_file, "\n")

  return(ses)
}

times_summary_all <- list()
for (ds in input_datasets) {

  # ---- Load dataset ----
  load(ds$rdata)

  #Create a summary of mean and standard deviation of times by column and save into a matrix
  times_summary <- matrix(nrow = ncol(times), ncol = 3)
  colnames(times_summary) <- c("mean", "sd", "conv")
  for (i in 1:ncol(times)) {
    times_summary[i, "mean"] <- mean(times[, i], na.rm = TRUE)
    times_summary[i, "sd"] <- sd(times[, i], na.rm = TRUE)
    times_summary[i, "conv"] <- conv[i]
  }

  rownames(times_summary) <- colnames(times)

  meta_string <- paste(unlist(ds$meta), collapse = "_")

  times_summary_all[[meta_string]] <- times_summary

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
    labs(title = paste("Model evaluation -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"),
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
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(data = meds_bottom4, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~score_type_label, scales = "free", ncol = 2, labeller = labeller(score_type_label = label_value)) +
    theme_bw() +
    labs(
      title = paste("Model evaluation (Bottom Four Score Types, outliers removed) -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"),
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
    labs(title = paste("Comparison of Models for x1, x2 -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"), y = "Value", x = "Model")

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

  # Remove extreme outliers for better visualization
  coef_whisker_limits <- results_long2 %>%
    group_by(variable, model) %>%
    summarize(
      Q1 = quantile(value, 0.25, na.rm = TRUE),
      Q3 = quantile(value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower = Q1 - 1.5 * IQR,
      upper = Q3 + 1.5 * IQR,
      .groups = "drop"
    )

  results_long2_no_outliers <- results_long2 %>%
    left_join(coef_whisker_limits, by = c("variable", "model")) %>%
    filter(value >= lower & value <= upper)

  coefplot=ggplot(results_long2_no_outliers, aes(x = model, y = value, fill = model)) +
    geom_boxplot(outlier.shape = NA) +
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
      title = paste("Coefficient Estimates and Standard Error for Each Model (10 Runs, outliers removed) -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"),
      y = "Value",
      x = "Model"
    ) +
    theme(legend.position = "none")
  print(coefplot)
  ggsave(plot=coefplot,file=paste("Charts/Coef_x1x2_only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 9, height = 8, dpi = 900)
}

times_summary_avg <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_avg)=rownames(times_summary_all[[1]])
rownames(times_summary_avg) <- names(times_summary_all)

times_summary_sd <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_sd)=rownames(times_summary_all[[1]])
rownames(times_summary_sd) <- names(times_summary_all)

for (i in 1:length(times_summary_all)) {
  times_summary_avg[i,]=times_summary_all[[i]][, "mean"]
  times_summary_sd[i,]=times_summary_all[[i]][, "sd"]
}

print(times_summary_avg)
