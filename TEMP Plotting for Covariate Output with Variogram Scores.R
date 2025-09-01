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

#selection_method <- readline("Enter choice (1, 2, or 3): ")
selection_method=2

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
      rdata= "Data/CoefSimData_PO_NA_0.2_5_5_5_1_0.01_100_10_2025-08-22.RData",
      meta=list(dist="PO", a=NA, b=.2, c=5, mu1=5, mu2=5, x1=1, x2=0.01, n=100)
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
  bic = "BIC",
  logliks = "-2*Log Likelihood"
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
    times_summary[i, "conv"] <- mean(conv[,i])
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
  
  # Calculate BIC using actual df values from df_values
  # BIC = -2*loglik + log(n)*df
  
  # Store original logliks matrix before processing
  original_logliks <- score_items$logliks
  
  # Calculate -2*loglik (standard transformation)
  score_items$logliks <- -2*score_items$logliks[, ref_models, drop = FALSE]
  
  # Calculate BIC using actual df values
  if (exists("df_values") && !is.null(df_values[["logliks"]])) {
    # Use actual df values from the data
    df_matrix <- df_values[["logliks"]][, ref_models, drop = FALSE]
    bic_matrix <- score_items$logliks + log(n) * df_matrix
  } else {
    # Fallback to default df if df_values not available
    cat("Warning: df_values not found, using default df = 4\n")
    bic_matrix <- score_items$logliks + log(n) * 4
  }
  
  score_items$bic <- bic_matrix
  
  # Set GEE log likelihood and BIC to NA (GEE doesn't have proper likelihood)
  if ("gee" %in% colnames(score_items$logliks)) {
    score_items$logliks[, "gee"] <- NA
  }
  if ("gee" %in% colnames(score_items$bic)) {
    score_items$bic[, "gee"] <- NA
  }
  
  # Remove vs2_wt_coronly 
  score_items$vs2_wt_coronly <- NULL

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

  score_type_order <- score_type_map[c("es", "vs1", "vs2", "vs2_wt", "bic", "logliks")]
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

  p=ggplot(score_long, aes(x = method_label, y = value, fill = method_label)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(data = meds, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~score_type_label, scales = "free", ncol = 2, labeller = labeller(score_type_label = label_value)) +
    theme_bw() +
    labs(title = paste("Model evaluation -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"),
         y = "Score Value") +
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())
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

  p2 = ggplot(score_no_outliers, aes(x = method_label, y = value, fill = method_label)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    geom_hline(data = meds_bottom4, aes(yintercept = min_median), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~score_type_label, scales = "free", ncol = 4, labeller = labeller(score_type_label = label_value)) +
    theme_bw() +
    labs(
      #title = paste("Model evaluation (Bottom Four Score Types, outliers removed) -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"),
      y = "Score Value"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x=element_blank()
    )

  print(p2)
  ggsave(paste("Charts/Variogram_p2only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),"_bottom4_no_outliers.png",sep=""), plot = p2, width = 10, height = 3, dpi = 900)

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
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank()) +
    labs(title = paste("Comparison of Models for x1, x2 -", dist, "dist (a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "x1:", x1, "x2:", x2, "n:", n, ")"), y = "Value")

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

  #dist_name=if(dist=="PO") {"Negative Binomial"} else if (dist=="NO") {"Normal"} else if (dist=="GA") {"Gamma"} else {"Bernoulli"}
  coefplot=ggplot(results_long2_no_outliers, aes(x = model, y = value, fill = model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ variable, scales = "free", ncol = 4,
               labeller = as_labeller(facet_labels)) +
    geom_hline(
      data = true_vals_df,
      aes(yintercept = true_value),
      colour = "red", linetype = "dashed", linewidth = 1
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank()) +
    labs(
      #title = paste(dist_name, ": X1, X2 Coef + SE (100 Simulations) -", "(a:", a, "b:", b, "c:", c, "μ1:", mu1, "μ2:", mu2, "n:", n, ")"),
      y = "Value"
    ) +
    theme(legend.position = "none")
  print(coefplot)
  ggsave(plot=coefplot,file=paste("Charts/Coef_x1x2_only_",paste(dist,a,b,c,mu1,mu2,x1,x2,n,Sys.Date(),sep="_"),".png",sep=""), width = 10, height = 3, dpi = 900)
}

# ---- Initialize summary matrices for mean, sd, and convergence ----
times_summary_avg <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_avg) = rownames(times_summary_all[[1]])
rownames(times_summary_avg) <- names(times_summary_all)

times_summary_sd <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_sd) = rownames(times_summary_all[[1]])
rownames(times_summary_sd) <- names(times_summary_all)

times_summary_conv <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_conv) = rownames(times_summary_all[[1]])
rownames(times_summary_conv) <- names(times_summary_all)

times_summary_conv <- matrix(nrow = length(times_summary_all), ncol = nrow(times_summary_all[[1]]))
colnames(times_summary_conv)=rownames(times_summary_all[[1]])
rownames(times_summary_conv) <- names(times_summary_all)

for (i in 1:length(times_summary_all)) {
  times_summary_avg[i,]=times_summary_all[[i]][, "mean"]
  times_summary_sd[i,]=times_summary_all[[i]][, "sd"]
  times_summary_conv[i,]=times_summary_all[[i]][, "conv"]
}

print(times_summary_avg)
print(times_summary_sd)
print(times_summary_conv)

# ---- Save summary tables to CSV files with date stamps ----
# Create a timestamp for the files
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
date_stamp <- Sys.Date()

# Create Charts directory if it doesn't exist
if (!dir.exists("Charts")) {
  dir.create("Charts", recursive = TRUE)
}

# Save timing summary tables
cat("Saving timing and convergence summary tables...\n")

# Save average timing data
avg_filename <- paste0("Charts/timing_summary_avg_", date_stamp, ".csv")
write.csv(times_summary_avg, file = avg_filename, row.names = TRUE)
cat("Saved average timing summary to:", avg_filename, "\n")

# Save standard deviation timing data
sd_filename <- paste0("Charts/timing_summary_sd_", date_stamp, ".csv")
write.csv(times_summary_sd, file = sd_filename, row.names = TRUE)
cat("Saved timing standard deviation summary to:", sd_filename, "\n")

# Save convergence data
conv_filename <- paste0("Charts/convergence_summary_", date_stamp, ".csv")
write.csv(times_summary_conv, file = conv_filename, row.names = TRUE)
cat("Saved convergence summary to:", conv_filename, "\n")

# Also save a combined summary with detailed timestamp
combined_filename <- paste0("Charts/timing_convergence_combined_", timestamp, ".csv")

# Create a combined data frame with clear column names
combined_data <- data.frame(
  Dataset = rownames(times_summary_avg),
  stringsAsFactors = FALSE
)

# Add timing averages with clear column names
for (method in colnames(times_summary_avg)) {
  combined_data[[paste0(method, "_avg_time")]] <- times_summary_avg[, method]
  combined_data[[paste0(method, "_sd_time")]] <- times_summary_sd[, method]
  combined_data[[paste0(method, "_convergence")]] <- times_summary_conv[, method]
}

write.csv(combined_data, file = combined_filename, row.names = FALSE)
cat("Saved combined timing and convergence summary to:", combined_filename, "\n")

cat("\n=== Summary of exported files ===\n")
cat("1. Average timing:", avg_filename, "\n")
cat("2. Timing std dev:", sd_filename, "\n") 
cat("3. Convergence rates:", conv_filename, "\n")
cat("4. Combined data:", combined_filename, "\n")
