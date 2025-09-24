# Evaluation Metric Plots - VS2 Weighted, AIC, BIC, -2*Log Likelihood Analysis
# This script creates evaluation metric plots with four subplots (one per distribution)
# following the same legend style and colors as the coefficient plots

cat("=== Evaluation Metric Plotting Script ===\n")
cat("Loading required data from main plotting script...\n")

# Check if required data objects exist (from main plotting script)
required_objects <- c("vs2_matrix", "vs2_wt_matrix", "loglik_matrix", "bic_matrix", "true_params_matrix", 
                      "skew_matrix", "model_colors", "rename_model", "models_to_plot")

missing_objects <- c()
for(obj in required_objects) {
  if(!exists(obj)) {
    missing_objects <- c(missing_objects, obj)
  }
}

if(length(missing_objects) > 0) {
  cat("ERROR: Missing required data objects:", paste(missing_objects, collapse = ", "), "\n")
  cat("Please run 'Plotting for Broad Models.r' first to generate the required data.\n")
  stop("Required data objects not found")
} else {
  cat("All required data objects found. Proceeding with evaluation metric plotting...\n")
}

# Load required libraries
library(ggplot2)
library(gridExtra)
library(grid)
library(latex2exp)

# OUTLIER REMOVAL FOR EVALUATION METRICS ----
cat("\n=== APPLYING OUTLIER REMOVAL FOR EVALUATION METRICS ===\n")

# Function to remove outliers using aggressive 1.5*IQR method for evaluation metrics
remove_outliers_eval <- function(matrix_data, matrix_name = "matrix", iqr_multiplier = 3, 
                                 remove_global_outliers = TRUE, verbose = TRUE) {
  
  if(verbose) cat("Processing", matrix_name, "for outlier removal...\n")
  
  # Make a copy to work with
  cleaned_matrix <- matrix_data
  total_removed <- 0
  
  # Step 1: Remove global outliers (computational failures)
  if(remove_global_outliers) {
    # Find infinite, NaN, or extremely large values
    # For evaluation metrics, use median absolute deviation approach
    global_outliers <- is.infinite(cleaned_matrix) | is.nan(cleaned_matrix)
    
    # Add extreme value detection based on median absolute deviation
    for(j in 1:ncol(cleaned_matrix)) {
      col_data <- cleaned_matrix[, j]
      if(all(is.na(col_data))) next
      
      median_val <- median(col_data, na.rm = TRUE)
      mad_val <- mad(col_data, na.rm = TRUE)
      
      # Only apply extreme value filter if MAD is not zero
      if(!is.na(mad_val) && mad_val > 0) {
        extreme_outliers <- !is.na(col_data) & abs(col_data - median_val) > 1000 * mad_val
        global_outliers[, j] <- global_outliers[, j] | extreme_outliers
      }
    }
    
    global_count <- sum(global_outliers, na.rm = TRUE)
    if(global_count > 0) {
      cleaned_matrix[global_outliers] <- NA
      total_removed <- total_removed + global_count
      if(verbose) cat("  Removed", global_count, "global outliers (inf/nan/extreme values)\n")
    }
  }
  
  # Step 2: Distribution-specific outlier removal using 1.5*IQR
  # Get unique distributions
  distributions <- unique(true_params_matrix[, "dist"])
  distributions <- distributions[!is.na(distributions)]
  
  for(dist in distributions) {
    if(verbose) cat("  Processing distribution:", dist, "\n")
    
    # Get rows for this distribution
    dist_rows <- which(true_params_matrix[, "dist"] == dist)
    
    if(length(dist_rows) == 0) next
    
    # Process each model column separately
    for(j in 1:ncol(cleaned_matrix)) {
      model_name <- colnames(cleaned_matrix)[j]
      
      # Get data for this distribution and model
      dist_model_data <- cleaned_matrix[dist_rows, j]
      
      # Skip if all missing
      if(all(is.na(dist_model_data))) next
      
      # Calculate IQR-based bounds
      q1 <- quantile(dist_model_data, 0.25, na.rm = TRUE)
      q3 <- quantile(dist_model_data, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      
      # Skip if IQR is 0 (no variation)
      if(iqr == 0 || is.na(iqr)) next
      
      # Define outlier bounds using 1.5*IQR (aggressive)
      lower_bound <- q1 - iqr_multiplier * iqr
      upper_bound <- q3 + iqr_multiplier * iqr
      
      # Identify outliers
      outliers <- !is.na(dist_model_data) & 
                 (dist_model_data < lower_bound | dist_model_data > upper_bound)
      
      outlier_count <- sum(outliers)
      
      if(outlier_count > 0) {
        # Remove outliers
        cleaned_matrix[dist_rows[outliers], j] <- NA
        total_removed <- total_removed + outlier_count
        
        if(verbose) {
          cat("    ", model_name, "- removed", outlier_count, "outliers",
              "(bounds: [", round(lower_bound, 3), ",", round(upper_bound, 3), "])\n")
        }
      }
    }
  }
  
  if(verbose) {
    remaining_values <- sum(!is.na(cleaned_matrix))
    original_values <- sum(!is.na(matrix_data))
    cat("  Total outliers removed:", total_removed, "\n")
    cat("  Remaining values:", remaining_values, "/", original_values, 
        "(", round(remaining_values/original_values*100, 1), "%)\n")
  }
  
  return(cleaned_matrix)
}

# Apply outlier removal to evaluation matrices
cat("Applying outlier removal to evaluation matrices...\n")

# Clean VS2 (unweighted) matrix
if(exists("vs2_matrix")) {
  vs2_matrix_clean <- remove_outliers_eval(vs2_matrix, "vs2_matrix", verbose = TRUE)
  cat("VS2 (unweighted) matrix cleaned.\n")
} else {
  cat("Warning: vs2_matrix not found, skipping outlier removal\n")
  vs2_matrix_clean <- NULL
}

# Clean VS2 weighted matrix
if(exists("vs2_wt_matrix")) {
  vs2_wt_matrix_clean <- remove_outliers_eval(vs2_wt_matrix, "vs2_wt_matrix", verbose = TRUE)
  cat("VS2 weighted matrix cleaned.\n")
} else {
  cat("Warning: vs2_wt_matrix not found, skipping outlier removal\n")
  vs2_wt_matrix_clean <- NULL
}

# Clean log likelihood matrix
if(exists("loglik_matrix")) {
  loglik_matrix_clean <- remove_outliers_eval(loglik_matrix, "loglik_matrix", verbose = TRUE)
  cat("Log likelihood matrix cleaned.\n")
} else {
  cat("Warning: loglik_matrix not found, skipping outlier removal\n")
  loglik_matrix_clean <- NULL
}

# Clean BIC matrix
if(exists("bic_matrix")) {
  bic_matrix_clean <- remove_outliers_eval(bic_matrix, "bic_matrix", verbose = TRUE)
  cat("BIC matrix cleaned.\n")
} else {
  cat("Warning: bic_matrix not found, skipping outlier removal\n")
  bic_matrix_clean <- NULL
}

# Clean AIC matrix if it exists
if(exists("aic_matrix") && !is.null(aic_matrix)) {
  aic_matrix_clean <- remove_outliers_eval(aic_matrix, "aic_matrix", verbose = TRUE)
  cat("AIC matrix cleaned.\n")
} else {
  cat("AIC matrix not found, skipping outlier removal\n")
  aic_matrix_clean <- NULL
}

cat("=== OUTLIER REMOVAL COMPLETED FOR EVALUATION METRICS ===\n")

# EVALUATION METRIC PLOTTING FUNCTIONS ----

# Define the evaluation score models (subset of all models that have evaluation scores)
# These are the models that have vs2, vs2_wt, loglik, and bic data
available_eval_models <- c("glm", "gee", "re_nosig", "lme4", "gamm", "cop", "cop_n", "cop_j", "cop_g", "cop_f")

# Filter to only models that are in both available_eval_models and models_to_plot
eval_models_to_plot <- intersect(models_to_plot, available_eval_models)

cat("Models available for evaluation plots:", paste(eval_models_to_plot, collapse = ", "), "\n")
cat("Model labels:", paste(sapply(eval_models_to_plot, rename_model), collapse = ", "), "\n")

# Function to prepare evaluation metric data for a specific distribution and metric
prepare_eval_data <- function(dist_name, metric_matrix, metric_name, models_to_use = NULL, use_cleaned_data = TRUE) {
  
  # If no models specified, use the evaluation models
  if(is.null(models_to_use)) {
    models_to_use <- eval_models_to_plot
  }
  
  # Use cleaned matrix if available and requested
  if(use_cleaned_data) {
    if(metric_name == "VS2" && exists("vs2_matrix_clean") && !is.null(vs2_matrix_clean)) {
      metric_matrix <- vs2_matrix_clean
      cat("Using cleaned VS2 (unweighted) data for", dist_name, "\n")
    } else if(metric_name == "VS2_Weighted" && exists("vs2_wt_matrix_clean") && !is.null(vs2_wt_matrix_clean)) {
      metric_matrix <- vs2_wt_matrix_clean
      cat("Using cleaned VS2 weighted data for", dist_name, "\n")
    } else if(metric_name == "NegLogLik" && exists("loglik_matrix_clean") && !is.null(loglik_matrix_clean)) {
      metric_matrix <- loglik_matrix_clean
      cat("Using cleaned log likelihood data for", dist_name, "\n")
    } else if(metric_name == "BIC" && exists("bic_matrix_clean") && !is.null(bic_matrix_clean)) {
      metric_matrix <- bic_matrix_clean
      cat("Using cleaned BIC data for", dist_name, "\n")
    } else if(metric_name == "AIC" && exists("aic_matrix_clean") && !is.null(aic_matrix_clean)) {
      metric_matrix <- aic_matrix_clean
      cat("Using cleaned AIC data for", dist_name, "\n")
    }
  }
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No", metric_name, "data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get metric values for this distribution and selected models
  metric_data <- metric_matrix[dist_rows, models_to_use, drop = FALSE]
  
  # Create long format data
  n_models <- length(models_to_use)
  n_obs <- length(dist_rows)
  
  plot_data <- data.frame(
    corr_vals = rep(corr_vals, n_models),
    metric_value = as.vector(metric_data),
    model = rep(models_to_use, each = n_obs),
    model_label = rep(sapply(models_to_use, rename_model), each = n_obs),
    dist = dist_name
  )
  
  # Remove rows with missing metric values
  plot_data <- plot_data[complete.cases(plot_data$metric_value), ]
  
  # Debug output
  cat(metric_name, "plot data for", dist_name, ":\n")
  cat("  Total observations:", nrow(plot_data), "\n")
  cat("  Models:", paste(unique(plot_data$model), collapse = ", "), "\n")
  cat("  Correlation range:", round(range(plot_data$corr_vals, na.rm = TRUE), 3), "\n")
  cat("  Metric range:", round(range(plot_data$metric_value, na.rm = TRUE), 3), "\n")
  
  return(plot_data)
}

# Function to create individual distribution plot for a metric
create_eval_dist_plot <- function(dist_name, metric_matrix, metric_name, y_label, 
                                  models_to_use = NULL, show_legend = TRUE, 
                                  y_limits = NULL, lower_is_better = TRUE) {
  
  plot_data <- prepare_eval_data(dist_name, metric_matrix, metric_name, models_to_use)
  
  if(is.null(plot_data) || nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(plot_data$model_label)
  cat("Unique models in", dist_name, metric_name, "data:", paste(unique_models, collapse = ", "), "\n")
  
  # Create the plot with consistent color mapping
  p <- ggplot(plot_data, aes(x = corr_vals, y = metric_value, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Use consistent color mapping, only showing models present in this plot
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    coord_cartesian(xlim = c(.2, .7)) +
    labs(
      title = if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name,
      x = "Kendall's τ",
      y = y_label,
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = if(show_legend) "bottom" else "none"
    )
  
  # Apply y-limits if specified
  if(!is.null(y_limits)) {
    p <- p + coord_cartesian(xlim = c(.2, .7), ylim = y_limits)
  }
  
  if(DEBUG_LEGENDS && show_legend) {
    cat("DEBUG: create_eval_dist_plot for", dist_name, metric_name, "- plot created with legend position:", if(show_legend) "bottom" else "none", "\n")
    cat("DEBUG: create_eval_dist_plot for", dist_name, metric_name, "- unique models in plot:", paste(unique_models, collapse = ", "), "\n")
    if(exists("model_colors") && length(model_colors) > 0) {
      available_colors <- intersect(unique_models, names(model_colors))
      cat("DEBUG: create_eval_dist_plot for", dist_name, metric_name, "- available colors:", length(available_colors), "of", length(unique_models), "\n")
    } else {
      cat("DEBUG: create_eval_dist_plot for", dist_name, metric_name, "- model_colors not available!\n")
    }
  }
  
  return(p)
}

# Function to create a 4-panel plot for a specific metric
create_eval_metric_plot <- function(metric_matrix, metric_name, main_title, y_label, 
                                   models_to_use = NULL, y_limits = NULL, lower_is_better = TRUE) {
  
  cat("\n=== Creating", metric_name, "plots ===\n")
  
  # Create individual plots for each distribution in the specified order
  distributions <- c("NO", "PO", "GA", "LO")  # Normal, Negative Binomial, Gamma, Bernoulli
  individual_plots <- list()
  
  for(i in seq_along(distributions)) {
    dist <- distributions[i]
    cat("Processing", dist, "distribution for", metric_name, "...\n")
    
    # Only show legend on the first plot for legend extraction
    show_legend_this_plot <- (i == 1)
    
    plot <- create_eval_dist_plot(
      dist_name = dist,
      metric_matrix = metric_matrix,
      metric_name = metric_name,
      y_label = y_label,
      models_to_use = models_to_use,
      show_legend = show_legend_this_plot,
      y_limits = y_limits,
      lower_is_better = lower_is_better
    )
    
    if(!is.null(plot)) {
      individual_plots[[length(individual_plots) + 1]] <- plot
    }
  }
  
  if(length(individual_plots) == 0) {
    cat("No valid", metric_name, "plots could be created\n")
    return(NULL)
  }
  
  # Extract legend using the same approach as coefficient plots
  cat("Extracting legend for", metric_name, "plots...\n")
  
  # Find the first distribution that has data for legend extraction
  legend_dist <- NULL
  for(dist in distributions) {  # Use the same order: Normal, Negative Binomial, Gamma, Bernoulli
    test_plot <- create_eval_dist_plot(
      dist_name = dist,
      metric_matrix = metric_matrix,
      metric_name = metric_name,
      y_label = y_label,
      models_to_use = models_to_use,
      show_legend = TRUE,
      y_limits = y_limits,
      lower_is_better = lower_is_better
    )
    
    if(!is.null(test_plot)) {
      legend_dist <- dist
      break
    }
  }
  
  if(is.null(legend_dist)) {
    cat("Error: No valid distribution found for", metric_name, "legend extraction\n")
    # Force manual legend creation (same as coefficient plots approach)
    shared_legend <- validate_legend(NULL, paste("Chart", metric_name), force_manual = TRUE)
  } else {
    cat("Using", legend_dist, "distribution for", metric_name, "legend extraction\n")
    
    # Extract legend from the working distribution with improved error handling (same as coefficient plots)
    tryCatch({
      legend_plot_temp <- create_eval_dist_plot(
        dist_name = legend_dist,
        metric_matrix = metric_matrix,
        metric_name = metric_name,
        y_label = y_label,
        models_to_use = models_to_use,
        show_legend = TRUE,
        y_limits = y_limits,
        lower_is_better = lower_is_better
      ) + 
        theme(
          legend.position = "bottom", 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 11),
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10)
        )
      
      # Extract the legend with improved method (same as coefficient plots)
      shared_legend <- get_legend(legend_plot_temp, force_manual = FALSE)
      
      # Validate and potentially replace with manual legend (same as coefficient plots)
      shared_legend <- validate_legend(shared_legend, paste("Chart", metric_name), force_manual = FALSE)
      
      cat("Legend extraction completed for", metric_name, "plots\n")
      
    }, error = function(e) {
      cat("Error extracting", metric_name, "legend:", e$message, "\n")
      shared_legend <- validate_legend(NULL, paste("Chart", metric_name), force_manual = TRUE)
    })
  }
  
  # Remove legends from individual plots (same as coefficient plots)
  plots_no_legend <- lapply(individual_plots, function(p) p + theme(legend.position = "none"))
  
  # Arrange plots in a single row (1 row x 4 columns)
  plots_grid <- do.call(arrangeGrob, c(plots_no_legend, ncol = 4, nrow = 1))
  
  # Combine the plots grid with the shared legend at the bottom (same as coefficient plots)
  combined_plot <- arrangeGrob(
    plots_grid, 
    shared_legend, 
    ncol = 1, 
    heights = c(10, 1),
    top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  return(combined_plot)
}

# CREATE THE EVALUATION METRIC PLOTS ----

cat("\n=== CREATING EVALUATION METRIC PLOTS ===\n")

# 1. VS2 (Unweighted) Plot
cat("\n--- VS2 (Unweighted) Analysis ---\n")
if(exists("vs2_matrix_clean") && !is.null(vs2_matrix_clean)) {
  vs2_plot <- create_eval_metric_plot(
    metric_matrix = vs2_matrix_clean,
    metric_name = "VS2",
    main_title = "Variogram Score (VS2) by Distribution",
    y_label = "VS2 Score",
    models_to_use = eval_models_to_plot,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
} else {
  cat("Warning: Cleaned VS2 (unweighted) matrix not available, using original data\n")
  if(exists("vs2_matrix")) {
    vs2_plot <- create_eval_metric_plot(
      metric_matrix = vs2_matrix,
      metric_name = "VS2",
      main_title = "Variogram Score (VS2) by Distribution",
      y_label = "VS2 Score",
      models_to_use = eval_models_to_plot,
      y_limits = NULL,  # Auto-scale
      lower_is_better = TRUE
    )
  } else {
    cat("Error: VS2 matrix not found\n")
    vs2_plot <- NULL
  }
}

if(!is.null(vs2_plot)) {
  # Display the plot
  grid.draw(vs2_plot)
  
  # Save the plot
  if(!dir.exists("Charts")) {
    dir.create("Charts")
  }
  
  chart_filename_vs2 <- paste0("Charts/VS2_by_Distribution_", Sys.Date(), ".png")
  # Save the plot with adjusted dimensions for single row layout
  ggsave(chart_filename_vs2, plot = vs2_plot, width = 20, height = 6, dpi = 600)
  cat("VS2 (unweighted) plot saved to:", chart_filename_vs2, "\n")
} else {
  cat("Could not create VS2 (unweighted) plot\n")
}

# 2. VS2 Weighted Plot
cat("\n--- VS2 Weighted Analysis ---\n")
if(exists("vs2_wt_matrix_clean") && !is.null(vs2_wt_matrix_clean)) {
  vs2_wt_plot <- create_eval_metric_plot(
    metric_matrix = vs2_wt_matrix_clean,
    metric_name = "VS2_Weighted",
    main_title = "Weighted Variogram Score (VS2) by Distribution",
    y_label = "VS2 Weighted Score",
    models_to_use = eval_models_to_plot,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
} else {
  cat("Warning: Cleaned VS2 weighted matrix not available, using original data\n")
  vs2_wt_plot <- create_eval_metric_plot(
    metric_matrix = vs2_wt_matrix,
    metric_name = "VS2_Weighted",
    main_title = "Weighted Variogram Score (VS2) by Distribution",
    y_label = "VS2 Weighted Score",
    models_to_use = eval_models_to_plot,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
}

if(!is.null(vs2_wt_plot)) {
  # Display the plot
  grid.draw(vs2_wt_plot)
  
  # Save the plot
  if(!dir.exists("Charts")) {
    dir.create("Charts")
  }
  
  chart_filename_vs2_wt <- paste0("Charts/VS2_Weighted_by_Distribution_", Sys.Date(), ".png")
  # Save the plot with adjusted dimensions for single row layout
  ggsave(chart_filename_vs2_wt, plot = vs2_wt_plot, width = 20, height = 6, dpi = 600)
  cat("VS2 Weighted plot saved to:", chart_filename_vs2_wt, "\n")
} else {
  cat("Could not create VS2 Weighted plot\n")
}

# 3. -2*Log Likelihood Plot
cat("\n--- -2*Log Likelihood Analysis ---\n")
# Exclude GEE and GAMM models for log likelihood plots
loglik_models_to_use <- setdiff(eval_models_to_plot, c("gee", "gamm"))
cat("Models for log likelihood plots (excluding GEE and GAMM):", paste(loglik_models_to_use, collapse = ", "), "\n")

if(exists("loglik_matrix_clean") && !is.null(loglik_matrix_clean)) {
  loglik_plot <- create_eval_metric_plot(
    metric_matrix = loglik_matrix_clean,
    metric_name = "NegLogLik",
    main_title = "-2×Log Likelihood by Distribution",
    y_label = "-2×Log Likelihood",
    models_to_use = loglik_models_to_use,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
} else {
  cat("Warning: Cleaned log likelihood matrix not available, using original data\n")
  loglik_plot <- create_eval_metric_plot(
    metric_matrix = loglik_matrix,
    metric_name = "NegLogLik",
    main_title = "-2×Log Likelihood by Distribution",
    y_label = "-2×Log Likelihood",
    models_to_use = loglik_models_to_use,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
}

if(!is.null(loglik_plot)) {
  # Display the plot
  grid.draw(loglik_plot)
  
  # Save the plot
  chart_filename_loglik <- paste0("Charts/NegLogLikelihood_by_Distribution_", Sys.Date(), ".png")
  # Save the plot with adjusted dimensions for single row layout
  ggsave(chart_filename_loglik, plot = loglik_plot, width = 20, height = 6, dpi = 600)
  cat("-2*Log Likelihood plot saved to:", chart_filename_loglik, "\n")
} else {
  cat("Could not create -2*Log Likelihood plot\n")
}

# 4. BIC Plot
cat("\n--- BIC Analysis ---\n")
# Exclude GEE and GAMM models for BIC plots
bic_models_to_use <- setdiff(eval_models_to_plot, c("gee", "gamm"))
cat("Models for BIC plots (excluding GEE and GAMM):", paste(bic_models_to_use, collapse = ", "), "\n")

if(exists("bic_matrix_clean") && !is.null(bic_matrix_clean)) {
  bic_plot <- create_eval_metric_plot(
    metric_matrix = bic_matrix_clean,
    metric_name = "BIC",
    main_title = "Bayesian Information Criterion (BIC) by Distribution",
    y_label = "BIC",
    models_to_use = bic_models_to_use,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
} else {
  cat("Warning: Cleaned BIC matrix not available, using original data\n")
  bic_plot <- create_eval_metric_plot(
    metric_matrix = bic_matrix,
    metric_name = "BIC",
    main_title = "Bayesian Information Criterion (BIC) by Distribution",
    y_label = "BIC",
    models_to_use = bic_models_to_use,
    y_limits = NULL,  # Auto-scale
    lower_is_better = TRUE
  )
}

if(!is.null(bic_plot)) {
  # Display the plot
  grid.draw(bic_plot)
  
  # Save the plot
  chart_filename_bic <- paste0("Charts/BIC_by_Distribution_", Sys.Date(), ".png")
  # Save the plot with adjusted dimensions for single row layout
  ggsave(chart_filename_bic, plot = bic_plot, width = 20, height = 6, dpi = 600)
  cat("BIC plot saved to:", chart_filename_bic, "\n")
} else {
  cat("Could not create BIC plot\n")
}

# 5. AIC Plot (derived from BIC and log-likelihood if needed)
# Note: AIC = -2*loglik + 2*df, while BIC = -2*loglik + df*log(n)
# If we have both loglik and BIC, we can compute AIC
cat("\n--- AIC Analysis (Computing from available data) ---\n")

# Try to create AIC matrix if it doesn't exist
if(!exists("aic_matrix") || is.null(aic_matrix)) {
  cat("AIC matrix not found, attempting to compute from BIC and log-likelihood...\n")
  
  # We need df_values to compute AIC properly
  # For now, we'll skip AIC if we don't have the degrees of freedom information
  if(exists("df_values") && !is.null(df_values)) {
    # Compute AIC = -2*loglik + 2*df
    # This would require accessing df_values from the original data
    cat("Computing AIC from degrees of freedom data...\n")
    # This is complex without the original structure, so we'll skip it for now
    cat("AIC computation requires access to degrees of freedom data from original results.\n")
    cat("Skipping AIC plot for now.\n")
    aic_plot <- NULL
  } else {
    cat("Cannot compute AIC without degrees of freedom information.\n")
    cat("Skipping AIC plot.\n")
    aic_plot <- NULL
  }
  } else {
    # AIC matrix exists, use cleaned version if available
    if(exists("aic_matrix_clean") && !is.null(aic_matrix_clean)) {
      aic_plot <- create_eval_metric_plot(
        metric_matrix = aic_matrix_clean,
        metric_name = "AIC",
        main_title = "Akaike Information Criterion (AIC) by Distribution",
        y_label = "AIC",
        models_to_use = eval_models_to_plot,
        y_limits = NULL,  # Auto-scale
        lower_is_better = TRUE
      )
    } else {
      cat("Warning: Cleaned AIC matrix not available, using original data\n")
      aic_plot <- create_eval_metric_plot(
        metric_matrix = aic_matrix,
        metric_name = "AIC",
        main_title = "Akaike Information Criterion (AIC) by Distribution",
        y_label = "AIC",
        models_to_use = eval_models_to_plot,
        y_limits = NULL,  # Auto-scale
        lower_is_better = TRUE
      )
    }
    
    if(!is.null(aic_plot)) {
      # Display the plot
      grid.draw(aic_plot)
    
    # Save the plot
    chart_filename_aic <- paste0("Charts/AIC_by_Distribution_", Sys.Date(), ".png")
    # Save the plot with adjusted dimensions for single row layout
    ggsave(chart_filename_aic, plot = aic_plot, width = 20, height = 6, dpi = 600)
    cat("AIC plot saved to:", chart_filename_aic, "\n")
  } else {
    cat("Could not create AIC plot\n")
  }
}

# SUMMARY ----
cat("\n=== EVALUATION METRIC PLOTS SUMMARY ===\n")
cat("Plots created for the following metrics:\n")
if(!is.null(vs2_wt_plot)) cat("✓ VS2 Weighted\n") else cat("✗ VS2 Weighted\n")
if(!is.null(loglik_plot)) cat("✓ -2×Log Likelihood\n") else cat("✗ -2×Log Likelihood\n")
if(!is.null(bic_plot)) cat("✓ BIC\n") else cat("✗ BIC\n")
if(exists("aic_plot") && !is.null(aic_plot)) cat("✓ AIC\n") else cat("✗ AIC (skipped - requires df data)\n")

cat("\nAll evaluation metric plots use the same legend style and colors as coefficient plots.\n")
cat("Each plot contains 4 subplots in a single row: Normal, Negative Binomial, Gamma, and Bernoulli distributions.\n")
cat("Charts saved to the 'Charts' directory.\n")
cat("\n=== EVALUATION METRIC PLOTTING COMPLETED ===\n")