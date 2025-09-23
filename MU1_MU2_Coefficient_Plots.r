# Check if computed true SEs are available
if(exists("true_ses_matrix") && !is.null(true_ses_matrix)) {
  cat("Using newly computed true SEs from true_ses_matrix\n")
  cat("True SEs matrix dimensions:", dim(true_ses_matrix), "\n")
  cat("Available columns:", colnames(true_ses_matrix), "\n")
} else {
  cat("Using old pre-computed true SEs from true_SE_out\n")
  cat("Note: Run the true SE computation section at the end to use improved SEs\n")
}

# Create plots for estimation bias vs correlation for each distribution

models_to_plot <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")

# OUTLIER REMOVAL FUNCTION ----
cat("\n=== APPLYING OUTLIER REMOVAL FOR MU1/MU2 PLOTS ===\n")

# Function to remove outliers using conservative 3*IQR method
remove_outliers <- function(matrix_data, matrix_name = "matrix", iqr_multiplier = 3, 
                           remove_global_outliers = TRUE, verbose = TRUE) {
  
  if(verbose) cat("Processing", matrix_name, "for outlier removal...\n")
  
  # Make a copy to work with
  cleaned_matrix <- matrix_data
  total_removed <- 0
  
  # Step 1: Remove global outliers (computational failures)
  if(remove_global_outliers) {
    # Find infinite, NaN, or extremely large values
    global_outliers <- is.infinite(cleaned_matrix) | is.nan(cleaned_matrix) | 
                      abs(cleaned_matrix) > 1000 * median(abs(cleaned_matrix), na.rm = TRUE)
    
    global_count <- sum(global_outliers, na.rm = TRUE)
    if(global_count > 0) {
      cleaned_matrix[global_outliers] <- NA
      total_removed <- total_removed + global_count
      if(verbose) cat("  Removed", global_count, "global outliers (inf/nan/extreme values)\n")
    }
  }
  
  # Step 2: Distribution-specific outlier removal using 3*IQR
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
      
      # Define outlier bounds using 3*IQR (conservative)
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
    cat("  Total outliers removed from", matrix_name, ":", total_removed, "\n")
    cat("  Percentage removed:", round(100 * total_removed / length(matrix_data), 2), "%\n")
  }
  
  return(cleaned_matrix)
}

# Apply outlier removal to the coefficient and SE matrices
cat("\nApplying outlier removal to coefficient matrices...\n")
mu1_coef_matrix_clean <- remove_outliers(mu1_coef_matrix, "mu1_coef_matrix")
mu2_coef_matrix_clean <- remove_outliers(mu2_coef_matrix, "mu2_coef_matrix")

cat("\nApplying outlier removal to standard error matrices...\n")  
mu1_se_matrix_clean <- remove_outliers(mu1_se_matrix, "mu1_se_matrix")
mu2_se_matrix_clean <- remove_outliers(mu2_se_matrix, "mu2_se_matrix")

# Summary of outlier removal
cat("\n=== OUTLIER REMOVAL SUMMARY ===\n")

matrices_info <- list(
  "mu1_coef" = list(original = mu1_coef_matrix, cleaned = mu1_coef_matrix_clean),
  "mu2_coef" = list(original = mu2_coef_matrix, cleaned = mu2_coef_matrix_clean),
  "mu1_se" = list(original = mu1_se_matrix, cleaned = mu1_se_matrix_clean),
  "mu2_se" = list(original = mu2_se_matrix, cleaned = mu2_se_matrix_clean)
)

for(matrix_name in names(matrices_info)) {
  orig <- matrices_info[[matrix_name]]$original
  clean <- matrices_info[[matrix_name]]$cleaned
  
  orig_missing <- sum(is.na(orig))
  clean_missing <- sum(is.na(clean))
  additional_missing <- clean_missing - orig_missing
  
  cat(matrix_name, ": originally", orig_missing, "missing,", 
      additional_missing, "additional removed,", 
      clean_missing, "total missing\n")
}

cat("\nOutlier removal complete! Using cleaned matrices for plotting.\n")

# Create a consistent color mapping for all models
all_model_labels <- sapply(models_to_plot, rename_model)
model_colors <- setNames(scales::hue_pal()(length(all_model_labels)), all_model_labels)

cat("Model color mapping:\n")
for(i in 1:length(model_colors)) {
  cat(sprintf("  %s: %s\n", names(model_colors)[i], model_colors[i]))
}

# First, we need to prepare the data by matching parameter sets between 
# true_params_matrix, mu1_coef_matrix, and skew_matrix using row indices

# Initialize exclusion tracking
exclusion_summary_global <- data.frame(
  distribution = character(),
  model = character(),
  model_label = character(),
  total_values = integer(),
  excluded_values = integer(),
  exclusion_rate = numeric(),
  stringsAsFactors = FALSE
)

# MU1 functions ####
# Create a function to prepare plot data for a specific distribution
prepare_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu1_coef_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get true mu1 values for this distribution

  # Make link function adjustments. Apply log for PO and GA and logit for LO
  if(dist_name %in% c("GA")) {
    true_mu1 <- (as.numeric(true_params_matrix$mu1[dist_rows])*as.numeric(true_params_matrix$a[dist_rows]))
  } else if (dist_name == "PO") {
    true_mu1 <- (as.numeric(true_params_matrix$mu1[dist_rows])*as.numeric(true_params_matrix$b[dist_rows])*as.numeric(true_params_matrix$c[dist_rows]))
  } else if(dist_name == "LO") {
    true_mu1 <- as.numeric(true_params_matrix$mu1[dist_rows])
  } else {
    true_mu1 <- as.numeric(true_params_matrix$mu1[dist_rows])
  }

  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated mu1 values for this distribution and selected models
  if(dist_name %in% c("GA")) {
    est_mu1_matrix <- exp(mu1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])
  } else if (dist_name == "PO") {
    est_mu1_matrix <- exp(mu1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])
  } else if(dist_name == "LO") {
    est_mu1_matrix <- exp(mu1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])/(1+exp(mu1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE]))
  } else {
    est_mu1_matrix <- mu1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  }

  # Apply exclusion logic only for LO distribution
  est_mu1_matrix_filtered <- est_mu1_matrix
  exclusion_mask <- matrix(FALSE, nrow = nrow(est_mu1_matrix), ncol = ncol(est_mu1_matrix))
  
  if(dist_name == "LO") {
    # Create exclusion mask for values more than 10x different from true values
    # Calculate exclusion statistics but do not actually apply the exclusion
    for(i in 1:length(true_mu1)) {
      if(!is.na(true_mu1[i])) {
        for(j in 1:ncol(est_mu1_matrix)) {
          if(!is.na(est_mu1_matrix[i, j])) {
            ratio <- abs(est_mu1_matrix[i, j] / true_mu1[i])
            #if(ratio > 10 || ratio < 0.1) { # Using 10x and 0.1x threshold
              exclusion_mask[i, j] <- FALSE
            #}
          }
        }
      }
    }
    
    # Do NOT apply exclusion mask - keep all data for plotting
    # est_mu1_matrix_filtered[exclusion_mask] <- NA
    
    # Store exclusion info for summary (global variable) - only for LO
    if(!exists("exclusion_summary_global")) {
      exclusion_summary_global <<- data.frame(
        distribution = character(),
        model = character(),
        model_label = character(),
        total_values = integer(),
        excluded_values = integer(),
        exclusion_rate = numeric(),
        stringsAsFactors = FALSE
      )
    }
    
    # Calculate exclusion statistics for each model (only for LO)
    for(j in 1:length(models_to_plot)) {
      model <- models_to_plot[j]
      model_label <- rename_model(model)
      
      total_non_na <- sum(!is.na(est_mu1_matrix[, j]))
      excluded_count <- sum(exclusion_mask[, j], na.rm = TRUE)
      
      exclusion_info <- data.frame(
        distribution = dist_name,
        model = model,
        model_label = model_label,
        total_values = total_non_na,
        excluded_values = excluded_count,
        exclusion_rate = if(total_non_na > 0) excluded_count / total_non_na else 0,
        stringsAsFactors = FALSE
      )
      
      exclusion_summary_global <<- rbind(exclusion_summary_global, exclusion_info)
    }
  }
  
  # Calculate relative bias using filtered data (estimated - true)
  bias_matrix <- est_mu1_matrix_filtered / rep(true_mu1, each = 1) -1
  
  # Convert to long format for ggplot
  plot_data <- data.frame(
    corr_vals = rep(corr_vals, length(models_to_plot)),
    bias = as.vector(bias_matrix),
    model = rep(models_to_plot, each = length(true_mu1)),
    model_label = rep(sapply(models_to_plot, rename_model), each = length(true_mu1)),
    dist = dist_name
  )
  
  # Remove rows with missing values
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  return(plot_data)
}

# Create individual plots for each distribution
create_dist_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  plot_data <- prepare_plot_data(dist_name, models_to_plot)
  
  if(is.null(plot_data) || nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(plot_data$model_label)
  cat("Unique models in", dist_name, "data:", paste(unique_models, collapse = ", "), "\n")
  
  # Create the plot with consistent color mapping
  p <- ggplot(plot_data, aes(x = corr_vals, y = bias, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    # Use consistent color mapping, only showing models
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    coord_cartesian(ylim = c(-1,1),xlim=c(.2,.7)) +
    labs(
      title = if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name,
      x = "Kendall's τ",
      y = TeX("$(\\hat{\\mu_1}/\\mu_1)-1$"),
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
  
  if(DEBUG_LEGENDS && show_legend) {
    cat("DEBUG: create_dist_plot for", dist_name, "- plot created with legend position:", if(show_legend) "bottom" else "none", "\n")
    cat("DEBUG: create_dist_plot for", dist_name, "- unique models in plot:", paste(unique_models, collapse = ", "), "\n")
    if(exists("model_colors") && length(model_colors) > 0) {
      available_colors <- intersect(unique_models, names(model_colors))
      cat("DEBUG: create_dist_plot for", dist_name, "- available colors:", length(available_colors), "of", length(unique_models), "\n")
    } else {
      cat("DEBUG: create_dist_plot for", dist_name, "- model_colors not available!\n")
    }
  }
  
  return(p)
}

# Create a function to prepare SE plot data for a specific distribution
prepare_se_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu1_se_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No SE data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated SE values for this distribution and selected models
  est_se_matrix <- mu1_se_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  
  # NEW CODE: Get theoretical SE values from true_ses_matrix (computed at end of script)
  # Check if true_ses_matrix is available (computed by get_true_ses function)
  if(exists("true_ses_matrix") && !is.null(true_ses_matrix)) {
    # Debug: Check structure of true_ses_matrix
    cat("Debug - true_ses_matrix structure for", dist_name, ":\n")
    cat("  Class:", class(true_ses_matrix), "\n")
    
    # true_ses_matrix is actually a list - extract the actual matrix
    if(is.list(true_ses_matrix) && "true_ses_matrix" %in% names(true_ses_matrix)) {
      actual_ses_matrix <- true_ses_matrix$true_ses_matrix
      cat("  Found true_ses_matrix inside list, class:", class(actual_ses_matrix), "\n")
      if(is.matrix(actual_ses_matrix)) {
        cat("  Matrix dimensions:", dim(actual_ses_matrix), "\n")
        cat("  Column names:", colnames(actual_ses_matrix), "\n")
        cat("  Number of dist_rows:", length(dist_rows), "\n")
      }
    } else {
      actual_ses_matrix <- true_ses_matrix
      cat("  Using true_ses_matrix directly, class:", class(actual_ses_matrix), "\n")
      if(is.matrix(actual_ses_matrix)) {
        cat("  Matrix dimensions:", dim(actual_ses_matrix), "\n")
        cat("  Column names:", colnames(actual_ses_matrix), "\n")
        cat("  Number of dist_rows:", length(dist_rows), "\n")
      }
    }

    # Use t1_se from the actual matrix for mu1 theoretical SEs with proper indexing
    if(is.matrix(actual_ses_matrix) && "t1_se" %in% colnames(actual_ses_matrix)) {
      # Check if dist_rows are within bounds
      if(max(dist_rows) <= nrow(actual_ses_matrix)) {
        true_se_vals <- actual_ses_matrix[dist_rows, "t1_se"]
        cat("Using computed true SEs from true_ses_matrix for", dist_name, "distribution\n")
      } else {
        cat("Error: dist_rows out of bounds for", dist_name, "- max dist_row:", max(dist_rows), "matrix rows:", nrow(actual_ses_matrix), "\n")
        true_se_vals <- rep(NA, length(dist_rows))
      }
    } else if(is.data.frame(actual_ses_matrix) && "t1_se" %in% colnames(actual_ses_matrix)) {
      if(max(dist_rows) <= nrow(actual_ses_matrix)) {
        true_se_vals <- actual_ses_matrix$t1_se[dist_rows]
        cat("Using computed true SEs from data.frame for", dist_name, "distribution\n")
      } else {
        cat("Error: dist_rows out of bounds for", dist_name, "- max dist_row:", max(dist_rows), "data.frame rows:", nrow(actual_ses_matrix), "\n")
        true_se_vals <- rep(NA, length(dist_rows))
      }
    } else {
      cat("Warning: t1_se column not found in true_ses_matrix for", dist_name, "\n")
      if(is.matrix(actual_ses_matrix) || is.data.frame(actual_ses_matrix)) {
        cat("  Available columns:", colnames(actual_ses_matrix), "\n")
      }
      true_se_vals <- rep(NA, length(dist_rows))
    }
  } else {
    # Fallback to old method if true_ses_matrix not available
    cat("Warning: true_ses_matrix not available, using old true_SE_out method for", dist_name, "\n")
    true_se_vals <- rep(NA, length(dist_rows))
  }
  
  # Create long format data combining estimated and theoretical SEs
  n_models <- length(models_to_plot)
  n_obs <- length(dist_rows)
  
  # Estimated SEs
  estimated_data <- data.frame(
    corr_vals = rep(corr_vals, n_models),
    se_value = as.vector(est_se_matrix),
    model = rep(models_to_plot, each = n_obs),
    model_label = rep(sapply(models_to_plot, rename_model), each = n_obs),
    type = "Estimated",
    dist = dist_name
  )
  
  # Theoretical SEs - only add if we have non-NA values
  if(any(!is.na(true_se_vals))) {
    theoretical_data <- data.frame(
      corr_vals = corr_vals,
      se_value = true_se_vals,
      model = "theoretical",
      model_label = "Theoretical",
      type = "Theoretical",
      dist = dist_name
    )
    
    # Combine estimated and theoretical data
    combined_se_data <- rbind(estimated_data, theoretical_data)
  } else {
    # Only use estimated data if no theoretical SEs available
    combined_se_data <- estimated_data
    cat("No theoretical SEs available for", dist_name, ", using only estimated SEs\n")
  }
  
  # Remove rows with missing se_value
  combined_se_data <- combined_se_data[complete.cases(combined_se_data$se_value), ]
  
  return(combined_se_data)
}

# Create SE plots for each distribution
create_se_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  se_data <- prepare_se_plot_data(dist_name, models_to_plot)
  
  if(is.null(se_data) || nrow(se_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(se_data$model_label)
  cat("Unique models in", dist_name, "SE data:", paste(unique_models, collapse = ", "), "\n")
  
  # Get unique model labels from the data and estimated model labels
  all_models_in_data <- unique(se_data$model_label)
  estimated_models <- unique(se_data$model_label[se_data$type == "Estimated"])
  cat("Models in", dist_name, "SE data:", paste(all_models_in_data, collapse = ", "), "\n")
  
  # Create the plot with consistent color mapping
  p <- ggplot(se_data, aes(x = corr_vals, y = se_value, color = model_label, linetype = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Manual color scale - only include models present in this plot
    scale_color_manual(values = c("Theoretical" = "black", model_colors[estimated_models]), 
                       breaks = all_models_in_data,
                       limits = all_models_in_data) +
    # Manual linetype scale - set theoretical to dashed, others to solid
    scale_linetype_manual(values = c("Theoretical" = "dashed",
                                    setNames(rep("solid", length(estimated_models)),
                                            estimated_models)),
                          breaks = all_models_in_data,
                          limits = all_models_in_data) +
    coord_cartesian(ylim = c(0, 0.4),xlim=c(.2,.7)) +
    labs(
      title = paste(if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name, "- SE"),
      x = "Kendall's τ",
      y = TeX("SE$(\\hat{\\mu_1})$"),
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
  
  return(p)
}

# MU2 PLOTTING FUNCTIONS ####

# Function to prepare mu2 plot data for a specific distribution
prepare_mu2_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu2_coef_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get true mu2 values for this distribution
  # Make link function adjustments. Apply log for PO and GA and logit for LO
  if(dist_name %in% c("GA")) {
    true_mu2 <- log(as.numeric(true_params_matrix$mu2[dist_rows])/as.numeric(true_params_matrix$mu1[dist_rows]))
  } else if (dist_name == "PO") {
    true_mu2 <- log(as.numeric(true_params_matrix$mu2[dist_rows])/as.numeric(true_params_matrix$mu1[dist_rows]))
  } else if(dist_name == "LO") {
    true_mu2 <- logit(as.numeric(true_params_matrix$mu2[dist_rows])) - logit(as.numeric(true_params_matrix$mu1[dist_rows]))
  } else {
    true_mu2 <- as.numeric(true_params_matrix$mu2[dist_rows])- as.numeric(true_params_matrix$mu1[dist_rows])
  }

  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated mu2 values for this distribution and selected models
  if(dist_name %in% c("GA")) {
    est_mu2_matrix <- (mu2_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])
  } else if (dist_name == "PO") {
    est_mu2_matrix <- (mu2_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])
  } else if(dist_name == "LO") {
    est_mu2_matrix <- (mu2_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE])
  } else {
    est_mu2_matrix <- mu2_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  }

  # Apply exclusion logic only for LO distribution
  est_mu2_matrix_filtered <- est_mu2_matrix
  exclusion_mask <- matrix(FALSE, nrow = nrow(est_mu2_matrix), ncol = ncol(est_mu2_matrix))
  
  if(dist_name == "LO") {
    # Create exclusion mask for values more than 10x different from true values
    for(i in 1:length(true_mu2)) {
      if(!is.na(true_mu2[i])) {
        for(j in 1:ncol(est_mu2_matrix)) {
          if(!is.na(est_mu2_matrix[i, j])) {
            ratio <- abs(est_mu2_matrix[i, j] / true_mu2[i])
            #if(ratio > 10 || ratio < 0.1) { # Using 10x and 0.1x threshold
              exclusion_mask[i, j] <- FALSE
            #}
          }
        }
      }
    }
    
    # Store exclusion info for summary (global variable) - only for LO
    if(!exists("exclusion_summary_mu2_global")) {
      exclusion_summary_mu2_global <<- data.frame(
        distribution = character(),
        model = character(),
        model_label = character(),
        total_values = integer(),
        excluded_values = integer(),
        exclusion_rate = numeric(),
        stringsAsFactors = FALSE
      )
    }
    
    # Calculate exclusion statistics for each model (only for LO)
    for(j in 1:length(models_to_plot)) {
      model_name <- models_to_plot[j]
      model_label <- rename_model(model_name)
      
      total_vals <- sum(!is.na(est_mu2_matrix[, j]))
      excluded_vals <- sum(exclusion_mask[, j], na.rm = TRUE)
      exclusion_rate <- if(total_vals > 0) excluded_vals / total_vals else 0
      
      exclusion_summary_mu2_global <<- rbind(exclusion_summary_mu2_global, data.frame(
        distribution = dist_name,
        model = model_name,
        model_label = model_label,
        total_values = total_vals,
        excluded_values = excluded_vals,
        exclusion_rate = exclusion_rate,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Calculate bias: (estimated - true) / true
  bias_matrix <- matrix(NA, nrow = nrow(est_mu2_matrix_filtered), ncol = ncol(est_mu2_matrix_filtered))
  
  for(i in 1:length(true_mu2)) {
    if(!is.na(true_mu2[i]) && true_mu2[i] != 0) {
      for(j in 1:ncol(est_mu2_matrix_filtered)) {
        if(!is.na(est_mu2_matrix_filtered[i, j])) {
          bias_matrix[i, j] <- (est_mu2_matrix_filtered[i, j] / true_mu2[i]) - 1
        }
      }
    }
  }
  
  # Create long format data
  n_models <- length(models_to_plot)
  n_obs <- nrow(bias_matrix)
  
  plot_data <- data.frame(
    corr_vals = rep(corr_vals, n_models),
    bias = as.vector(bias_matrix),
    model = rep(models_to_plot, each = n_obs),
    model_label = rep(sapply(models_to_plot, rename_model), each = n_obs),
    dist = dist_name
  )
  
  # Remove rows with missing bias values
  plot_data <- plot_data[complete.cases(plot_data$bias), ]
  
  # Print summary statistics
  cat("Mu2 plot data for", dist_name, ":\n")
  cat("  Total observations:", nrow(plot_data), "\n")
  cat("  Models:", paste(unique(plot_data$model), collapse = ", "), "\n")
  cat("  Correlation range:", round(range(plot_data$corr_vals, na.rm = TRUE), 3), "\n")
  cat("  Bias range:", round(range(plot_data$bias, na.rm = TRUE), 3), "\n")
  
  return(plot_data)
}

# Function to create mu2 bias plot for a specific distribution
create_mu2_dist_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  plot_data <- prepare_mu2_plot_data(dist_name, models_to_plot)
  
  if(is.null(plot_data) || nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(plot_data$model_label)
  cat("Unique models in", dist_name, "mu2 data:", paste(unique_models, collapse = ", "), "\n")
  
  # Create the plot with consistent color mapping
  p <- ggplot(plot_data, aes(x = corr_vals, y = bias, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    # Use consistent color mapping, only showing models present in this plot
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    coord_cartesian(ylim = c(-.25,.25),xlim=c(.2,.7)) +
    labs(
      title = if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name,
      x = "Kendall's τ",
      y = TeX("$(\\hat{\\mu_2}/\\mu_2)-1$"),
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
  
  return(p)
}

# Function to prepare mu2 SE plot data
prepare_mu2_se_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu2_se_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No mu2 SE data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated SE values for this distribution and selected models
  est_se_matrix <- mu2_se_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  
  # Get theoretical SE values from true_ses_matrix (computed at end of script)
  # Check if true_ses_matrix is available (computed by get_true_ses function)
  if(exists("true_ses_matrix") && !is.null(true_ses_matrix)) {
    # Debug: Check structure of true_ses_matrix
    cat("Debug - true_ses_matrix structure for mu2", dist_name, ":\n")
    cat("  Class:", class(true_ses_matrix), "\n")
    
    # true_ses_matrix is actually a list - extract the actual matrix
    if(is.list(true_ses_matrix) && "true_ses_matrix" %in% names(true_ses_matrix)) {
      actual_ses_matrix <- true_ses_matrix$true_ses_matrix
      cat("  Found true_ses_matrix inside list, class:", class(actual_ses_matrix), "\n")
      if(is.matrix(actual_ses_matrix)) {
        cat("  Matrix dimensions:", dim(actual_ses_matrix), "\n")
        cat("  Column names:", colnames(actual_ses_matrix), "\n")
        cat("  Number of dist_rows:", length(dist_rows), "\n")
      }
    } else {
      actual_ses_matrix <- true_ses_matrix
      cat("  Using true_ses_matrix directly, class:", class(actual_ses_matrix), "\n")
      if(is.matrix(actual_ses_matrix)) {
        cat("  Matrix dimensions:", dim(actual_ses_matrix), "\n")
        cat("  Column names:", colnames(actual_ses_matrix), "\n")
        cat("  Number of dist_rows:", length(dist_rows), "\n")
      }
    }

    # Use t2_se from the actual matrix for mu2 theoretical SEs with proper indexing
    if(is.matrix(actual_ses_matrix) && "t2_se" %in% colnames(actual_ses_matrix)) {
      # Check if dist_rows are within bounds
      if(max(dist_rows) <= nrow(actual_ses_matrix)) {
        true_se_vals <- actual_ses_matrix[dist_rows, "t2_se"]
        cat("Using computed true SEs from true_ses_matrix for mu2", dist_name, "distribution\n")
      } else {
        cat("Error: dist_rows out of bounds for mu2", dist_name, "- max dist_row:", max(dist_rows), "matrix rows:", nrow(actual_ses_matrix), "\n")
        true_se_vals <- rep(NA, length(dist_rows))
      }
    } else if(is.data.frame(actual_ses_matrix) && "t2_se" %in% colnames(actual_ses_matrix)) {
      if(max(dist_rows) <= nrow(actual_ses_matrix)) {
        true_se_vals <- actual_ses_matrix$t2_se[dist_rows]
        cat("Using computed true SEs from data.frame for mu2", dist_name, "distribution\n")
      } else {
        cat("Error: dist_rows out of bounds for mu2", dist_name, "- max dist_row:", max(dist_rows), "data.frame rows:", nrow(actual_ses_matrix), "\n")
        true_se_vals <- rep(NA, length(dist_rows))
      }
    } else {
      cat("Warning: t2_se column not found in true_ses_matrix for mu2", dist_name, "\n")
      if(is.matrix(actual_ses_matrix) || is.data.frame(actual_ses_matrix)) {
        cat("  Available columns:", colnames(actual_ses_matrix), "\n")
      }
      true_se_vals <- rep(NA, length(dist_rows))
    }
  } else {
    # Fallback to old method if true_ses_matrix not available
    cat("Warning: true_ses_matrix not available, using old true_SE_out method for mu2", dist_name, "\n")
    true_se_vals <- rep(NA, length(dist_rows))
  }
  
  # Create long format data combining estimated and theoretical SEs
  n_models <- length(models_to_plot)
  n_obs <- length(dist_rows)
  
  # Estimated SEs
  estimated_data <- data.frame(
    corr_vals = rep(corr_vals, n_models),
    se_value = as.vector(est_se_matrix),
    model = rep(models_to_plot, each = n_obs),
    model_label = rep(sapply(models_to_plot, rename_model), each = n_obs),
    type = "Estimated",
    dist = dist_name
  )
  
  # Theoretical SEs - only add if we have non-NA values
  if(any(!is.na(true_se_vals))) {
    theoretical_data <- data.frame(
      corr_vals = corr_vals,
      se_value = true_se_vals,
      model = "theoretical",
      model_label = "Theoretical",
      type = "Theoretical",
      dist = dist_name
    )
    
    # Combine estimated and theoretical data
    combined_se_data <- rbind(estimated_data, theoretical_data)
  } else {
    # Only use estimated data if no theoretical SEs available
    combined_se_data <- estimated_data
    cat("No theoretical SEs available for mu2", dist_name, ", using only estimated SEs\n")
  }
  
  # Remove rows with missing se_value
  combined_se_data <- combined_se_data[complete.cases(combined_se_data$se_value), ]
  
  return(combined_se_data)
}

# Function to create mu2 SE plot for a specific distribution
create_mu2_se_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  se_data <- prepare_mu2_se_plot_data(dist_name, models_to_plot)
  
  if(is.null(se_data) || nrow(se_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(se_data$model_label)
  cat("Unique models in", dist_name, "mu2 SE data:", paste(unique_models, collapse = ", "), "\n")
  
  # Get unique model labels from the data and estimated model labels
  all_models_in_data <- unique(se_data$model_label)
  estimated_models <- unique(se_data$model_label[se_data$type == "Estimated"])
  cat("Models in", dist_name, "mu2 SE data:", paste(all_models_in_data, collapse = ", "), "\n")
  
  # Create the plot with consistent color mapping
  p <- ggplot(se_data, aes(x = corr_vals, y = se_value, color = model_label, linetype = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Manual color scale - only include models present in this plot
    scale_color_manual(values = c("Theoretical" = "black", model_colors[estimated_models]), 
                       breaks = all_models_in_data,
                       limits = all_models_in_data) +
    # Manual linetype scale - set theoretical to dashed, others to solid
    scale_linetype_manual(values = c("Theoretical" = "dashed",
                                    setNames(rep("solid", length(estimated_models)),
                                            estimated_models)),
                          breaks = all_models_in_data,
                          limits = all_models_in_data) +
    coord_cartesian(ylim = c(0, 0.4),xlim=c(.2,.7)) +
    labs(
      title = paste(if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name, "- SE"),
      x = "Kendall's τ",
      y = TeX("SE$(\\hat{\\mu_1})$"),
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
  
  return(p)
}

# Create mu1 BIAS AND SE CHARTS VS CORRELATION ----

# Create plots for each distribution
cat("Creating mu1 bias vs correlation plots for each distribution...\n")

# Debug: Check true_ses_matrix structure before plotting
cat("\n=== DEBUG: Checking true_ses_matrix structure ===\n")
if(exists("true_ses_matrix")) {
  cat("true_ses_matrix exists\n")
  cat("Class:", class(true_ses_matrix), "\n")
  if(is.matrix(true_ses_matrix)) {
    cat("Dimensions:", dim(true_ses_matrix), "\n")
    cat("Column names:", colnames(true_ses_matrix), "\n")
    cat("Row names (first 5):", head(rownames(true_ses_matrix), 5), "\n")
  } else if(is.data.frame(true_ses_matrix)) {
    cat("Dimensions:", dim(true_ses_matrix), "\n")
    cat("Column names:", colnames(true_ses_matrix), "\n")
  } else if(is.list(true_ses_matrix)) {
    cat("List structure with names:", names(true_ses_matrix), "\n")
  } else {
    cat("Length:", length(true_ses_matrix), "\n")
    cat("Names:", names(true_ses_matrix), "\n")
  }
} else {
  cat("true_ses_matrix does not exist\n")
}
cat("=== END DEBUG ===\n\n")

# Create bias plots without individual legends (we'll add a shared legend)
plot_NO_bias <- create_dist_plot("NO", models_to_plot, show_legend = FALSE)
plot_PO_bias <- create_dist_plot("PO", models_to_plot, show_legend = FALSE) 
plot_GA_bias <- create_dist_plot("GA", models_to_plot, show_legend = FALSE)
plot_LO_bias <- create_dist_plot("LO", models_to_plot, show_legend = FALSE)

# Create SE plots without individual legends
plot_NO_se <- tryCatch(create_se_plot("NO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating NO SE plot:", e$message, "\n"); NULL })
plot_PO_se <- tryCatch(create_se_plot("PO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating PO SE plot:", e$message, "\n"); NULL })
plot_GA_se <- tryCatch(create_se_plot("GA", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating GA SE plot:", e$message, "\n"); NULL })
plot_LO_se <- tryCatch(create_se_plot("LO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating LO SE plot:", e$message, "\n"); NULL })

# Create a plot with legend to extract it for mu1
legend_plot <- create_dist_plot("NO", models_to_plot, show_legend = TRUE)

# Remove NULL plots and organize in pairs
bias_plots <- list(plot_NO_bias, plot_PO_bias, plot_GA_bias, plot_LO_bias)
se_plots <- list(plot_NO_se, plot_PO_se, plot_GA_se, plot_LO_se)

# Filter out NULL plots but keep the pairing
valid_indices <- which(!sapply(bias_plots, is.null) & !sapply(se_plots, is.null))
bias_plots <- bias_plots[valid_indices]
se_plots <- se_plots[valid_indices]

if(length(bias_plots) > 0 && length(se_plots) > 0) {
  # Use gridExtra for plot arrangement with shared legend
  library(gridExtra)
  library(grid)
  
  # First, check what distributions are available
  available_dists <- unique(true_params_matrix$dist)
  cat("Available distributions:", paste(available_dists, collapse = ", "), "\n")
  
  # Find the first distribution that has data for legend extraction
  legend_dist <- NULL
  for(dist in c("GA", "NO", "PO", "LO")) {
    if(dist %in% available_dists) {
      test_plot <- create_dist_plot(dist, models_to_plot, show_legend = TRUE)

      if(!is.null(test_plot)) {
        legend_dist <- dist
        break
      }
    }
  }
  
  if(is.null(legend_dist)) {
    cat("Error: No valid distribution found for legend extraction\n")
    # Force manual legend creation
    shared_legend <- validate_legend(NULL, "Chart 2 - Mu1", force_manual = TRUE)
  } else {
    cat("Using", legend_dist, "distribution for legend extraction\n")
    
    # Extract legend from the working distribution with improved error handling
    tryCatch({
      legend_plot_temp <- create_dist_plot(legend_dist, models_to_plot, show_legend = TRUE) + 
        theme(
          legend.position = "bottom", 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 11),
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10)
        )
      
      # Extract the legend with improved method
      shared_legend <- get_legend(legend_plot_temp, force_manual = FALSE)
      
      # Validate and potentially replace with manual legend
      shared_legend <- validate_legend(shared_legend, "Chart 2 - Mu1", force_manual = FALSE)
      
      cat("Legend extraction completed for Mu1 plots\n")
      
    }, error = function(e) {
      cat("Error extracting legend:", e$message, "\n")
      shared_legend <- validate_legend(NULL, "Chart 2 - Mu1", force_manual = TRUE)
    })
  }
  
  # Interleave bias and SE plots for each distribution
  all_plots <- vector("list", length(bias_plots) * 2)
  for(i in 1:length(bias_plots)) {
    all_plots[2*i-1] <- bias_plots[i]  # Bias plot (left column)
    all_plots[2*i] <- se_plots[i]      # SE plot (right column)
  }
  
  # Arrange plots in a 4x2 grid (4 rows, 2 columns)
  plots_grid <- do.call(arrangeGrob, c(all_plots, ncol = 2))
  
  # Combine the plots grid with the shared legend at the bottom
  combined_plot <- arrangeGrob(plots_grid, shared_legend, 
                              ncol = 1, 
                              heights = c(10, 1))
  
  # Display the combined plot
  grid.draw(combined_plot)
  
  # Save the combined plot to Charts directory
  chart_filename <- paste0("Charts/Mu1_Bias_and_SE_vs_Correlation_Combined_", Sys.Date(), ".png")
  
  ggsave(chart_filename, plot = combined_plot, width = 15, height = 14, dpi = 600)
  cat("Combined bias and SE plot saved to:", chart_filename, "\n")
  
  cat("Mu1 bias and SE vs correlation plots created successfully!\n")
  cat("Number of bias plots created:", length(bias_plots), "\n")
  cat("Number of SE plots created:", length(se_plots), "\n")
  
} else {
  cat("No plots could be created - check data availability\n")
}

# Print some diagnostic information
cat("\nDiagnostic information:\n")
cat("Distribution counts in true_params_matrix:\n")
print(table(true_params_matrix$dist))
cat("\nSkew matrix dimensions:", dim(skew_matrix), "\n")
cat("Available models in mu1_coef_matrix:\n")
print(colnames(mu1_coef_matrix))
cat("\nNon-missing values in mu1_coef_matrix by model:\n")
for(model in models_to_plot) {
  if(model %in% colnames(mu1_coef_matrix_clean)) {
    non_missing <- sum(!is.na(mu1_coef_matrix_clean[, model]))
    model_label <- rename_model(model)
    cat(model_label, "(", model, "):", non_missing, "non-missing values\n")
  } else {
    cat(model, ": not found in data\n")
  }
}

# Print exclusion summary (only for LO distribution)
if(exists("exclusion_summary_global") && nrow(exclusion_summary_global) > 0) {
  cat("\n=== EXCLUSION SUMMARY (LO Distribution Only) ===\n")
  cat("Values that would be excluded for being >10x or <0.1x the true value:\n")
  cat("(Note: Exclusions calculated but NOT applied - all data used in plots)\n\n")
  
  # Print by distribution and model
  for(dist in unique(exclusion_summary_global$distribution)) {
    cat("Distribution:", dist, "\n")
    dist_data <- exclusion_summary_global[exclusion_summary_global$distribution == dist, ]
    
    for(i in 1:nrow(dist_data)) {
      row <- dist_data[i, ]
      cat(sprintf("  %s (%s): %d/%d excluded (%.1f%%)\n", 
                  row$model_label, row$model, row$excluded_values, 
                  
                  row$total_values, row$exclusion_rate * 100))
    }
    cat("\n")
  }
  
  # Overall summary
  total_excluded <- sum(exclusion_summary_global$excluded_values)
  total_values <- sum(exclusion_summary_global$total_values)
  cat(sprintf("OVERALL: %d/%d values excluded (%.1f%%)\n", 
              total_excluded, total_values, 
              if(total_values > 0) total_excluded / total_values * 100 else 0))
} else {
  cat("\nNo exclusion data available (exclusions only applied to LO distribution).\n")
}

# Create mu2 bias and SE plots for each distribution
cat("Creating mu2 bias and SE plots for each distribution...\n")

# Debug: Check true_ses_matrix structure before mu2 plotting
cat("\n=== DEBUG: Checking true_ses_matrix structure for mu2 ===\n")
if(exists("true_ses_matrix")) {
  cat("true_ses_matrix exists\n")
  cat("Class:", class(true_ses_matrix), "\n")
  if(is.matrix(true_ses_matrix)) {
    cat("Dimensions:", dim(true_ses_matrix), "\n")
    cat("Column names:", colnames(true_ses_matrix), "\n")
    cat("Row names (first 5):", head(rownames(true_ses_matrix), 5), "\n")
  } else if(is.data.frame(true_ses_matrix)) {
    cat("Dimensions:", dim(true_ses_matrix), "\n")
    cat("Column names:", colnames(true_ses_matrix), "\n")
  } else if(is.list(true_ses_matrix)) {
    cat("List structure with names:", names(true_ses_matrix), "\n")
  } else {
    cat("Length:", length(true_ses_matrix), "\n")
    cat("Names:", names(true_ses_matrix), "\n")
  }
} else {
  cat("true_ses_matrix does not exist\n")
}
cat("=== END DEBUG for mu2 ===\n\n")

# Create mu2 bias plots without individual legends (we'll add a shared legend)
plot_NO_mu2_bias <- create_mu2_dist_plot("NO", models_to_plot, show_legend = FALSE)
plot_PO_mu2_bias <- create_mu2_dist_plot("PO", models_to_plot, show_legend = FALSE) 
plot_GA_mu2_bias <- create_mu2_dist_plot("GA", models_to_plot, show_legend = FALSE)
plot_LO_mu2_bias <- create_mu2_dist_plot("LO", models_to_plot, show_legend = FALSE)

# Create mu2 SE plots without individual legends
plot_NO_mu2_se <- tryCatch(create_mu2_se_plot("NO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating NO mu2 SE plot:", e$message, "\n"); NULL })
plot_PO_mu2_se <- tryCatch(create_mu2_se_plot("PO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating PO mu2 SE plot:", e$message, "\n"); NULL })
plot_GA_mu2_se <- tryCatch(create_mu2_se_plot("GA", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating GA mu2 SE plot:", e$message, "\n"); NULL })
plot_LO_mu2_se <- tryCatch(create_mu2_se_plot("LO", models_to_plot, show_legend = FALSE), error = function(e) { cat("Error creating LO mu2 SE plot:", e$message, "\n"); NULL })

# Create a plot with legend to extract it for mu2
legend_plot_mu2 <- create_mu2_dist_plot("NO", models_to_plot, show_legend = TRUE)

# Remove NULL plots and organize in pairs
bias_plots_mu2 <- list(plot_NO_mu2_bias, plot_PO_mu2_bias, plot_GA_mu2_bias, plot_LO_mu2_bias)
se_plots_mu2 <- list(plot_NO_mu2_se, plot_PO_mu2_se, plot_GA_mu2_se, plot_LO_mu2_se)

# Filter out NULL plots but keep the pairing
valid_indices_mu2 <- which(!sapply(bias_plots_mu2, is.null) & !sapply(se_plots_mu2, is.null))
bias_plots_mu2 <- bias_plots_mu2[valid_indices_mu2]
se_plots_mu2 <- se_plots_mu2[valid_indices_mu2]

if(length(bias_plots_mu2) > 0 && length(se_plots_mu2) > 0) {
  # Use gridExtra for plot arrangement with shared legend
  library(gridExtra)
  library(grid)
  
  # First, check what distributions are available for mu2
  available_dists_mu2 <- unique(true_params_matrix$dist)
  cat("Available distributions for mu2:", paste(available_dists_mu2, collapse = ", "), "\n")
  
  # Find the first distribution that has data for legend extraction
  legend_dist_mu2 <- NULL
  for(dist in c("GA", "NO", "PO", "LO")) {
    if(dist %in% available_dists_mu2) {
      test_plot_mu2 <- create_mu2_dist_plot(dist, models_to_plot, show_legend = TRUE)

      if(!is.null(test_plot_mu2)) {
        legend_dist_mu2 <- dist
        break
      }
    }
  }
  
  if(is.null(legend_dist_mu2)) {
    cat("Error: No valid distribution found for mu2 legend extraction\n")
    # Force manual legend creation
    shared_legend_mu2 <- validate_legend(NULL, "Chart 3 - Mu2", force_manual = TRUE)
  } else {
    cat("Using", legend_dist_mu2, "distribution for mu2 legend extraction\n")
    
    # Extract legend from the working distribution with improved error handling
    tryCatch({
      legend_plot_temp_mu2 <- create_mu2_dist_plot(legend_dist_mu2, models_to_plot, show_legend = TRUE) + 
        theme(
          legend.position = "bottom", 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 11),
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10)
        )
      
      # Extract the legend with improved method
      shared_legend_mu2 <- get_legend(legend_plot_temp_mu2, force_manual = FALSE)
      
      # Validate and potentially replace with manual legend
      shared_legend_mu2 <- validate_legend(shared_legend_mu2, "Chart 3 - Mu2", force_manual = FALSE)
      
      cat("Legend extraction completed for Mu2 plots\n")
      
    }, error = function(e) {
      cat("Error extracting mu2 legend:", e$message, "\n")
      shared_legend_mu2 <- validate_legend(NULL, "Chart 3 - Mu2", force_manual = TRUE)
    })
  }
  
  # Interleave bias and SE plots for each distribution
  all_plots_mu2 <- vector("list", length(bias_plots_mu2) * 2)
  for(i in 1:length(bias_plots_mu2)) {
    all_plots_mu2[2*i-1] <- bias_plots_mu2[i]  # Bias plot (left column)
    all_plots_mu2[2*i] <- se_plots_mu2[i]      # SE plot (right column)
  }
  
  # Arrange plots in a 4x2 grid (4 rows, 2 columns)
  plots_grid_mu2 <- do.call(arrangeGrob, c(all_plots_mu2, ncol = 2))
  
  # Combine the plots grid with the shared legend at the bottom
  combined_plot_mu2 <- arrangeGrob(plots_grid_mu2, shared_legend_mu2, 
                              ncol = 1, 
                              heights = c(10, 1))
  
  # Display the combined plot
  grid.draw(combined_plot_mu2)
  
  # Save the combined plot to Charts directory
  chart_filename_mu2 <- paste0("Charts/Mu2_Bias_and_SE_vs_Correlation_Combined_", Sys.Date(), ".png")
  
  ggsave(chart_filename_mu2, plot = combined_plot_mu2, width = 15, height = 14, dpi = 600)
  cat("Combined mu2 bias and SE plot saved to:", chart_filename_mu2, "\n")
  
  cat("Mu2 bias and SE vs correlation plots created successfully!\n")
  cat("Number of mu2 bias plots created:", length(bias_plots_mu2), "\n")
  cat("Number of mu2 SE plots created:", length(se_plots_mu2), "\n")
  
} else {
  cat("No mu2 plots could be created - check data availability\n")
}

# Print some diagnostic information for mu2
cat("\nDiagnostic information for mu2:\n")
cat("Distribution counts in true_params_matrix:\n")
print(table(true_params_matrix$dist))
cat("\nSkew matrix dimensions:", dim(skew_matrix), "\n")
cat("Available models in mu2_coef_matrix:\n")
print(colnames(mu2_coef_matrix))
cat("\nNon-missing values in mu2_coef_matrix by model:\n")
for(model in models_to_plot) {
  if(model %in% colnames(mu2_coef_matrix_clean)) {
    non_missing <- sum(!is.na(mu2_coef_matrix_clean[, model]))
    model_label <- rename_model(model)
    cat(model_label, "(", model, "):", non_missing, "non-missing values\n")
  } else {
    cat(model, ": not found in mu2 data\n")
  }
}

# Print exclusion summary for mu2 (only for LO distribution)
if(exists("exclusion_summary_mu2_global") && nrow(exclusion_summary_mu2_global) > 0) {
  cat("\n=== MU2 EXCLUSION SUMMARY (LO Distribution Only) ===\n")
  cat("Values that would be excluded for being >10x or <0.1x the true value:\n")
  cat("(Note: Exclusions calculated but NOT applied - all data used in plots)\n\n")
  
  # Print by distribution and model
  for(dist in unique(exclusion_summary_mu2_global$distribution)) {
    cat("Distribution:", dist, "\n")
    dist_data <- exclusion_summary_mu2_global[exclusion_summary_mu2_global$distribution == dist, ]
    
    for(i in 1:nrow(dist_data)) {
      row <- dist_data[i, ]
      cat(sprintf("  %s (%s): %d/%d excluded (%.1f%%)\n", 
                  row$model_label, row$model, row$excluded_values, 
                  
                  row$total_values, row$exclusion_rate * 100))
    }
    cat("\n")
  }
  
  # Overall summary for mu2
  total_excluded_mu2 <- sum(exclusion_summary_mu2_global$excluded_values)
  total_values_mu2 <- sum(exclusion_summary_mu2_global$total_values)
  cat(sprintf("OVERALL MU2: %d/%d values excluded (%.1f%%)\n", 
              total_excluded_mu2, total_values_mu2, 
              if(total_values_mu2 > 0) total_excluded_mu2 / total_values_mu2 * 100 else 0))
} else {
  cat("\nNo mu2 exclusion data available (exclusions only applied to LO distribution).\n")
}

# End of script - summarize legend debugging if enabled
