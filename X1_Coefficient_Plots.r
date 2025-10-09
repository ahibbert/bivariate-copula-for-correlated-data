# X1 Coefficient Plots - Bias and SE Analysis
# This script creates bias and standard error plots for x1 coefficients
# following the same structure as the mu1 plots in "Plotting for Broad Models.r"

cat("=== X1 Coefficient Plotting Script ===\n")
cat("Loading required data from main plotting script...\n")

# Check if required data objects exist (from main plotting script)
required_objects <- c("x1_coef_matrix", "x1_se_matrix", "true_params_matrix", 
                      "skew_matrix", "model_colors", "models_to_plot", "true_ses_matrix")

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
  cat("All required data objects found. Proceeding with x1 plotting...\n")
}

# OUTLIER REMOVAL FOR X1 PLOTS ----
cat("\n=== APPLYING OUTLIER REMOVAL FOR X1 PLOTS ===\n")

# Function to remove outliers using conservative 3*IQR method
remove_outliers_x1 <- function(matrix_data, matrix_name = "matrix", iqr_multiplier = 3, 
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

# Apply outlier removal to X1 coefficient and SE matrices
cat("\nApplying outlier removal to X1 coefficient and SE matrices...\n")
x1_coef_matrix_clean <- remove_outliers_x1(x1_coef_matrix, "x1_coef_matrix")
x1_se_matrix_clean <- remove_outliers_x1(x1_se_matrix, "x1_se_matrix")

# Summary of outlier removal for X1
cat("\n=== X1 OUTLIER REMOVAL SUMMARY ===\n")

matrices_info_x1 <- list(
  "x1_coef" = list(original = x1_coef_matrix, cleaned = x1_coef_matrix_clean),
  "x1_se" = list(original = x1_se_matrix, cleaned = x1_se_matrix_clean)
)

for(matrix_name in names(matrices_info_x1)) {
  orig <- matrices_info_x1[[matrix_name]]$original
  clean <- matrices_info_x1[[matrix_name]]$cleaned
  
  orig_missing <- sum(is.na(orig))
  clean_missing <- sum(is.na(clean))
  additional_missing <- clean_missing - orig_missing
  
  cat(matrix_name, ": originally", orig_missing, "missing,", 
      additional_missing, "additional removed,", 
      clean_missing, "total missing\n")
}

cat("\nX1 outlier removal complete! Using cleaned matrices for plotting.\n")

# Load required libraries
library(ggplot2)
library(latex2exp)  # For TeX labels
library(gridExtra)
library(grid)
library(gtable)     # For legend creation
library(cowplot)    # For get_legend

# Source common functions for model renaming
if(file.exists("common_functions.R")) {
  source("common_functions.R")
} else {
  cat("Warning: common_functions.R not found. Using basic model renaming.\n")
  # Basic model renaming function
  rename_model <- function(model_name) {
    switch(model_name,
           "glm" = "GLM",
           "gee" = "GEE",
           "re_nosig" = "GAMLSS (Sig)",
           "re_np" = "GAMLSS (NP)",
           "lme4" = "LME4",
           "gamm" = "GAMM",
           "cop_n" = "Copula (N)",
           "cop" = "Copula",
           model_name)  # Return original if no match
  }
}

# Use the legend handling functions from main plotting file
# (get_legend and validate_legend are defined in the main plotting script)

# Function to prepare x1 bias plot data for a specific distribution
prepare_x1_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    non_missing_counts <- colSums(!is.na(x1_coef_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No x1 bias data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get true x1 values for this distribution
  true_x1 <- as.numeric(true_params_matrix$x1[dist_rows])
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated x1 values for this distribution and selected models
  # Note: x1 coefficients should not need link function transformations like mu1
  est_x1_matrix <- x1_coef_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  
  # Calculate relative bias using (estimated - true) / true
  bias_matrix <- est_x1_matrix / rep(true_x1, each = 1) - 1
  
  # Convert to long format for ggplot
  plot_data <- data.frame(
    corr_vals = rep(corr_vals, length(models_to_plot)),
    bias = as.vector(bias_matrix),
    model = rep(models_to_plot, each = length(true_x1)),
    model_label = rep(sapply(models_to_plot, rename_model), each = length(true_x1)),
    dist = dist_name
  )
  
  # Remove rows with missing values
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  return(plot_data)
}

# Create individual x1 bias plots for each distribution
create_x1_dist_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  plot_data <- prepare_x1_plot_data(dist_name, models_to_plot)
  
  if(is.null(plot_data) || nrow(plot_data) == 0) {
    cat("No data available for x1 bias plot:", dist_name, "\n")
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(plot_data$model_label)
  cat("Unique models in x1", dist_name, "data:", paste(unique_models, collapse = ", "), "\n")
  
  # Add line type mapping for different models
  plot_data$line_type <- ifelse(plot_data$model_label %in% c("GEE", "GLM"), "dotted",
                         ifelse(plot_data$model_label %in% c("GAMLSS", "GAMLSS NP", "LME4", "GAMM"), "dashed", "solid"))
  
  # Create the plot with consistent color mapping and line types for GJRM models
  p <- ggplot(plot_data, aes(x = corr_vals, y = bias, color = model_label, linetype = line_type)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    # Use consistent color mapping, only showing models present in this plot
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    # Add line type scale
    scale_linetype_identity() +
    coord_cartesian(xlim = c(0.2, 0.7),ylim=if(dist_name=="LO"){c(-0.2,1.2)}else{c(-0.2,.2)}) +
    labs(
      title = paste(if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name, "- Bias"),
      x = "Kendall's τ",
      y = TeX("$(\\hat{\\beta_{x_1}}/\\beta_{x_1})-1$"),
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
    cat("DEBUG: create_x1_dist_plot for", dist_name, "- plot created with legend position:", if(show_legend) "bottom" else "none", "\n")
    cat("DEBUG: create_x1_dist_plot for", dist_name, "- unique models in plot:", paste(unique_models, collapse = ", "), "\n")
    if(exists("model_colors") && length(model_colors) > 0) {
      available_colors <- intersect(unique_models, names(model_colors))
      cat("DEBUG: create_x1_dist_plot for", dist_name, "- available colors:", length(available_colors), "of", length(unique_models), "\n")
    } else {
      cat("DEBUG: create_x1_dist_plot for", dist_name, "- model_colors not available!\n")
    }
  }
  
  return(p)
}

# Function to prepare x1 SE plot data for a specific distribution
prepare_x1_se_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    non_missing_counts <- colSums(!is.na(x1_se_matrix_clean))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  
  if(length(dist_rows) == 0) {
    cat("No x1 SE data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated SE values for this distribution and selected models
  est_se_matrix <- x1_se_matrix_clean[dist_rows, models_to_plot, drop = FALSE]
  
  # Get theoretical SE values from true_ses_matrix if available
  true_se_vals <- rep(NA, length(dist_rows))
  
  if(exists("true_ses_matrix") && !is.null(true_ses_matrix)) {
    # Handle different structures of true_ses_matrix
    if(is.list(true_ses_matrix) && "true_ses_matrix" %in% names(true_ses_matrix)) {
      # If it's a list containing the matrix
      actual_ses_matrix <- true_ses_matrix$true_ses_matrix
    } else if(is.matrix(true_ses_matrix)) {
      # If it's already a matrix
      actual_ses_matrix <- true_ses_matrix
    } else {
      actual_ses_matrix <- NULL
    }
    
    if(!is.null(actual_ses_matrix) && is.matrix(actual_ses_matrix) && "x1_se" %in% colnames(actual_ses_matrix)) {
      true_se_vals <- actual_ses_matrix[dist_rows, "x1_se"]
      cat("Found theoretical x1 SEs for", dist_name, ":", sum(!is.na(true_se_vals)), "non-missing values\n")
    } else {
      cat("Warning: x1_se column not found in true_ses_matrix for", dist_name, "\n")
    }
  } else {
    cat("Warning: true_ses_matrix not available, using estimated SEs only for x1", dist_name, "\n")
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
    cat("No theoretical SEs available for x1", dist_name, ", using only estimated SEs\n")
  }
  
  # Remove rows with missing se_value
  combined_se_data <- combined_se_data[complete.cases(combined_se_data$se_value), ]
  
  return(combined_se_data)
}

# Function to create x1 SE plot for a specific distribution
create_x1_se_plot <- function(dist_name, models_to_plot = NULL, show_legend = TRUE) {
  
  se_data <- prepare_x1_se_plot_data(dist_name, models_to_plot)
  
  if(is.null(se_data) || nrow(se_data) == 0) {
    cat("No data available for x1 SE plot:", dist_name, "\n")
    return(NULL)
  }
  
  # Get unique model labels from the data
  all_models_in_data <- unique(se_data$model_label)
  estimated_models <- unique(se_data$model_label[se_data$type == "Estimated"])
  cat("Models in x1", dist_name, "SE data:", paste(all_models_in_data, collapse = ", "), "\n")
  
  # Add line type mapping for different models based on evaluation metrics approach
  se_data$line_type <- ifelse(se_data$type == "Theoretical", "dashed",
                       ifelse(se_data$model_label %in% c("GEE", "GLM"), "dotted",
                       ifelse(se_data$model_label %in% c("GAMLSS", "GAMLSS NP", "LME4", "GAMM"), "dashed", "solid")))
  
  # Create the plot with consistent color mapping and line types for GJRM models
  p <- ggplot(se_data, aes(x = corr_vals, y = se_value, color = model_label, linetype = line_type)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Manual color scale - only include models present in this plot
    scale_color_manual(values = c("Theoretical" = "black", model_colors[estimated_models]), 
                       breaks = all_models_in_data,
                       limits = all_models_in_data) +
    # Add line type scale
    scale_linetype_identity() +
    coord_cartesian(xlim = c(0.2, 0.7)) +
    labs(
      title = paste(if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name, "- SE"),
      x = "Kendall's τ",
      y = TeX("SE$(\\hat{\\beta_{x_1}})$"),
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

# MAIN PLOTTING SECTION ----

cat("\n=== Creating X1 Bias and SE Plots ===\n")

# Create x1 bias plots without individual legends (we'll add a shared legend)
plot_NO_x1_bias <- create_x1_dist_plot("NO", models_to_plot, show_legend = FALSE)
plot_PO_x1_bias <- create_x1_dist_plot("PO", models_to_plot, show_legend = FALSE) 
plot_GA_x1_bias <- create_x1_dist_plot("GA", models_to_plot, show_legend = FALSE)
plot_LO_x1_bias <- create_x1_dist_plot("LO", models_to_plot, show_legend = FALSE)

# Create x1 SE plots without individual legends
plot_NO_x1_se <- tryCatch(create_x1_se_plot("NO", models_to_plot, show_legend = FALSE), 
                          error = function(e) { cat("Error creating NO x1 SE plot:", e$message, "\n"); NULL })
plot_PO_x1_se <- tryCatch(create_x1_se_plot("PO", models_to_plot, show_legend = FALSE), 
                          error = function(e) { cat("Error creating PO x1 SE plot:", e$message, "\n"); NULL })
plot_GA_x1_se <- tryCatch(create_x1_se_plot("GA", models_to_plot, show_legend = FALSE), 
                          error = function(e) { cat("Error creating GA x1 SE plot:", e$message, "\n"); NULL })
plot_LO_x1_se <- tryCatch(create_x1_se_plot("LO", models_to_plot, show_legend = FALSE), 
                          error = function(e) { cat("Error creating LO x1 SE plot:", e$message, "\n"); NULL })

# Remove NULL plots and organize in pairs
bias_plots_x1 <- list(plot_NO_x1_bias, plot_PO_x1_bias, plot_GA_x1_bias, plot_LO_x1_bias)
se_plots_x1 <- list(plot_NO_x1_se, plot_PO_x1_se, plot_GA_x1_se, plot_LO_x1_se)

# Filter out NULL plots but keep the pairing
valid_indices_x1 <- which(!sapply(bias_plots_x1, is.null) & !sapply(se_plots_x1, is.null))
bias_plots_x1 <- bias_plots_x1[valid_indices_x1]
se_plots_x1 <- se_plots_x1[valid_indices_x1]

if(length(bias_plots_x1) > 0 && length(se_plots_x1) > 0) {
  # Use gridExtra for plot arrangement with shared legend (same as MU1/MU2 plots)
  library(gridExtra)
  library(grid)
  
  # Check what distributions are available for x1
  available_dists_x1 <- unique(true_params_matrix$dist)
  cat("Available distributions for x1:", paste(available_dists_x1, collapse = ", "), "\n")
  
  # Find the first distribution that has data for legend extraction (same approach as MU1/MU2)
  legend_dist_x1 <- NULL
  for(dist in c("GA", "NO", "PO", "LO")) {
    if(dist %in% available_dists_x1) {
      test_plot_x1 <- create_x1_dist_plot(dist, models_to_plot, show_legend = TRUE)
      
      if(!is.null(test_plot_x1)) {
        legend_dist_x1 <- dist
        break
      }
    }
  }
  
  if(is.null(legend_dist_x1)) {
    cat("Error: No valid distribution found for x1 legend extraction\n")
    # Force manual legend creation (same as MU1/MU2 approach)
    shared_legend_x1 <- validate_legend(NULL, "Chart X1", force_manual = TRUE)
  } else {
    cat("Using", legend_dist_x1, "distribution for x1 legend extraction\n")
    
    # Extract legend from the working distribution with improved error handling (same as MU1/MU2)
    tryCatch({
      legend_plot_temp_x1 <- create_x1_dist_plot(legend_dist_x1, models_to_plot, show_legend = TRUE) + 
        theme(
          legend.position = "bottom", 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 11),
          legend.box = "horizontal",
          legend.margin = margin(t = 10, b = 10)
        )
      
      # Extract the legend with improved method (same as MU1/MU2)
      shared_legend_x1 <- get_legend(legend_plot_temp_x1, force_manual = FALSE)
      
      # Validate and potentially replace with manual legend (same as MU1/MU2)
      shared_legend_x1 <- validate_legend(shared_legend_x1, "Chart X1", force_manual = FALSE)
      
      cat("Legend extraction completed for X1 plots\n")
      
    }, error = function(e) {
      cat("Error extracting x1 legend:", e$message, "\n")
      shared_legend_x1 <- validate_legend(NULL, "Chart X1", force_manual = TRUE)
    })
  }
  
  # Interleave bias and SE plots for each distribution
  all_plots_x1 <- vector("list", length(bias_plots_x1) * 2)
  for(i in 1:length(bias_plots_x1)) {
    all_plots_x1[2*i-1] <- bias_plots_x1[i]  # Bias plot (left column)
    all_plots_x1[2*i] <- se_plots_x1[i]      # SE plot (right column)
  }
  
  # Arrange plots in a grid (rows, 2 columns)
  plots_grid_x1 <- do.call(arrangeGrob, c(all_plots_x1, ncol = 2))
  
  # Combine the plots grid with the shared legend at the bottom
  combined_plot_x1 <- arrangeGrob(plots_grid_x1, shared_legend_x1, 
                                  ncol = 1, 
                                  heights = c(10, 1))
  
  # Display the combined plot
  grid.draw(combined_plot_x1)
  
  # Save the combined plot to Charts directory
  chart_filename_x1 <- paste0("Charts/X1_Bias_and_SE_vs_Correlation_Combined_", Sys.Date(), ".png")
  
  # Create Charts directory if it doesn't exist
  if(!dir.exists("Charts")) {
    dir.create("Charts")
  }
  
  ggsave(chart_filename_x1, plot = combined_plot_x1, width = 15, height = 14, dpi = 600)
  cat("Combined x1 bias and SE plot saved to:", chart_filename_x1, "\n")
  
  cat("X1 bias and SE vs correlation plots created successfully!\n")
  cat("Number of x1 bias plots created:", length(bias_plots_x1), "\n")
  cat("Number of x1 SE plots created:", length(se_plots_x1), "\n")
  
} else {
  cat("No x1 plots could be created - check data availability\n")
}

# Print diagnostic information for x1
cat("\n=== X1 DIAGNOSTIC INFORMATION ===\n")
cat("Distribution counts in true_params_matrix:\n")
print(table(true_params_matrix$dist))
cat("\nSkew matrix dimensions:", dim(skew_matrix), "\n")
cat("Available models in x1_coef_matrix:\n")
print(colnames(x1_coef_matrix))
cat("\nNon-missing values in x1_coef_matrix by model:\n")
for(model in models_to_plot) {
  if(model %in% colnames(x1_coef_matrix_clean)) {
    non_missing <- sum(!is.na(x1_coef_matrix_clean[, model]))
    model_label <- rename_model(model)
    cat(model_label, "(", model, "):", non_missing, "non-missing values\n")
  } else {
    cat(model, ": not found in x1 data\n")
  }
}

cat("\nNon-missing values in x1_se_matrix by model:\n")
for(model in models_to_plot) {
  if(model %in% colnames(x1_se_matrix_clean)) {
    non_missing <- sum(!is.na(x1_se_matrix_clean[, model]))
    model_label <- rename_model(model)
    cat(model_label, "(", model, "):", non_missing, "non-missing SE values\n")
  } else {
    cat(model, ": not found in x1 SE data\n")
  }
}

# Check true x1 parameter values
cat("\nTrue x1 parameter values:\n")
x1_values <- unique(true_params_matrix$x1)
cat("Unique true x1 values:", paste(sort(x1_values), collapse = ", "), "\n")

# Check if theoretical SEs are available
if(exists("true_ses_matrix") && !is.null(true_ses_matrix)) {
  if(is.list(true_ses_matrix) && "true_ses_matrix" %in% names(true_ses_matrix)) {
    actual_matrix <- true_ses_matrix$true_ses_matrix
  } else {
    actual_matrix <- true_ses_matrix
  }
  
  if(is.matrix(actual_matrix) && "x1_se" %in% colnames(actual_matrix)) {
    x1_se_availability <- sum(!is.na(actual_matrix[, "x1_se"]))
    cat("Theoretical x1 SEs available:", x1_se_availability, "out of", nrow(actual_matrix), "parameter combinations\n")
  } else {
    cat("No x1_se column found in theoretical SE matrix\n")
  }
} else {
  cat("No theoretical SE matrix available\n")
}

cat("\n=== X1 PLOTTING COMPLETED ===\n")
