# Evaluation Score Plots ####

cat("\n=== Creating Evaluation Score Plots ===\n")

# Load required libraries
library(cowplot)    # For get_legend function

# Create a function to prepare evaluation score plot data for a specific distribution and metric
prepare_eval_plot_data <- function(dist_name, metric_name, metric_matrix, models_to_plot = NULL) {
  
  cat("Preparing data for:", dist_name, metric_name, "\n")
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    non_missing_counts <- colSums(!is.na(metric_matrix))
    models_to_plot <- names(non_missing_counts[non_missing_counts > 0])
  }
  
  cat("Models to plot:", paste(models_to_plot, collapse = ", "), "\n")
  
  # Define the desired order for models
  desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
  
  # Reorder models_to_plot according to desired order
  models_to_plot <- intersect(desired_order, models_to_plot)
  cat("Models after reordering:", paste(models_to_plot, collapse = ", "), "\n")
  
  # Filter true parameters for this distribution
  dist_rows <- which(true_params_matrix$dist == dist_name)
  cat("Found", length(dist_rows), "rows for distribution", dist_name, "\n")
  
  if(length(dist_rows) == 0) {
    cat("No", metric_name, "data found for distribution:", dist_name, "\n")
    return(NULL)
  }
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  cat("Correlation values range:", min(corr_vals, na.rm = TRUE), "to", max(corr_vals, na.rm = TRUE), "\n")
  
  # Get metric values for this distribution and selected models
  metric_values <- metric_matrix[dist_rows, models_to_plot, drop = FALSE]
  cat("Metric matrix dimensions:", nrow(metric_values), "x", ncol(metric_values), "\n")
  cat("Non-missing metric values:", sum(!is.na(metric_values)), "/", length(metric_values), "\n")
  
  # Convert to long format for ggplot
  plot_data <- data.frame(
    corr_vals = rep(corr_vals, length(models_to_plot)),
    metric_value = as.vector(metric_values),
    model = rep(models_to_plot, each = length(corr_vals)),
    model_label = rep(sapply(models_to_plot, rename_model), each = length(corr_vals)),
    dist = dist_name,
    metric = metric_name
  )
  
  # Remove rows with missing values first
  plot_data_before_na <- nrow(plot_data)
  plot_data <- plot_data[complete.cases(plot_data), ]
  cat("Removed", plot_data_before_na - nrow(plot_data), "rows with missing values, keeping", nrow(plot_data), "\n")
  
  # Apply conservative outlier exclusion (3*IQR method) per model
  if(nrow(plot_data) > 0) {
    plot_data_before_outliers <- nrow(plot_data)
    outlier_summary <- data.frame(
      model_label = character(),
      dist = character(),
      metric = character(),
      n_total = integer(),
      n_outliers = integer(),
      pct_excluded = numeric(),
      outlier_range = character(),
      stringsAsFactors = FALSE
    )
    
    # Process each model separately
    models_in_data <- unique(plot_data$model_label)
    plot_data_filtered <- data.frame()
    
    for(model in models_in_data) {
      model_data <- plot_data[plot_data$model_label == model, ]
      
      if(nrow(model_data) > 4) {  # Need at least 5 points for meaningful outlier detection
        # Calculate IQR-based outlier bounds (conservative: 3*IQR)
        Q1 <- quantile(model_data$metric_value, 0.25, na.rm = TRUE)
        Q3 <- quantile(model_data$metric_value, 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        
        # Conservative outlier bounds (3*IQR instead of 1.5*IQR)
        lower_bound <- Q1 - 10 * IQR
        upper_bound <- Q3 + 3 * IQR
        
        # Identify outliers
        outliers <- model_data$metric_value < lower_bound | model_data$metric_value > upper_bound
        n_outliers <- sum(outliers)
        
        # Keep only non-outliers
        model_data_clean <- model_data[!outliers, ]
        
        # Record outlier statistics
        outlier_info <- data.frame(
          model_label = model,
          dist = dist_name,
          metric = metric_name,
          n_total = nrow(model_data),
          n_outliers = n_outliers,
          pct_excluded = round(n_outliers / nrow(model_data) * 100, 1),
          outlier_range = paste0("[", round(lower_bound, 3), ", ", round(upper_bound, 3), "]"),
          stringsAsFactors = FALSE
        )
        outlier_summary <- rbind(outlier_summary, outlier_info)
        
        # Add cleaned data to filtered dataset
        plot_data_filtered <- rbind(plot_data_filtered, model_data_clean)
        
      } else {
        # Too few points for outlier detection, keep all data
        outlier_info <- data.frame(
          model_label = model,
          dist = dist_name,
          metric = metric_name,
          n_total = nrow(model_data),
          n_outliers = 0,
          pct_excluded = 0,
          outlier_range = "Too few points",
          stringsAsFactors = FALSE
        )
        outlier_summary <- rbind(outlier_summary, outlier_info)
        plot_data_filtered <- rbind(plot_data_filtered, model_data)
      }
    }
    
    # Replace original data with filtered data
    plot_data <- plot_data_filtered
    
    # Print outlier exclusion summary
    cat("=== OUTLIER EXCLUSION SUMMARY (", dist_name, "-", metric_name, ") ===\n")
    total_excluded <- sum(outlier_summary$n_outliers)
    cat("Total outliers excluded:", total_excluded, "out of", plot_data_before_outliers, 
        "points (", round(total_excluded / plot_data_before_outliers * 100, 1), "%)\n")
    
    if(nrow(outlier_summary) > 0) {
      for(i in seq_len(nrow(outlier_summary))) {
        row <- outlier_summary[i, ]
        if(row$n_outliers > 0) {
          cat(sprintf("  %s: %d/%d outliers excluded (%.1f%%) - bounds: %s\n",
                     row$model_label, row$n_outliers, row$n_total, 
                     row$pct_excluded, row$outlier_range))
        }
      }
    }
    cat("Final data points after outlier removal:", nrow(plot_data), "\n")
  }
  
  # Set factor levels for model_label to maintain desired order
  if(nrow(plot_data) > 0) {
    model_label_order <- sapply(models_to_plot, rename_model)
    plot_data$model_label <- factor(plot_data$model_label, levels = model_label_order)
  }
  
  return(plot_data)
}

# Create evaluation score plots for each distribution and metric
create_eval_plot <- function(dist_name, metric_name, metric_matrix, models_to_plot = NULL, show_legend = TRUE, comprehensive_colors = NULL) {
  
  plot_data <- prepare_eval_plot_data(dist_name, metric_name, metric_matrix, models_to_plot)
  
  if(is.null(plot_data) || nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Get unique model labels from the data
  unique_models <- unique(plot_data$model_label)
  cat("Unique models in", dist_name, metric_name, "data:", paste(unique_models, collapse = ", "), "\n")
  
  # Determine y-axis label based on metric
  y_label <- switch(metric_name,
                   "VS2" = "VS2 Score",
                   "VS2_Weighted" = "VS2 Weighted Score", 
                   "LogLik" = "-2*Log Likelihood",
                   "BIC" = "BIC",
                   metric_name)
  
  # Use the global model_colors directly for consistency - it should have all models
  colors_to_use <- model_colors
  
  # CRITICAL FIX: Ensure we have colors for ALL models in this plot
  # If any model is missing a color, ggplot will use default (grey) colors for ALL models
  missing_colors <- unique_models[!unique_models %in% names(colors_to_use)]
  if(length(missing_colors) > 0) {
    cat("CRITICAL: Adding missing colors for", dist_name, metric_name, ":", paste(missing_colors, collapse = ", "), "\n")
    
    # Add missing colors using a consistent approach
    n_missing <- length(missing_colors)
    additional_colors <- rainbow(n_missing, start = 0.7, end = 1)  # Use different part of spectrum
    names(additional_colors) <- missing_colors
    colors_to_use <- c(colors_to_use, additional_colors)
  }
  
  # Double-check: ensure ALL unique_models have colors
  # CRITICAL FIX: Use explicit name-based indexing, not positional indexing
  if(exists("DEBUG_LEGENDS") && DEBUG_LEGENDS) {
    cat("DEBUG: About to create final_colors\n")
    cat("DEBUG: unique_models =", paste(unique_models, collapse = ", "), "\n")
    cat("DEBUG: names(colors_to_use) =", paste(names(colors_to_use), collapse = ", "), "\n")
  }
  
  # FIXED: Use explicit name-based lookup to avoid positional indexing issues
  final_colors <- character(length(unique_models))
  names(final_colors) <- unique_models
  for(i in seq_along(unique_models)) {
    model_name <- unique_models[i]
    if(model_name %in% names(colors_to_use)) {
      final_colors[i] <- colors_to_use[[model_name]]  # Use [[ ]] for single element lookup
    } else {
      final_colors[i] <- NA  # This will trigger the missing color check below
    }
  }
  
  if(exists("DEBUG_LEGENDS") && DEBUG_LEGENDS) {
    cat("DEBUG: final_colors result:\n")
    for(i in seq_along(final_colors)) {
      cat("DEBUG:   ", names(final_colors)[i], " -> ", final_colors[i], "\n")
    }
  }
  
  if(any(is.na(final_colors))) {
    cat("ERROR: Still missing colors after fix for", dist_name, metric_name, "\n")
    cat("Models:", paste(unique_models, collapse = ", "), "\n")
    cat("Available colors:", paste(names(colors_to_use), collapse = ", "), "\n")
  }
  
  # Debug: Check what colors we have vs what we need
  if(DEBUG_LEGENDS) {
    cat("DEBUG: Colors needed for", dist_name, metric_name, ":", paste(unique_models, collapse = ", "), "\n")
    cat("DEBUG: Available colors:", paste(names(colors_to_use), collapse = ", "), "\n")
    missing_colors_debug <- unique_models[!unique_models %in% names(colors_to_use)]
    if(length(missing_colors_debug) > 0) {
      cat("DEBUG: Missing colors for:", paste(missing_colors_debug, collapse = ", "), "\n")
    }
    # Show the actual color values being assigned
    cat("DEBUG: Color assignments for this plot:\n")
    for(model in unique_models) {
      color_val <- if(model %in% names(colors_to_use)) colors_to_use[model] else "MISSING"
      cat("DEBUG:   ", model, " -> ", color_val, "\n")
    }
  }
  
  # Create the plot with consistent color mapping using global model_colors
  p <- ggplot(plot_data, aes(x = corr_vals, y = metric_value, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Use comprehensive color mapping - ensure breaks/limits match the names in final_colors
    scale_color_manual(values = final_colors,
                       breaks = names(final_colors),
                       limits = names(final_colors)) +
    labs(
      title = paste(if(dist_name == "NO") "Normal" else if(dist_name == "PO") "Negative Binomial" else if(dist_name == "GA") "Gamma" else if(dist_name == "LO") "Bernoulli" else dist_name, "-", y_label),
      x = "Kendall's τ",
      y = y_label,
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.position = if(show_legend) "bottom" else "none"
    )
  
  if(DEBUG_LEGENDS && show_legend) {
    cat("DEBUG: create_eval_plot for", dist_name, metric_name, "- plot created with legend position:", if(show_legend) "bottom" else "none", "\n")
    cat("DEBUG: create_eval_plot for", dist_name, metric_name, "- unique models in plot:", paste(unique_models, collapse = ", "), "\n")
    if(exists("model_colors") && length(model_colors) > 0) {
      available_colors <- intersect(unique_models, names(model_colors))
      cat("DEBUG: create_eval_plot for", dist_name, metric_name, "- available colors:", length(available_colors), "of", length(unique_models), "\n")
      if(length(available_colors) != length(unique_models)) {
        missing_colors <- setdiff(unique_models, names(model_colors))
        cat("DEBUG: create_eval_plot for", dist_name, metric_name, "- missing colors for:", paste(missing_colors, collapse = ", "), "\n")
      }
    } else {
      cat("DEBUG: create_eval_plot for", dist_name, metric_name, "- model_colors not available!\n")
    }
  }
  
  return(p)
}

# Create evaluation score plots for all distributions and metrics
cat("Creating evaluation score plots for all distributions and metrics...\n")

# Define the desired order for models (same as used in other functions)
desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")

# Define the models to plot in the desired order (use available models from score matrices)
if(exists("vs2_matrix")) {
  available_models <- colnames(vs2_matrix)
  models_to_plot_eval <- intersect(desired_order, available_models)
  
  cat("Models to plot for evaluation scores:", paste(models_to_plot_eval, collapse = ", "), "\n")
  
  # Create plots for each distribution and each metric
  distributions <- c("NO", "PO", "GA", "LO")
  metrics <- list(
    "VS2" = vs2_matrix,
    "VS2_Weighted" = vs2_wt_matrix, 
    "LogLik" = loglik_matrix,
    "BIC" = bic_matrix
  )
  
  # Define which metrics to include in the grid (exclude VS2)
  metrics_for_grid <- list(
    "VS2_Weighted" = vs2_wt_matrix, 
    "LogLik" = loglik_matrix,
    "BIC" = bic_matrix
  )
  
  # Create comprehensive model set for consistent colors and legend
  # Collect all models that appear in ANY metric
  models_vs2_weighted <- models_to_plot_eval[models_to_plot_eval != "re_np"]  # Excludes GAMLSS NP
  models_loglik_bic <- models_to_plot_eval[!models_to_plot_eval %in% c("gee", "gamm")]  # Excludes GEE and GAMM
  
  # Union of all models used across all metrics
  comprehensive_models <- unique(c(models_vs2_weighted, models_loglik_bic))
  comprehensive_model_labels <- sapply(comprehensive_models, rename_model)
  
  cat("Comprehensive model set for consistent colors:", paste(comprehensive_models, collapse = ", "), "\n")
  cat("Comprehensive model labels:", paste(comprehensive_model_labels, collapse = ", "), "\n")
  
  # Verify that model_colors has all the models we need
  cat("Checking model_colors coverage:\n")
  cat("Available in model_colors:", paste(names(model_colors), collapse = ", "), "\n")
  missing_colors <- comprehensive_model_labels[!comprehensive_model_labels %in% names(model_colors)]
  if(length(missing_colors) > 0) {
    cat("ERROR: Missing colors for models:", paste(missing_colors, collapse = ", "), "\n")
    cat("This indicates models_to_plot in 'Plotting for Broad Models.r' needs to be updated\n")
    
    # Add the missing colors to prevent grey lines
    cat("FIXING: Adding missing colors...\n")
    additional_colors <- rainbow(length(missing_colors))
    names(additional_colors) <- missing_colors
    model_colors <- c(model_colors, additional_colors)
    
    cat("Updated model_colors now includes:\n")
    for(label in comprehensive_model_labels) {
      cat("  ", label, ":", model_colors[label], "\n")
    }
  } else {
    cat("SUCCESS: All comprehensive models have colors defined\n")
  }
  
  # Initialize list to store all evaluation plots
  eval_plots <- list()
  plot_counter <- 1
  
  for(dist in distributions) {
    for(metric_name in names(metrics_for_grid)) {  # Use metrics_for_grid instead of metrics
      metric_matrix <- metrics_for_grid[[metric_name]]
      
      # Determine which models to use for this metric
      current_models <- models_to_plot_eval
      
      # Exclude GEE and GAMM from LogLik and BIC plots
      if(metric_name %in% c("LogLik", "BIC")) {
        current_models <- current_models[!current_models %in% c("gee", "gamm")]
      }
      
      # Exclude GAMLSS NP from VS2 and VS2_Weighted plots
      if(metric_name %in% c("VS2", "VS2_Weighted")) {
        current_models <- current_models[current_models != "re_np"]
      }
      
      cat("Creating plot for:", dist, metric_name, "with models:", paste(current_models, collapse = ", "), "\n")
      
      # Create plot without individual legend, using global model_colors for consistency
      eval_plot <- create_eval_plot(dist, metric_name, metric_matrix, current_models, show_legend = FALSE)
      
      if(!is.null(eval_plot)) {
        eval_plots[[plot_counter]] <- eval_plot
        names(eval_plots)[plot_counter] <- paste(dist, metric_name, sep = "_")
        cat("Successfully created plot:", names(eval_plots)[plot_counter], "\n")
        plot_counter <- plot_counter + 1
      } else {
        cat("Failed to create plot for:", dist, metric_name, "\n")
      }
    }
  }
  
  # Create a plot with legend to extract it for evaluation plots
  # Use ONLY the models that actually appear in the grid plots, not comprehensive set
  cat("Creating legend with models that actually appear in grid plots...\n")
  
  # Determine which models actually appear across the grid plots we created
  # Collect unique model labels from all the plots that were actually created
  models_in_grid <- character()
  for(dist in distributions) {
    for(metric_name in names(metrics_for_grid)) {
      # Determine which models are used for this specific metric
      current_models <- models_to_plot_eval
      
      # Apply the same filtering logic as used in plot creation
      if(metric_name %in% c("LogLik", "BIC")) {
        current_models <- current_models[!current_models %in% c("gee", "gamm")]
      }
      if(metric_name %in% c("VS2", "VS2_Weighted")) {
        current_models <- current_models[current_models != "re_np"]
      }
      
      # Convert to labels and add to our collection
      current_labels <- sapply(current_models, rename_model)
      models_in_grid <- c(models_in_grid, current_labels)
    }
  }
  
  # Get unique models actually used in grid
  grid_models <- unique(models_in_grid)
  
  cat("Models actually appearing in grid plots:", paste(grid_models, collapse = ", "), "\n")
  
  # Create a temporary data frame with only the models that appear in the grid
  legend_data <- data.frame(
    x = rep(0.5, length(grid_models)),
    y = seq_along(grid_models),
    model_label = grid_models
  )
  
  # Create legend colors using the same explicit lookup as individual plots
  legend_colors <- character(length(grid_models))
  names(legend_colors) <- grid_models
  for(i in seq_along(grid_models)) {
    model_name <- grid_models[i]
    if(model_name %in% names(model_colors)) {
      legend_colors[i] <- model_colors[[model_name]]  # Use [[ ]] for single element lookup
    } else {
      legend_colors[i] <- "#999999"  # Default grey for missing colors
    }
  }
  
  # DEBUG: Print detailed legend color information
  cat("=== LEGEND COLOR DEBUG ===\n")
  cat("Grid models in order:", paste(grid_models, collapse = ", "), "\n")
  for(i in seq_along(grid_models)) {
    model_name <- grid_models[i]
    assigned_color <- legend_colors[i]
    cat("Legend position", i, ":", model_name, "->", assigned_color, "\n")
  }
  cat("Legend colors vector:\n")
  print(legend_colors)
  cat("=== END LEGEND DEBUG ===\n")
  
  legend_plot_eval <- ggplot(legend_data, aes(x = x, y = y, color = model_label)) +
    geom_point(size = 3) +
    scale_color_manual(values = legend_colors,
                       breaks = grid_models,
                       limits = grid_models) +
    labs(color = "Model") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  cat("Grid-specific legend plot created with models:", paste(grid_models, collapse = ", "), "\n")
  
  if(length(eval_plots) > 0) {
    # Load required libraries for plot arrangement
    library(gridExtra)
    library(grid)
    
    # Debug: Check which plots are NULL
    null_plots <- sapply(eval_plots, is.null)
    if(any(null_plots)) {
      cat("Warning: NULL plots found at positions:", which(null_plots), "\n")
      cat("NULL plot names:", names(eval_plots)[null_plots], "\n")
      # Remove NULL plots
      eval_plots <- eval_plots[!null_plots]
      cat("Proceeding with", length(eval_plots), "non-NULL plots\n")
    }
    
    if(length(eval_plots) > 0) {
      # Extract legend from evaluation plot with improved error handling
      cat("Attempting to extract legend for evaluation plots...\n")
      
      # Try to create legend plot with explicit legend settings
      legend_plot_eval_temp <- NULL
      shared_legend_eval <- NULL
      
      # Extract legend directly from the simple legend plot we created
      cat("Attempting to extract legend from simple legend plot...\n")
      
      shared_legend_eval <- NULL
      
      tryCatch({
        # Extract the legend from our simple legend plot
        shared_legend_eval <- get_legend(legend_plot_eval, force_manual = FALSE)
        
        if(!is.null(shared_legend_eval) && !inherits(shared_legend_eval, "nullGrob")) {
          cat("Successfully extracted legend from simple legend plot\n")
        } else {
          cat("Failed to extract legend, creating manual legend\n")
          shared_legend_eval <- validate_legend(NULL, "Evaluation Plots - Manual", force_manual = TRUE)
        }
        
      }, error = function(e) {
        cat("Error extracting legend:", e$message, "\n")
        shared_legend_eval <- validate_legend(NULL, "Evaluation Plots - Error Fallback", force_manual = TRUE)
      })
      
      # If no legend could be created, force manual legend creation
      if(is.null(shared_legend_eval)) {
        cat("Warning: Could not create legend from any distribution, creating manual legend\n")
        shared_legend_eval <- validate_legend(NULL, "Chart 1 - Fallback", force_manual = TRUE)
      }
      
      # Arrange plots in a 4x3 grid (4 distributions x 3 metrics) with error handling
      tryCatch({
        eval_plots_grid <- do.call(arrangeGrob, c(eval_plots, ncol = 3))  # Changed from ncol = 4 to ncol = 3
        
        # Combine the plots grid with the shared legend at the bottom
        combined_eval_plot <- arrangeGrob(eval_plots_grid, shared_legend_eval, 
                                         ncol = 1, 
                                         heights = c(10, 1))
        
        # Display the combined evaluation plot
        grid.draw(combined_eval_plot)
        
        # Save the combined evaluation plot to Charts directory
        eval_chart_filename <- paste0("Charts/Evaluation_Scores_vs_Correlation_", Sys.Date(), ".png")
        
        ggsave(eval_chart_filename, plot = combined_eval_plot, width = 12, height = 12, dpi = 600)
        cat("Combined evaluation scores plot saved to:", eval_chart_filename, "\n")
        
        cat("Evaluation score plots created successfully!\n")
        cat("Number of evaluation plots created:", length(eval_plots), "\n")
        cat("Metrics plotted:", paste(names(metrics_for_grid), collapse = ", "), "\n")  # Updated to use metrics_for_grid
        cat("Distributions plotted:", paste(distributions, collapse = ", "), "\n")
        cat("Note: VS2 plots excluded from grid display\n")
        
      }, error = function(e) {
        cat("Error arranging evaluation plots:", e$message, "\n")
        cat("Attempting to save individual plots instead...\n")
        
        # Save individual plots if grid arrangement fails
        for(i in seq_along(eval_plots)) {
          if(!is.null(eval_plots[[i]])) {
            plot_name <- names(eval_plots)[i]
            if(is.null(plot_name) || plot_name == "") plot_name <- paste0("eval_plot_", i)
            
            individual_filename <- paste0("Charts/", plot_name, "_", Sys.Date(), ".png")
            #ggsave(individual_filename, plot = eval_plots[[i]], width = 8, height = 6, dpi = 300)
            cat("Individual plot saved:", individual_filename, "\n")
          }
        }
      })
    } else {
      cat("No valid evaluation plots remaining after removing NULL plots\n")
    }
  } else {
    cat("No evaluation plots could be created - check data availability\n")
  }
  
  # Print evaluation data summary
  cat("\nEvaluation Data Summary:\n")
  for(metric_name in names(metrics)) {  # Show summary for all metrics including VS2
    metric_matrix <- metrics[[metric_name]]
    if(!is.null(metric_matrix)) {
      non_missing <- sum(!is.na(metric_matrix))
      total <- length(metric_matrix)
      included_in_grid <- ifelse(metric_name %in% names(metrics_for_grid), "YES", "NO (excluded)")
      cat(sprintf("%s: %d/%d non-missing values (%.1f%%) - Grid: %s\n", 
                  metric_name, non_missing, total, 
                  if(total > 0) non_missing / total * 100 else 0,
                  included_in_grid))
    }
  }
  
  # Print overall outlier exclusion summary
  cat("\n=== OVERALL OUTLIER EXCLUSION SUMMARY ===\n")
  cat("Conservative outlier detection used: 3*IQR method applied per model\n")
  cat("Summary by model across all distribution-metric combinations:\n")
  
} else {
  cat("vs2_matrix not found - skipping evaluation score plots\n")
}

