bt_mode <- FALSE  # Set to TRUE for BT files, FALSE for non-BT files

# EXTRACT ALL THE DATA FROM SAVED SIMULATIONS ----

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)

# Add debugging flag for legend issues
DEBUG_LEGENDS <- TRUE

cat("=== LEGEND DEBUGGING ENABLED ===\n")
cat("Libraries loaded successfully\n")
cat("cowplot available:", requireNamespace("cowplot", quietly = TRUE), "\n")

# Helper function to validate legend extraction and create manual legend if needed
validate_legend <- function(legend_grob, context = "", force_manual = FALSE) {
  if(DEBUG_LEGENDS) {
    cat("DEBUG:", context, "- validate_legend called, force_manual =", force_manual, "\n")
    cat("DEBUG:", context, "- Legend class:", paste(class(legend_grob), collapse = ", "), "\n")
  }
  
  # If force_manual is TRUE, create manual legend immediately
  if(force_manual) {
    if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Force manual legend creation\n")
    manual_legend <- create_model_legend()
    if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Manual legend created\n")
    return(manual_legend)
  }
  
  # Check if the legend is valid
  is_valid <- FALSE
  
  if(!is.null(legend_grob)) {
    if(inherits(legend_grob, "gtable") && nrow(legend_grob) > 0 && ncol(legend_grob) > 0) {
      is_valid <- TRUE
      if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Valid gtable with", nrow(legend_grob), "x", ncol(legend_grob), "dimensions\n")
    } else if(inherits(legend_grob, "grob") && !inherits(legend_grob, "zeroGrob") && !inherits(legend_grob, "nullGrob")) {
      is_valid <- TRUE
      if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Valid grob legend\n")
    } else {
      if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Invalid or empty legend (", paste(class(legend_grob), collapse = ", "), ")\n")
    }
  } else {
    if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Legend is NULL\n")
  }
  
  # If legend is invalid, create manual replacement
  if(!is_valid) {
    if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Creating manual replacement legend\n")
    manual_legend <- create_model_legend()
    return(manual_legend)
  }
  
  if(DEBUG_LEGENDS) cat("DEBUG:", context, "- Legend validation passed\n")
  return(legend_grob)
}

# Define robust get_legend function with manual legend creation capability
if(!exists("create_model_legend")) {
  # Load required libraries
  library(cowplot)
  library(ggplot2)
  library(grid)
  library(gtable)
  
  # Function to create a manual legend based on model colors and names
  create_model_legend <- function(model_names_input = NULL, model_colors_input = NULL) {
    tryCatch({
      if(DEBUG_LEGENDS) cat("DEBUG: Creating manual model legend\n")
      
      # Use provided inputs or fall back to global variables
      use_names <- model_names_input
      use_colors <- model_colors_input
      
      # If not provided, try to use global variables
      if(is.null(use_names) && exists("models_to_plot")) {
        use_names <- sapply(models_to_plot, function(x) if(exists("rename_model")) rename_model(x) else x)
      }
      if(is.null(use_colors) && exists("model_colors")) {
        use_colors <- model_colors
      }
      
      # Final fallback to common model names
      if(is.null(use_names)) {
        use_names <- c("GLM", "GEE", "GAMLSS", "LME4", "GAMM", "GJRM")
      }
      if(is.null(use_colors)) {
        use_colors <- rainbow(length(use_names))
        names(use_colors) <- use_names
      }
      
      if(DEBUG_LEGENDS) cat("DEBUG: Using", length(use_names), "models for legend:", paste(use_names, collapse = ", "), "\n")
      
      # Ensure we have colors for all models
      legend_colors <- rep("#000000", length(use_names))  # default black
      names(legend_colors) <- use_names
      
      # Map available colors
      for(i in seq_along(use_names)) {
        model_name <- use_names[i]
        if(!is.null(use_colors) && model_name %in% names(use_colors)) {
          legend_colors[i] <- use_colors[model_name]
        } else if(!is.null(use_colors) && length(use_colors) >= i) {
          legend_colors[i] <- use_colors[i]
        }
      }
      
      # Create legend using grid functions
      n_models <- length(use_names)
      
      # Create a simple legend grob
      if(requireNamespace("grid", quietly = TRUE)) {
        legend_grob <- grid::legendGrob(
          labels = use_names,
          pch = rep(16, n_models),  # filled circles
          gp = grid::gpar(col = legend_colors, fontsize = 9),
          byrow = TRUE,
          ncol = min(3, ceiling(n_models/2))  # Arrange in columns
        )
        
        if(DEBUG_LEGENDS) cat("DEBUG: Manual legend created successfully with grid::legendGrob\n")
        return(legend_grob)
      }
      
      # If grid legendGrob fails, create using textGrob
      legend_text <- paste("Models:", paste(use_names, collapse = ", "))
      return(grid::textGrob(legend_text, gp = grid::gpar(fontsize = 8)))
      
    }, error = function(e) {
      if(DEBUG_LEGENDS) cat("DEBUG: Manual legend creation failed:", e$message, "\n")
      return(grid::textGrob("Model Legend", gp = grid::gpar(fontsize = 8)))
    })
  }
  
  # Main get_legend function
  get_legend <- function(plot, force_manual = FALSE) {
    if(DEBUG_LEGENDS) cat("DEBUG: get_legend called, force_manual =", force_manual, "\n")
    
    # If force_manual is TRUE, skip extraction and create manual legend
    if(force_manual) {
      if(DEBUG_LEGENDS) cat("DEBUG: Forcing manual legend creation\n")
      return(create_model_legend())
    }
    
    tryCatch({
      # Ensure the plot exists
      if(is.null(plot)) {
        if(DEBUG_LEGENDS) cat("DEBUG: NULL plot provided\n")
        return(create_model_legend())
      }
      
      if(DEBUG_LEGENDS) cat("DEBUG: Plot class:", paste(class(plot), collapse = ", "), "\n")
      
      # Force legend to visible position for extraction
      plot_with_legend <- plot + theme(legend.position = "bottom")
      if(DEBUG_LEGENDS) cat("DEBUG: Set legend position to bottom\n")
      
      # Method 1: Try cowplot approach
      if(requireNamespace("cowplot", quietly = TRUE)) {
        if(DEBUG_LEGENDS) cat("DEBUG: Attempting cowplot::get_legend\n")
        legend_grob <- cowplot::get_legend(plot_with_legend)
        
        if(!is.null(legend_grob) && !inherits(legend_grob, "zeroGrob") && 
           inherits(legend_grob, "gtable") && nrow(legend_grob) > 0) {
          if(DEBUG_LEGENDS) cat("DEBUG: cowplot extraction successful\n")
          return(legend_grob)
        } else {
          if(DEBUG_LEGENDS) cat("DEBUG: cowplot returned empty/invalid legend\n")
        }
      }
      
      # Method 2: Try ggplotGrob approach
      if(DEBUG_LEGENDS) cat("DEBUG: Attempting ggplotGrob extraction\n")
      plot_grob <- ggplotGrob(plot_with_legend)
      
      # Look for guide elements
      guide_indices <- which(grepl("guide", plot_grob$layout$name, ignore.case = TRUE))
      if(DEBUG_LEGENDS) cat("DEBUG: Found", length(guide_indices), "guide elements\n")
      
      if(length(guide_indices) > 0) {
        for(idx in guide_indices) {
          legend_candidate <- plot_grob$grobs[[idx]]
          if(inherits(legend_candidate, "gtable") && nrow(legend_candidate) > 0) {
            if(DEBUG_LEGENDS) cat("DEBUG: ggplotGrob extraction successful at index", idx, "\n")
            return(legend_candidate)
          }
        }
      }
      
      # Method 3: Look for any gtable that might be a legend
      if(DEBUG_LEGENDS) cat("DEBUG: Searching all grobs for potential legends\n")
      for(i in seq_along(plot_grob$grobs)) {
        grob <- plot_grob$grobs[[i]]
        if(inherits(grob, "gtable") && nrow(grob) > 1 && ncol(grob) > 0) {
          # Check if it contains legend-like elements
          grob_names <- if(!is.null(names(grob$grobs))) names(grob$grobs) else ""
          if(any(grepl("legend|key|guide", grob_names, ignore.case = TRUE))) {
            if(DEBUG_LEGENDS) cat("DEBUG: Found potential legend grob at index", i, "\n")
            return(grob)
          }
        }
      }
      
      if(DEBUG_LEGENDS) cat("DEBUG: No legend found in plot, creating manual legend\n")
      
    }, error = function(e) {
      if(DEBUG_LEGENDS) cat("DEBUG: Error in legend extraction:", e$message, "\n")
    })
    
    # If all extraction methods fail, create manual legend
    return(create_model_legend())
  }
}

# Set working directory and paths
data_path <- ifelse(bt_mode, "Data/broad_sims_bt", "Data/broad_sims")
rdata_files <- list.files(data_path, pattern = "\\.RData$", full.names = TRUE)
cat("Loading", length(rdata_files), "files from:", data_path, "\n")

# Initialize containers for results
all_results <- list()
parameter_info <- list()

# Model names (should match the column names in your results)
model_names <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
                "cop", "cop_n", "cop_j", "cop_g", "cop_f",
                "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t")

cat("Loading all RData files...\n")

# Load all RData files and extract results
for (i in seq_along(rdata_files)) {
  file_path <- rdata_files[i]
  file_name <- basename(file_path)
  
  if (i %% 100 == 0) {
    cat("Processing file", i, "of", length(rdata_files), "\n")
  }
  
  # Load the RData file
  tryCatch({
    load(file_path)
    
    # Extract parameter information from filename
    # Expected format: CoefSimData_DIST_a_b_c_mu1_mu2_x1_x2_n_outer_sims_date.RData
    file_parts <- strsplit(gsub("\\.RData$", "", file_name), "_")[[1]]
    
    if (length(file_parts) >= 12) {
      param_info <- list(
        file_name = file_name,
        dist = file_parts[2],
        a = as.numeric(file_parts[3]),
        b = as.numeric(file_parts[4]),
        c = as.numeric(file_parts[5]),
        mu1 = as.numeric(file_parts[6]),
        mu2 = as.numeric(file_parts[7]),
        x1 = as.numeric(file_parts[8]),
        x2 = as.numeric(file_parts[9]),
        n = as.numeric(file_parts[10]),
        outer_sims = as.numeric(file_parts[11])
      )
      
      # Convert "NA" strings to actual NA values for specific parameters only
      for (param in c("a", "b", "c")) {
        if (is.character(param_info[[param]]) && param_info[[param]] == "NA") {
          param_info[[param]] <- NA
        }
      }
      
      parameter_info[[i]] <- param_info
      
      # Store the loaded results
      all_results[[i]] <- list(
        par_estimates = par_estimates,
        score_items = score_items,
        df_values = df_values,
        times = times,
        conv = conv,
        params = param_info
      )
      
    } else {
      cat("Warning: Could not parse filename:", file_name, "\n")
    }
    
  }, error = function(e) {
    cat("Error loading file", file_path, ":", e$message, "\n")
  })
}

# Remove any NULL entries
all_results <- all_results[!sapply(all_results, is.null)]
parameter_info <- parameter_info[!sapply(parameter_info, is.null)]

cat("Successfully loaded", length(all_results), "result files\n")

# Create row names for matrices (parameter set identifiers)
row_names <- sapply(parameter_info, function(p) {
  paste(p$dist, 
        ifelse(is.na(p$a), "NA", p$a),
        ifelse(is.na(p$b), "NA", p$b), 
        ifelse(is.na(p$c), "NA", p$c),
        p$mu1, p$mu2, p$x1, p$x2, p$n, 
        sep = "_")
})

# Initialize result matrices
n_models <- length(model_names)
n_param_sets <- length(all_results)

cat("Creating result matrices for", n_param_sets, "parameter sets and", n_models, "models\n")

# Evaluation Score Matrices
vs2_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                     dimnames = list(row_names, model_names))
vs2_wt_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                        dimnames = list(row_names, model_names))
loglik_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                        dimnames = list(row_names, model_names))
bic_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                     dimnames = list(row_names, model_names))

# Coefficient Matrices
mu1_coef_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                          dimnames = list(row_names, model_names))
mu2_coef_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                          dimnames = list(row_names, model_names))
x1_coef_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                         dimnames = list(row_names, model_names))
x2_coef_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                         dimnames = list(row_names, model_names))

# Standard Error Matrices
mu1_se_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                        dimnames = list(row_names, model_names))
mu2_se_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                        dimnames = list(row_names, model_names))
#TEMP FIX NEED TO FLOW BACK TO MAIN FITTING FUNCTION

x1_se_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                       dimnames = list(row_names, model_names))
x2_se_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                       dimnames = list(row_names, model_names))

# Timing Matrix
timing_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                        dimnames = list(row_names, model_names))

# Convergence Matrix
convergence_matrix <- matrix(NA, nrow = n_param_sets, ncol = n_models,
                             dimnames = list(row_names, model_names))

# True Parameter Values Matrix
param_names <- c("a", "b", "c", "mu1", "mu2", "x1", "x2", "n", "dist")
true_params_matrix <- matrix(NA, nrow = n_param_sets, ncol = length(param_names),
                             dimnames = list(row_names, param_names))

cat("Filling matrices with results...\n")

# Define which models actually have scores available (not coefficients)
available_score_models <- c("glm", "gee", "re_nosig", "lme4", "gamm",
                           "cop", "cop_n", "cop_j", "cop_g", "cop_f")

# Get indices of available models for scores in the full model_names vector
available_score_indices <- which(model_names %in% available_score_models)

# Fill matrices with data
for (i in seq_along(all_results)) {
  result <- all_results[[i]]
  
  if (i %% 100 == 0) {
    cat("Filling matrices: processed", i, "of", length(all_results), "results\n")
  }
  
  tryCatch({
    # Debug: Check structure of current result
    if (i == 1) {
      cat("Debug - First result structure:\n")
      cat("par_estimates$mu1 length:", length(result$par_estimates$mu1), "\n")
      cat("par_estimates$mu1 class:", class(result$par_estimates$mu1), "\n")
      if (is.matrix(result$par_estimates$mu1)) {
        cat("par_estimates$mu1 dimensions:", dim(result$par_estimates$mu1), "\n")
      }
    }
    
    # Evaluation scores - only fill available score model columns
    if ("vs2" %in% names(result$score_items)) {
      vs2_matrix[i, available_score_indices] <- result$score_items$vs2
    }
    if ("vs2_wt" %in% names(result$score_items)) {
      vs2_wt_matrix[i, available_score_indices] <- result$score_items$vs2_wt
    }
    if ("logliks" %in% names(result$score_items)) {
      loglik_matrix[i,] <- -2 * result$score_items$logliks
    }
    
    # Calculate BIC: -2*loglik + k*log(n) - only for available score models
    if ("logliks" %in% names(result$score_items) && !is.null(result$df_values$logliks)) {
      log_lik <- result$score_items$logliks
      df <- result$df_values$logliks
      n <- result$params$n
      bic_matrix[i, ] <- -2 * log_lik + df * log(n)
    }
    
    # Coefficients - handle both vector and matrix cases
    if ("mu1" %in% names(result$par_estimates)) {
      coef_data <- result$par_estimates$mu1
      if (is.matrix(coef_data)) {
        coef_data <- as.vector(coef_data[1, ])  # Take first row if matrix
      }
      if (length(coef_data) == length(model_names)) {
        mu1_coef_matrix[i, ] <- coef_data
      } else {
        cat("Warning: mu1 length mismatch for result", i, "- expected", length(model_names), "got", length(coef_data), "\n")
      }
    }
    if ("mu2" %in% names(result$par_estimates)) {
      coef_data <- result$par_estimates$mu2
      if (is.matrix(coef_data)) {
        coef_data <- as.vector(coef_data[1, ])
      }
      if (length(coef_data) == length(model_names)) {
        mu2_coef_matrix[i, ] <- coef_data
      }
    }
    if ("x1" %in% names(result$par_estimates)) {
      coef_data <- result$par_estimates$x1
      if (is.matrix(coef_data)) {
        coef_data <- as.vector(coef_data[1, ])
      }
      if (length(coef_data) == length(model_names)) {
        x1_coef_matrix[i, ] <- coef_data
      }
    }
    if ("x2" %in% names(result$par_estimates)) {
      coef_data <- result$par_estimates$x2
      if (is.matrix(coef_data)) {
        coef_data <- as.vector(coef_data[1, ])
      }
      if (length(coef_data) == length(model_names)) {
        x2_coef_matrix[i, ] <- coef_data
      }
    }
    
    # Standard errors - handle both vector and matrix cases
    if ("mu1_se" %in% names(result$par_estimates)) {
      se_data <- result$par_estimates$mu1_se
      if (is.matrix(se_data)) {
        se_data <- as.vector(se_data[1, ])
      }
      if (length(se_data) == length(model_names)) {
        mu1_se_matrix[i, ] <- se_data
      }
    }
    if ("mu2_se" %in% names(result$par_estimates)) {
      se_data <- result$par_estimates$mu2_se
      if (is.matrix(se_data)) {
        se_data <- as.vector(se_data[1, ])
      }
      if (length(se_data) == length(model_names)) {
        mu2_se_matrix[i, ] <- se_data
      }
    }
    if ("x1_se" %in% names(result$par_estimates)) {
      se_data <- result$par_estimates$x1_se
      if (is.matrix(se_data)) {
        se_data <- as.vector(se_data[1, ])
      }
      if (length(se_data) == length(model_names)) {
        x1_se_matrix[i, ] <- se_data
      }
    }
    if ("x2_se" %in% names(result$par_estimates)) {
      se_data <- result$par_estimates$x2_se
      if (is.matrix(se_data)) {
        se_data <- as.vector(se_data[1, ])
      }
      if (length(se_data) == length(model_names)) {
        x2_se_matrix[i, ] <- se_data
      }
    }
    
    # Timing - only available score models
    if (!is.null(result$times) && is.matrix(result$times)) {
      timing_matrix[i, ] <- as.numeric(result$times[1, ])
    }
    
    # Convergence - only available score models
    if (!is.null(result$conv) && is.matrix(result$conv)) {
      convergence_matrix[i, ] <- result$conv[1, ]
    }
    
    # Fill true parameter values from filename
    true_params_matrix[i, "a"] <- result$params$a
    true_params_matrix[i, "b"] <- result$params$b
    true_params_matrix[i, "c"] <- result$params$c
    true_params_matrix[i, "mu1"] <- result$params$mu1
    true_params_matrix[i, "mu2"] <- result$params$mu2
    true_params_matrix[i, "x1"] <- result$params$x1
    true_params_matrix[i, "x2"] <- result$params$x2
    true_params_matrix[i, "n"] <- result$params$n
    true_params_matrix[i, "dist"] <- result$params$dist

    
  }, error = function(e) {
    cat("Error processing result", i, ":", e$message, "\n")
  })
}
true_params_matrix=data.frame(true_params_matrix)
true_params_matrix

cat("Matrix creation complete!\n")

# Print summary information
cat("\n=== SUMMARY ===\n")
cat("Parameter sets loaded:", n_param_sets, "\n")
cat("Models per parameter set:", n_models, "\n")
cat("Total data points:", n_param_sets * n_models, "\n")

# Show distribution by distribution type
dist_table <- table(sapply(parameter_info, function(x) x$dist))
cat("\nDistribution breakdown:\n")
print(dist_table)

# Show sample size breakdown
n_table <- table(sapply(parameter_info, function(x) x$n))
cat("\nSample size breakdown:\n")
print(n_table)

# Check for missing data
cat("\nMissing data summary:\n")
cat("VS2 scores missing:", sum(is.na(vs2_matrix)), "/", length(vs2_matrix), "\n")
cat("VS2 weighted missing:", sum(is.na(vs2_wt_matrix)), "/", length(vs2_wt_matrix), "\n")
cat("Log-likelihood missing:", sum(is.na(loglik_matrix)), "/", length(loglik_matrix), "\n")
cat("BIC missing:", sum(is.na(bic_matrix)), "/", length(bic_matrix), "\n")
cat("Coefficient estimates missing:", sum(is.na(mu1_coef_matrix)), "/", length(mu1_coef_matrix), "\n")

cat("\nMatrices created and ready for plotting!\n")
cat("Available matrices:\n")
cat("- Evaluation: vs2_matrix, vs2_wt_matrix, loglik_matrix, bic_matrix\n")
cat("- Coefficients: mu1_coef_matrix, mu2_coef_matrix, x1_coef_matrix, x2_coef_matrix\n")
cat("- Standard Errors: mu1_se_matrix, mu2_se_matrix, x1_se_matrix, x2_se_matrix\n")
cat("- Other: timing_matrix, convergence_matrix\n")
cat("- True Parameters: true_params_matrix (columns: a, b, c, mu1, mu2, x1, x2, n)\n")
cat("- Parameter info: parameter_info (list with parameter details for each set)\n")
#Calculate skewness and correlation for each parameter set
set.seed(1)
skew_matrix=matrix(NA,nrow=nrow(true_params_matrix),ncol=2)
for (i in 1:nrow(true_params_matrix)) {

  print(c(i,nrow(true_params_matrix)))
  source("common_functions.R")
  library(e1071)
  library(MASS)
  
  parameters=true_params_matrix[i,]
  
  # Run 10 simulations and take the average
  skew_vals <- numeric(10)
  corr_vals <- numeric(10)
  
  for (sim in 1:10) {
    simData=generateBivDist_withCov(dist=parameters["dist"],
                            a=as.numeric(parameters["a"]),
                            b=as.numeric(parameters["b"]),
                            c=as.numeric(parameters["c"]),
                            mu1=as.numeric(parameters["mu1"]),
                            mu2=as.numeric(parameters["mu2"]),
                            x1=0, #Calculating skew without covariates
                            x2=0
                            ,n=as.numeric(parameters["n"]))

    t1=simData[simData$time==1,"random_variable"]
    t2=simData[simData$time==0,"random_variable"]

    corr_vals[sim] <- cor(t1,t2,method="kendall")
    skew1 <- skewness(t1)
    skew2 <- skewness(t2)
    skew_vals[sim] <- (skew1+skew2)/2
  }
  
  # Store the averages
  skew_matrix[i,] <- c(mean(skew_vals), mean(corr_vals))
}

load("input_params_list.RData")

# Convert input_params_list to matrix format
cat("Converting input_params_list to matrix format...\n")

# Get the parameter names from the first list element
param_names <- names(input_params_list[[1]])
n_params <- length(input_params_list)

cat("Found", n_params, "parameter sets with columns:", paste(param_names, collapse = ", "), "\n")

# Initialize matrix
input_params_matrix <- matrix(NA, nrow = n_params, ncol = length(param_names))
colnames(input_params_matrix) <- param_names

# Fill matrix row by row
for (i in 1:n_params) {
  for (j in 1:length(param_names)) {
    param_name <- param_names[j]
    value <- input_params_list[[i]][[param_name]]
    
    # Handle different data types
    if (is.null(value)) {
      input_params_matrix[i, j] <- NA
    } else if (is.character(value)) {
      # For character values (like "dist"), we'll convert to character later
      input_params_matrix[i, j] <- value
    } else {
      input_params_matrix[i, j] <- as.numeric(value)
    }
  }
}

# Convert to data frame for better handling of mixed data types
input_params_df <- data.frame(input_params_matrix, stringsAsFactors = FALSE)

# Convert numeric columns to proper numeric type
numeric_cols <- c("a", "b", "c", "mu1", "mu2", "x1", "x2", "n")
for (col in numeric_cols) {
  if (col %in% colnames(input_params_df)) {
    input_params_df[[col]] <- as.numeric(input_params_df[[col]])
  }
}

# Keep dist as character
if ("dist" %in% colnames(input_params_df)) {
  input_params_df$dist <- as.character(input_params_df$dist)
}

cat("Matrix conversion complete!\n")
cat("Dimensions:", nrow(input_params_df), "rows x", ncol(input_params_df), "columns\n")

# Show first few rows
cat("\nFirst 5 rows of input_params_df:\n")
print(head(input_params_df, 5))

# Show summary
cat("\nSummary of parameter distributions:\n")
if ("dist" %in% colnames(input_params_df)) {
  cat("Distribution types:\n")
  print(table(input_params_df$dist))
}

if ("n" %in% colnames(input_params_df)) {
  cat("\nSample sizes:\n")
  print(table(input_params_df$n))
}

cat("\nColumn types:\n")
print(sapply(input_params_df, class))

parameters=input_params_df

# Load true values from file if available ----
# Load newly computed true SEs (if available) to override the old true_SE_out
# This uses the get_true_ses function results computed at the end of this script
#tryCatch({
#  cat("Attempting to load computed true SEs...\n")
#  
#  # Try to find the most recent true_ses_complete file
#  cache_files <- list.files("Cache", pattern = "true_ses_complete_.*\\.rds$", full.names = TRUE)
#  if(length(cache_files) > 0) {
#    # Use the most recent file
#    cache_file <- cache_files[which.max(file.mtime(cache_files))]
#    cat("Loading true SEs from:", basename(cache_file), "\n")
#    
#    final_data <- readRDS(cache_file)
#    true_ses_matrix <- final_data$true_ses_matrix
#    cat("Successfully loaded true_ses_matrix with dimensions:", dim(true_ses_matrix), "\n")
#    cat("Available SE columns:", colnames(true_ses_matrix), "\n")
#    cat("Non-missing t1_se values:", sum(!is.na(true_ses_matrix[, "t1_se"])), "\n")
#  } else {
#    cat("No true_ses_complete files found in Cache directory\n")
#    true_ses_matrix <- NULL
#  }
#}, error = function(e) {
#  cat("Could not load computed true SEs:", e$message, "\n")
#  cat("Will use old true_SE_out method instead\n")
#  true_ses_matrix <- NULL
#})
# Compute true values (new simulations) ----

# Check if we already have a valid true_ses_matrix from the previous loading section
if(exists("true_ses_matrix") && !is.null(true_ses_matrix) && 
   is.matrix(true_ses_matrix) && sum(!is.na(true_ses_matrix)) > 0) {
  cat("\n=== Using Previously Loaded True Standard Errors ===\n")
  cat("Found valid true_ses_matrix with dimensions:", dim(true_ses_matrix), "\n")
  cat("Non-missing values:", sum(!is.na(true_ses_matrix)), "/", length(true_ses_matrix), "\n")
  cat("Skipping computation as valid data already exists.\n")
  
  # Apply the transformations that would normally be done at the end
  cat("Applying transformations to loaded true_ses_matrix...\n")
  true_ses_matrix <- sqrt(true_ses_matrix)
  true_ses_matrix[true_params_matrix[,"dist"]=="GA",] <- true_ses_matrix[true_params_matrix[,"dist"]=="GA",]/sqrt(1000)
  cat("Transformations applied successfully.\n")
  
} else {
  cat("\n=== Computing True Standard Errors for All Parameter Combinations ===\n")
  cat("No valid true_ses_matrix found, starting computation...\n")
  
  # Load the get_true_ses function
  source("common_functions.R")
  
  # Create a filename for saving progress
  progress_filename <- paste0("Cache/true_ses_progress_", Sys.Date(), ".rds")
  final_filename <- paste0("Cache/true_ses_complete_", Sys.Date(), ".rds")
  
  # Check if we have previous progress to restore
  if(file.exists(progress_filename)) {
    cat("Found previous progress file, loading...\n")
    progress_data <- readRDS(progress_filename)
    true_ses_matrix <- progress_data$true_ses_matrix
    completed_rows <- progress_data$completed_rows
    start_row <- max(completed_rows) + 1
    cat("Resuming from row", start_row, "of", nrow(true_params_matrix), "\n")
  } else {
    cat("Starting fresh computation...\n")
    # Initialize matrix to store true standard errors
    # get_true_ses returns SEs for: t1, t2, x1, x2
    true_ses_matrix <- matrix(NA, nrow = nrow(true_params_matrix), ncol = 4)
    colnames(true_ses_matrix) <- c("t1_se", "t2_se", "x1_se", "x2_se")
    rownames(true_ses_matrix) <- rownames(true_params_matrix)
    completed_rows <- integer(0)
    start_row <- 1
  }

cat("Total parameter combinations to process:", nrow(true_params_matrix), "\n")
cat("Starting from row:", start_row, "\n")

# Initialize timing variables for time estimation
computation_start_time <- Sys.time()
row_times <- numeric()

# Process each row
for(i in start_row:nrow(true_params_matrix)) {
  
  # Record start time for this row
  row_start_time <- Sys.time()
  
  cat("Processing row", i, "of", nrow(true_params_matrix), "...\n")
  
  # Get parameters for this row
  params <- true_params_matrix[i, ]
  
  # Print current parameter combination
  cat("  Distribution:", params$dist, 
      "| a:", params$a, 
      "| b:", params$b, 
      "| c:", params$c,
      "| mu1:", params$mu1, 
      "| mu2:", params$mu2, 
      "| n:", params$n, "\n")
  
  tryCatch({
    # Call get_true_ses function
    true_ses_result <- get_true_ses(
      n = as.numeric(params$n),
      a = as.numeric(params$a),
      b = as.numeric(params$b), 
      c = as.numeric(params$c),
      mu1 = as.numeric(params$mu1),
      mu2 = as.numeric(params$mu2),
      dist = params$dist,
      x1 = as.numeric(params$x1),
      x2 = as.numeric(params$x2)
    )
    
    # Store results in matrix
    if(is.vector(true_ses_result) && length(true_ses_result) >= 4) {
      true_ses_matrix[i, 1:4] <- true_ses_result[1:4]
    } else if(is.list(true_ses_result) || (is.vector(true_ses_result) && !is.null(names(true_ses_result)))) {
      # If result is a named vector or list, extract by names
      if("t1" %in% names(true_ses_result)) {
        true_ses_matrix[i, "t1_se"] <- true_ses_result[["t1"]]
      }
      if("t2" %in% names(true_ses_result)) {
        true_ses_matrix[i, "t2_se"] <- true_ses_result[["t2"]]
      }
      if("x1" %in% names(true_ses_result)) {
        true_ses_matrix[i, "x1_se"] <- true_ses_result[["x1"]]
      }
      if("x2" %in% names(true_ses_result)) {
        true_ses_matrix[i, "x2_se"] <- true_ses_result[["x2"]]
      }
    } else {
      cat("    Warning: Unexpected result format from get_true_ses\n")
      cat("    Result class:", class(true_ses_result), "\n")
      cat("    Result length:", length(true_ses_result), "\n")
      if(!is.null(names(true_ses_result))) {
        cat("    Result names:", paste(names(true_ses_result), collapse = ", "), "\n")
      }
    }
    
    # Update completed rows
    completed_rows <- c(completed_rows, i)
    
    cat("    Success! SEs:", 
        "t1_se =", round(true_ses_matrix[i, "t1_se"], 6),
        "| t2_se =", round(true_ses_matrix[i, "t2_se"], 6), 
        "| x1_se =", round(true_ses_matrix[i, "x1_se"], 6),
        "| x2_se =", round(true_ses_matrix[i, "x2_se"], 6), "\n")
    
  }, error = function(e) {
    cat("    Error computing true SEs:", e$message, "\n")
    # Keep NAs in the matrix for failed computations
  })
  
  # Record end time for this row and update timing estimates
  row_end_time <- Sys.time()
  row_duration <- as.numeric(row_end_time - row_start_time, units = "secs")
  row_times <- c(row_times, row_duration)
  
  # Keep only last 50 row times for moving average (more responsive to recent performance)
  if(length(row_times) > 50) {
    row_times <- tail(row_times, 50)
  }
  
  # Calculate time estimates
  avg_time_per_row <- mean(row_times)
  rows_completed <- i - start_row + 1
  rows_remaining <- nrow(true_params_matrix) - i
  estimated_time_remaining <- rows_remaining * avg_time_per_row
  elapsed_time <- as.numeric(Sys.time() - computation_start_time, units = "mins")
  
  # Display time estimate every 10 rows or for the first few rows
  if(i %% 10 == 0 || i <= start_row + 5) {
    cat("    Time: Row", i, "took", round(row_duration, 1), "sec |",
        "Avg:", round(avg_time_per_row, 1), "sec/row |",
        "Est. remaining:", 
        if(estimated_time_remaining < 60) paste(round(estimated_time_remaining, 1), "sec") 
        else if(estimated_time_remaining < 3600) paste(round(estimated_time_remaining/60, 1), "min")
        else paste(round(estimated_time_remaining/3600, 1), "hr"), "\n")
  }
  
  # Save progress every 50 rows
  if(i %% 50 == 0 || i == nrow(true_params_matrix)) {
    cat("  Saving progress at row", i, "...\n")
    progress_data <- list(
      true_ses_matrix = true_ses_matrix,
      completed_rows = completed_rows,
      last_completed_row = i,
      timestamp = Sys.time()
    )
    saveRDS(progress_data, progress_filename)
  }
  
  # Print progress every 100 rows
  if(i %% 100 == 0) {
    successful_computations <- sum(!is.na(true_ses_matrix[1:i, "t1_se"]))
    total_elapsed <- as.numeric(Sys.time() - computation_start_time, units = "mins")
    
    cat("=== PROGRESS UPDATE ===\n")
    cat("Completed:", i, "/", nrow(true_params_matrix), "rows (", 
        round(i/nrow(true_params_matrix)*100, 1), "%)\n")
    cat("Successful computations:", successful_computations, "/", i, "(", 
        round(successful_computations/i*100, 1), "%)\n")
    cat("Total elapsed time:", 
        if(total_elapsed < 60) paste(round(total_elapsed, 1), "min")
        else paste(round(total_elapsed/60, 1), "hr"), "\n")
    cat("Estimated total time:", 
        if(i > start_row) {
          total_estimated <- (total_elapsed / (i - start_row + 1)) * nrow(true_params_matrix)
          if(total_estimated < 60) paste(round(total_estimated, 1), "min")
          else paste(round(total_estimated/60, 1), "hr")
        } else "calculating...", "\n")
    cat("========================\n")
  }
}

# Final save
cat("Computation complete! Saving final results...\n")
total_computation_time <- as.numeric(Sys.time() - computation_start_time, units = "mins")
cat("Total computation time:", 
    if(total_computation_time < 60) paste(round(total_computation_time, 1), "minutes")
    else paste(round(total_computation_time/60, 1), "hours"), "\n")

final_data <- list(
  true_ses_matrix = true_ses_matrix,
  true_params_matrix = true_params_matrix,
  completed_rows = completed_rows,
  computation_time = Sys.time(),
  total_rows = nrow(true_params_matrix),
  total_computation_minutes = total_computation_time
)
saveRDS(final_data, final_filename)

# Print final summary
successful_computations <- sum(!is.na(true_ses_matrix[, "t1_se"]))
cat("\n=== FINAL SUMMARY ===\n")
cat("Total parameter combinations:", nrow(true_params_matrix), "\n")
cat("Successful computations:", successful_computations, "\n")
cat("Success rate:", round(successful_computations/nrow(true_params_matrix)*100, 1), "%\n")
cat("Results saved to:", final_filename, "\n")

# Show breakdown by distribution
cat("\nSuccess breakdown by distribution:\n")
for(dist in unique(true_params_matrix$dist)) {
  dist_rows <- which(true_params_matrix$dist == dist)
  dist_successes <- sum(!is.na(true_ses_matrix[dist_rows, "t1_se"]))
  cat("  ", dist, ":", dist_successes, "/", length(dist_rows), 
      "(", round(dist_successes/length(dist_rows)*100, 1), "%)\n")
}

# Clean up progress file if computation completed successfully
if(file.exists(progress_filename) && successful_computations > 0) {
  cat("Removing progress file as computation completed successfully.\n")
  file.remove(progress_filename)
}

cat("\nTrue standard errors matrix computation complete!\n")
cat("Matrix dimensions:", dim(true_ses_matrix), "\n")
cat("Available in variable: true_ses_matrix\n")

  # Apply transformations only if we computed the matrix
  cat("Applying transformations to computed true_ses_matrix...\n")
  true_ses_matrix <- sqrt(true_ses_matrix)
  true_ses_matrix[true_params_matrix[,"dist"]=="GA",] <- true_ses_matrix[true_params_matrix[,"dist"]=="GA",]/sqrt(1000)
  cat("Transformations applied successfully.\n")

} # End of computation block

#true_ses_matrix=sqrt(true_ses_matrix)

# PLOTTING SECTION ######

library(ggplot2)
library(gridExtra)
library(latex2exp)

# Define model renaming function and colors for all plots
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

# Define models to plot and create consistent color mapping
# Include all models that might appear in any of the evaluation plots
# This should match the model_names used to create the matrices
models_to_plot <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")
all_model_labels <- sapply(models_to_plot, rename_model)

# Load scales package for color generation
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(scales)

# Create consistent color mapping for all models that could appear
model_colors <- setNames(scales::hue_pal()(length(all_model_labels)), all_model_labels)

# Debug: Print the color mapping and check against matrix columns
cat("Model color mapping established:\n")
for(i in 1:length(model_colors)) {
  cat("  ", names(model_colors)[i], ":", model_colors[i], "\n")
}

# Check what models are actually available in the matrices
cat("\nModels available in evaluation matrices:\n")
if(exists("vs2_matrix")) {
  cat("VS2 matrix columns:", paste(colnames(vs2_matrix), collapse = ", "), "\n")
  # Show which matrix models have colors defined
  matrix_models <- colnames(vs2_matrix)
  matrix_model_labels <- sapply(matrix_models, rename_model)
  missing_in_colors <- matrix_model_labels[!matrix_model_labels %in% names(model_colors)]
  if(length(missing_in_colors) > 0) {
    cat("WARNING: Matrix models without colors:", paste(missing_in_colors, collapse = ", "), "\n")
    cat("Matrix model codes:", paste(matrix_models[!matrix_model_labels %in% names(model_colors)], collapse = ", "), "\n")
  }
}
if(exists("loglik_matrix")) {
  cat("LogLik matrix columns:", paste(colnames(loglik_matrix), collapse = ", "), "\n")
}

#Run the plotting files
source("Evaluation_plots.R")
source("MU1_MU2_Coefficient_Plots.R")
source("X1_Coefficient_Plots.R")
source("X2_Coefficient_Plots.R")
