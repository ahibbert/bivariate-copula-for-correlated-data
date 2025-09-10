# EXTRACT ALL THE DATA FROM SAVED SIMULATIONS ----

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Set working directory and paths
data_path <- "Data/broad_sims"
cat("Loading broad simulation results from:", data_path, "\n")

# Get list of all RData files in the broad_sims directory
rdata_files <- list.files(data_path, pattern = "\\.RData$", full.names = TRUE)
cat("Found", length(rdata_files), "RData files to process\n")

if (length(rdata_files) == 0) {
  stop("No RData files found in ", data_path, ". Please check the path.")
}

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

# Load true values (old simulations) - to delete ----
numDerivResults <- readRDS(file="Data/numDerivResults_20231127.rds")
se_mles <- readRDS("Data/se_mles_20231127_n100sims20_ALL.rds")

se_sim=se_mles[,c(1,2,3)]/sqrt(10)
se_nd=numDerivResults[,c(1,2,5)]

se_final=se_nd
for (i in 1:3) {
  is.nan(se_final[,i])
  se_final[is.nan(se_final[,i]),i]=se_sim[is.nan(se_final[,i]),i]
}

true_SE_GA=cbind(parameters[parameters$dist=="GA",],se_final)

load(file="Data/lo_mle")
true_SE_LO_in=lo_mle
trueSE_LO=merge(parameters[parameters$dist=="LO",],true_SE_LO_in,by=c("a","b","c","mu1","mu2"))[,c(1:9,13:15)]
colnames(trueSE_LO)[10:12]=c("mu1_se","mu2_se_B2","mu2_se_Bt")
trueSE_LO[,c("mu1_se","mu2_se_B2","mu2_se_Bt")]=sqrt(trueSE_LO[,c("mu1_se","mu2_se_B2","mu2_se_Bt")])
trueSE_NO<-t(rbind((parameters[,"a"]/sqrt(parameters[,"n"]))
                ,(parameters[,"b"]/sqrt(parameters[,"n"]))
                ,sqrt((parameters[,"a"]^2)+(parameters[,"b"]^2)-2*parameters[,"a"]*parameters[,"b"]*parameters[,"c"])/sqrt(parameters[,"n"])))[parameters$dist=="NO",]
colnames(trueSE_NO)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")

trueSE_NO_final=cbind(parameters[parameters$dist=="NO",],trueSE_NO)

load("Data/nb_mle")
trueSE_PO<-matrix(ncol=3,nrow=nrow(parameters[parameters$dist=="PO",]))
#Come back to this
mu1=parameters[,"mu1"]
mu2=parameters[,"mu2"]
a=parameters[,"a"]
b=parameters[,"b"]
c=parameters[,"c"]
n=parameters[,"n"]

  e_x1 = mu1*c*b
  e_x2 = mu2*c*b
  v_x1 = (((mu1^2)*(c*b^2)+(mu1*c*b))/((mu1*c*b)^2))
  v_x2 = (((mu2^2)*(c*b^2)+(mu2*c*b))/((mu2*c*b)^2))
  se_1=sqrt(v_x1)     /sqrt(n)
  se_2=sqrt(v_x2)     /sqrt(n)
  se_3=sqrt(
      (v_x2 + v_x1)
      - log((mu1*mu2*c)/(e_x1*e_x2))
    ) /sqrt(n)

trueSE_PO[,3]=sqrt(merge(parameters,nb_mle,by=c("a","b","c","mu1","mu2"))[,"var_Bt"])

trueSE_PO[,1]=se_1[parameters$dist=="PO"]
trueSE_PO[,2]=se_2[parameters$dist=="PO"]

colnames(trueSE_PO)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")

trueSE_PO_final=cbind(parameters[parameters$dist=="PO",],trueSE_PO)
true_SE_out=rbind(true_SE_GA,trueSE_NO_final,trueSE_PO_final,trueSE_LO)

# Load true values from file if available ----
# Load newly computed true SEs (if available) to override the old true_SE_out
# This uses the get_true_ses function results computed at the end of this script
tryCatch({
  cat("Attempting to load computed true SEs...\n")
  
  # Try to find the most recent true_ses_complete file
  cache_files <- list.files("Cache", pattern = "true_ses_complete_.*\\.rds$", full.names = TRUE)
  if(length(cache_files) > 0) {
    # Use the most recent file
    cache_file <- cache_files[which.max(file.mtime(cache_files))]
    cat("Loading true SEs from:", basename(cache_file), "\n")
    
    final_data <- readRDS(cache_file)
    true_ses_matrix <- final_data$true_ses_matrix
    cat("Successfully loaded true_ses_matrix with dimensions:", dim(true_ses_matrix), "\n")
    cat("Available SE columns:", colnames(true_ses_matrix), "\n")
    cat("Non-missing t1_se values:", sum(!is.na(true_ses_matrix[, "t1_se"])), "\n")
  } else {
    cat("No true_ses_complete files found in Cache directory\n")
    true_ses_matrix <- NULL
  }
}, error = function(e) {
  cat("Could not load computed true SEs:", e$message, "\n")
  cat("Will use old true_SE_out method instead\n")
  true_ses_matrix <- NULL
})

# Compute true values (new simulations) ----

cat("\n=== Computing True Standard Errors for All Parameter Combinations ===\n")

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

# Plotting ####

# NOTE: This plotting section can now use either:
# 1. NEW: true_ses_matrix (computed by get_true_ses function at end of script) - PREFERRED
# 2. OLD: true_SE_out (pre-computed theoretical SEs) - FALLBACK
# The code will automatically detect which is available and use the appropriate source

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
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Define model renaming function (moved here to be available for all plots)
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

models_to_plot <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")

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

# Create a function to prepare plot data for a specific distribution
prepare_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu1_coef_matrix))
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

  #true_mu1 <- as.numeric(true_params_matrix$mu1[dist_rows])
  
  # Get correlation values for this distribution (from skew_matrix column 2)
  corr_vals <- skew_matrix[dist_rows, 2]
  
  # Get estimated mu1 values for this distribution and selected models
  if(dist_name %in% c("GA")) {
    est_mu1_matrix <- exp(mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE])
  } else if (dist_name == "PO") {
    est_mu1_matrix <- exp(mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE])
  } else if(dist_name == "LO") {
    est_mu1_matrix <- exp(mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE])/(1+exp(mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE]))
  } else {
    est_mu1_matrix <- mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE]
  }

  #est_mu1_matrix <- mu1_coef_matrix[dist_rows, models_to_plot, drop = FALSE]
  
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
            if(ratio > 10 || ratio < 0.1) { # Using 10x and 0.1x threshold
              exclusion_mask[i, j] <- TRUE
            }
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
    # Use consistent color mapping, only showing models present in this plot
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    coord_cartesian( ylim=c(-1,1)) +
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
  
  return(p)
}

# Create a function to prepare SE plot data for a specific distribution
prepare_se_plot_data <- function(dist_name, models_to_plot = NULL) {
  
  # If no models specified, use all available models with non-missing data
  if(is.null(models_to_plot)) {
    # Find models with at least some non-missing data
    non_missing_counts <- colSums(!is.na(mu1_se_matrix))
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
  est_se_matrix <- mu1_se_matrix[dist_rows, models_to_plot, drop = FALSE]
  
  # OLD CODE: Get theoretical SE values from true_SE_out by matching parameter sets
  # true_se_vals <- rep(NA, length(dist_rows))
  # 
  # # Match each row in true_params_matrix with corresponding row in true_SE_out
  # for(i in 1:length(dist_rows)) {
  #   row_idx <- dist_rows[i]
  #   params <- true_params_matrix[row_idx, ]
  #   
  #   # Find matching row in true_SE_out
  #   matching_rows <- which(
  #     true_SE_out$dist == params$dist &
  #     true_SE_out$mu1 == as.numeric(params$mu1) &
  #     true_SE_out$mu2 == as.numeric(params$mu2) &
  #     as.numeric(true_SE_out$n) == as.numeric(params$n)
  #   )
  #   
  #   # Handle a, b, c parameters which can be NA
  #   if(length(matching_rows) > 0) {
  #     for(param in c("a", "b", "c")) {
  #       if(!is.na(params[[param]])) {
  #         matching_rows <- matching_rows[
  #           abs(as.numeric(true_SE_out[[param]][matching_rows]) - as.numeric(params[[param]])) < 1e-10
  #         ]
  #       } else {
  #         matching_rows <- matching_rows[is.na(true_SE_out[[param]][matching_rows])]
  #       }
  #     }
  #   }
  #   
  #   if(length(matching_rows) > 0) {
  #     true_se_vals[i] <- true_SE_out$mu1_se[matching_rows[1]]
  #   }
  # }
  
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
    }
    
    # Use t1_se from the actual matrix for mu1 theoretical SEs with proper indexing
    if(is.matrix(actual_ses_matrix) && "t1_se" %in% colnames(actual_ses_matrix)) {
      true_se_vals <- actual_ses_matrix[dist_rows, "t1_se"]
      cat("Using computed true SEs from true_ses_matrix for", dist_name, "distribution\n")
    } else if(is.data.frame(actual_ses_matrix) && "t1_se" %in% colnames(actual_ses_matrix)) {
      true_se_vals <- actual_ses_matrix$t1_se[dist_rows]
      cat("Using computed true SEs from data.frame for", dist_name, "distribution\n")
    } else {
      cat("Warning: t1_se column not found in true_ses_matrix for", dist_name, "\n")
      true_se_vals <- rep(NA, length(dist_rows))
    }
  } else {
    # Fallback to old method if true_ses_matrix not available
    cat("Warning: true_ses_matrix not available, using old true_SE_out method for", dist_name, "\n")
    true_se_vals <- rep(NA, length(dist_rows))
    
    # Match each row in true_params_matrix with corresponding row in true_SE_out
    for(i in 1:length(dist_rows)) {
      row_idx <- dist_rows[i]
      params <- true_params_matrix[row_idx, ]
      
      # Find matching row in true_SE_out
      matching_rows <- which(
        true_SE_out$dist == params$dist &
        true_SE_out$mu1 == as.numeric(params$mu1) &
        true_SE_out$mu2 == as.numeric(params$mu2) &
        as.numeric(true_SE_out$n) == as.numeric(params$n)
      )
      
      # Handle a, b, c parameters which can be NA
      if(length(matching_rows) > 0) {
        for(param in c("a", "b", "c")) {
          if(!is.na(params[[param]])) {
            matching_rows <- matching_rows[
              abs(as.numeric(true_SE_out[[param]][matching_rows]) - as.numeric(params[[param]])) < 1e-10
            ]
          } else {
            matching_rows <- matching_rows[is.na(true_SE_out[[param]][matching_rows])]
          }
        }
      }
      
      if(length(matching_rows) > 0) {
        true_se_vals[i] <- true_SE_out$mu1_se[matching_rows[1]]
      }
    }
  }
  
  # Convert estimated SEs to long format for ggplot
  se_plot_data <- data.frame(
    corr_vals = rep(corr_vals, length(models_to_plot)),
    se_value = as.vector(est_se_matrix),
    model = rep(models_to_plot, each = length(corr_vals)),
    model_label = rep(sapply(models_to_plot, rename_model), each = length(corr_vals)),
    dist = dist_name,
    type = "Estimated"
  )
  
  # Add theoretical SE data
  true_se_data <- data.frame(
    corr_vals = corr_vals,
    se_value = true_se_vals,
    model = "true",
    model_label = "Theoretical",
    dist = dist_name,
    type = "Theoretical"
  )
  
  # Combine estimated and theoretical data
  combined_se_data <- rbind(se_plot_data, true_se_data)
  
  # Remove rows with missing values
  combined_se_data <- combined_se_data[complete.cases(combined_se_data), ]
  
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
    coord_cartesian(ylim = c(0, 0.4)) +
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

# Create a plot with legend to extract it
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
  
  # Extract legend from one of the plots
  legend_plot_temp <- create_dist_plot("NO", models_to_plot, show_legend = TRUE) + 
    theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.title = element_text(size = 11))
  
  # Function to extract legend from ggplot
  get_legend <- function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  # Extract the legend
  shared_legend <- get_legend(legend_plot_temp)
  
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
  if(model %in% colnames(mu1_coef_matrix)) {
    non_missing <- sum(!is.na(mu1_coef_matrix[, model]))
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

# PLOTTING SECTION 

## Evaluation Score Plots

cat("\n=== Creating Evaluation Score Plots ===\n")

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
  
  # Remove rows with missing values
  plot_data_before <- nrow(plot_data)
  plot_data <- plot_data[complete.cases(plot_data), ]
  cat("Removed", plot_data_before - nrow(plot_data), "rows with missing values, keeping", nrow(plot_data), "\n")
  
  # Set factor levels for model_label to maintain desired order
  if(nrow(plot_data) > 0) {
    model_label_order <- sapply(models_to_plot, rename_model)
    plot_data$model_label <- factor(plot_data$model_label, levels = model_label_order)
  }
  
  return(plot_data)
}

# Create evaluation score plots for each distribution and metric
create_eval_plot <- function(dist_name, metric_name, metric_matrix, models_to_plot = NULL, show_legend = TRUE) {
  
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
  
  # Create the plot with consistent color mapping
  p <- ggplot(plot_data, aes(x = corr_vals, y = metric_value, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    # Use consistent color mapping, only showing models present in this plot
    scale_color_manual(values = model_colors[unique_models], 
                       breaks = unique_models,
                       limits = unique_models) +
    coord_cartesian(xlim = c(0.2, 0.8)) +
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
  
  return(p)
}

# Create evaluation score plots for all distributions and metrics
cat("Creating evaluation score plots for all distributions and metrics...\n")

# Define the desired order for models (same as used in other functions)
desired_order <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm", "cop_n", "cop")

# Define the models to plot in the desired order (use available models from score matrices)
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

# Initialize list to store all evaluation plots
eval_plots <- list()
plot_counter <- 1

for(dist in distributions) {
  for(metric_name in names(metrics)) {
    metric_matrix <- metrics[[metric_name]]
    
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
    
    # TEMPORARY: Exclude GAMM from VS2 and VS2_Weighted plots
    # To revert: comment out or delete the following 3 lines
    #if(metric_name %in% c("VS2", "VS2_Weighted")) {
    #  current_models <- current_models[current_models != "gamm"]
    #}
    
    cat("Creating plot for:", dist, metric_name, "with models:", paste(current_models, collapse = ", "), "\n")
    
    # Create plot without individual legend
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
legend_plot_eval <- create_eval_plot("NO", "VS2", vs2_matrix, models_to_plot_eval, show_legend = TRUE)

# Function to extract legend from ggplot (needed for evaluation plots)
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

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
    # Extract legend from evaluation plot with error handling
    tryCatch({
      legend_plot_eval_temp <- create_eval_plot("NO", "VS2", vs2_matrix, models_to_plot_eval, show_legend = TRUE)
      
      if(!is.null(legend_plot_eval_temp)) {
        legend_plot_eval_temp <- legend_plot_eval_temp + 
          theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_text(size = 9))
        
        # Extract the legend
        shared_legend_eval <- get_legend(legend_plot_eval_temp)
        
        if(is.null(shared_legend_eval)) {
          cat("Warning: Failed to extract legend, proceeding without shared legend\n")
          shared_legend_eval <- grid::nullGrob()  # Create empty grob
        }
      } else {
        cat("Warning: Failed to create legend plot, proceeding without shared legend\n")
        shared_legend_eval <- grid::nullGrob()  # Create empty grob
      }
    }, error = function(e) {
      cat("Error creating legend:", e$message, "\n")
      shared_legend_eval <- grid::nullGrob()  # Create empty grob
    })
    
    # Arrange plots in a 4x4 grid (4 distributions x 4 metrics) with error handling
    tryCatch({
      eval_plots_grid <- do.call(arrangeGrob, c(eval_plots, ncol = 4))
      
      # Combine the plots grid with the shared legend at the bottom
      combined_eval_plot <- arrangeGrob(eval_plots_grid, shared_legend_eval, 
                                       ncol = 1, 
                                       heights = c(10, 1))
      
      # Display the combined evaluation plot
      grid.draw(combined_eval_plot)
      
      # Save the combined evaluation plot to Charts directory
      eval_chart_filename <- paste0("Charts/Evaluation_Scores_vs_Correlation_", Sys.Date(), ".png")
      
      ggsave(eval_chart_filename, plot = combined_eval_plot, width = 16, height = 12, dpi = 600)
      cat("Combined evaluation scores plot saved to:", eval_chart_filename, "\n")
      
      cat("Evaluation score plots created successfully!\n")
      cat("Number of evaluation plots created:", length(eval_plots), "\n")
      cat("Metrics plotted:", paste(names(metrics), collapse = ", "), "\n")
      cat("Distributions plotted:", paste(distributions, collapse = ", "), "\n")
      
    }, error = function(e) {
      cat("Error arranging evaluation plots:", e$message, "\n")
      cat("Attempting to save individual plots instead...\n")
      
      # Save individual plots if grid arrangement fails
      for(i in 1:length(eval_plots)) {
        if(!is.null(eval_plots[[i]])) {
          plot_name <- names(eval_plots)[i]
          if(is.null(plot_name) || plot_name == "") plot_name <- paste0("eval_plot_", i)
          
          individual_filename <- paste0("Charts/", plot_name, "_", Sys.Date(), ".png")
          ggsave(individual_filename, plot = eval_plots[[i]], width = 8, height = 6, dpi = 300)
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
for(metric_name in names(metrics)) {
  metric_matrix <- metrics[[metric_name]]
  non_missing <- sum(!is.na(metric_matrix))
  total <- length(metric_matrix)
  cat(sprintf("%s: %d/%d non-missing values (%.1f%%)\n", 
              metric_name, non_missing, total, 
              if(total > 0) non_missing / total * 100 else 0))
}

# Note: true_ses_matrix loading has been moved to earlier in the script (after true_SE_out creation)
# The plotting functions will automatically use true_ses_matrix if available, otherwise fall back to true_SE_out
