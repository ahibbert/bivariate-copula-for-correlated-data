#### EXTRACT ALL THE DATA FROM SAVED SIMULATIONS ####

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

#Extract theoretical valus for errors 
load(file="Data/numDerivResults_20231127.rds")
load("Data/se_mles_20231127_n100sims20_ALL.rds") #se_mles

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

########### OK NOW PLOTTING TIME #############

# Create plots for estimation bias vs correlation for each distribution
library(ggplot2)
library(gridExtra)
library(latex2exp)

# Define model renaming function (copied from TEMP Plotting file)
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
  
  # Create the plot with default ggplot2 colors (like TEMP file)
  p <- ggplot(plot_data, aes(x = corr_vals, y = bias, color = model_label)) +
    geom_smooth(method = "loess", se = TRUE, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    coord_cartesian(ylim=c(-1,1)) +
    labs(
      title = paste("Distribution:", dist_name),
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

# Define models to include in the plots (all models except unwanted GJRM variants)
# Keep all non-GJRM models and only GJRM (C) and GJRM (N)
all_available_models <- colnames(mu1_coef_matrix)
gjrm_models_to_exclude <- c("cop_j", "cop_g", "cop_f", "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t")
models_to_plot <- all_available_models[!all_available_models %in% gjrm_models_to_exclude]

# Filter to only models with some non-missing data
non_missing_counts <- colSums(!is.na(mu1_coef_matrix[, models_to_plot, drop = FALSE]))
models_to_plot <- models_to_plot[non_missing_counts > 0]

cat("Models to be plotted:", paste(models_to_plot, collapse = ", "), "\n")
cat("Model labels:", paste(sapply(models_to_plot, rename_model), collapse = ", "), "\n")

# Create plots for each distribution
cat("Creating mu1 bias vs correlation plots for each distribution...\n")

# Create plots without individual legends (we'll add a shared legend)
plot_NO <- create_dist_plot("NO", models_to_plot, show_legend = FALSE)
plot_PO <- create_dist_plot("PO", models_to_plot, show_legend = FALSE) 
plot_GA <- create_dist_plot("GA", models_to_plot, show_legend = FALSE)
plot_LO <- create_dist_plot("LO", models_to_plot, show_legend = FALSE)

# Create a plot with legend to extract it
legend_plot <- create_dist_plot("NO", models_to_plot, show_legend = TRUE)

# Remove NULL plots
plots_list <- list(plot_NO, plot_PO, plot_GA, plot_LO)
plots_list <- plots_list[!sapply(plots_list, is.null)]

if(length(plots_list) > 0) {
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
  
  # Create plots without legends
  plots_no_legend <- plots_list
  
  # Remove NULL plots from the no-legend list
  plots_no_legend <- plots_no_legend[!sapply(plots_no_legend, is.null)]
  
  # Arrange plots in a 2x2 grid without legends
  plots_grid <- do.call(arrangeGrob, c(plots_no_legend, ncol = 2))
  
  # Combine the plots grid with the shared legend at the bottom
  combined_plot <- arrangeGrob(plots_grid, shared_legend, 
                              ncol = 1, 
                              heights = c(10, 1))
  
  # Display the combined plot
  grid.draw(combined_plot)
  
  # Save the combined plot to Charts directory
  chart_filename <- paste0("Charts/Mu1_Bias_vs_Correlation_Combined_", Sys.Date(), ".png")
  
  ggsave(chart_filename, plot = combined_plot, width = 12, height = 11, dpi = 300)
  cat("Combined plot saved to:", chart_filename, "\n")
  
  cat("Mu1 bias vs correlation plots created successfully!\n")
  cat("Number of plots created:", length(plots_list), "\n")
  
  # Also create individual plots for closer inspection
  if(!is.null(plot_NO)) {
    cat("NO distribution plot available\n")
  }
  if(!is.null(plot_PO)) {
    cat("PO distribution plot available\n") 
  }
  if(!is.null(plot_GA)) {
    cat("GA distribution plot available\n")
  }
  if(!is.null(plot_LO)) {
    cat("LO distribution plot available\n")
  }
  
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

