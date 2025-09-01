# Load required libraries and source dependencies
library(callr)
library(parallel)
library(future)
library(future.apply)
source("common_functions.R")

# Enhanced parallel processing setup with resource management
# Detect available cores (leave more free for system stability)
n_cores <- max(1, min(detectCores() - 2, 8))  # Cap at 6 cores for stability
cat(sprintf("Setting up parallel processing with %d cores\n", n_cores))

# Set up future plan with resource limits
plan(multisession, 
     workers = n_cores,
     # Set memory limits and timeouts
     rscript_args = c("--max-ppsize=500000000")  # 500MB per process
)

# Configure future options for better error handling
options(future.globals.maxSize = 1024 * 1024^2)  # 1GB max for globals
options(future.rng.onMisuse = "ignore")

# Define the list of input parameter sets as a list of lists
input_params_list <- list(
  list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  list(dist="NO", a=1, b=1, c=0.5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  list(dist="LO",a=NA,b=NA,c=.25,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
  list(dist="LO",a=NA,b=NA,c=.5,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
  list(dist="GA", a=.5, b=2, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  list(dist="GA", a=.75, b=1.5, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  list(dist="PO", a=NA, b=.2, c=5, mu1=1, mu2=1, x1=1, x2=0.01, n=1000), #
  #list(dist="PO", a=NA, b=5, c=5, mu1=5, mu2=5, x1=1, x2=0.01, n=1000),
  #list(dist="PO", a=NA, b=2, c=2, mu1=2, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="PO", a=NA, b=1, c=2, mu1=1, mu2=1, x1=1, x2=0.01, n=1000),
  #list(dist="PO", a=NA, b=.2, c=1, mu1=1, mu2=1, x1=1, x2=0.01, n=1000),
  list(dist="PO", a=NA, b=1, c=1, mu1=1, mu2=1, x1=1, x2=0.01, n=1000)
  #list(dist="PO", a=NA, b=5, c=.2, mu1=1, mu2=1, x1=1, x2=0.01, n=1000),

  # All parameters different
  #list(dist="PO", a=NA, b=.2, c=1, mu1=2, mu2=5, x1=1, x2=0.01, n=1000),
  #list(dist="PO", a=NA, b=5, c=2, mu1=1, mu2=.2, x1=1, x2=0.01, n=1000),

  # Reverse pairing
  #list(dist="PO", a=NA, b=1, c=.2, mu1=2, mu2=5, x1=1, x2=0.01, n=1000),

  # Edge cases and more extreme values
  #list(dist="PO", a=NA, b=.2, c=5, mu1=5, mu2=.2, x1=1, x2=0.01, n=1000),
  #list(dist="PO", a=NA, b=5, c=.2, mu1=.2, mu2=5, x1=1, x2=0.01, n=1000)
)

outer_sims <- 100  # Number of simulations per input set

# Enhanced function to run a single simulation with better error handling
run_single_simulation <- function(params, sim_index, param_set_idx) {
  max_retries <- 5  # Reduced retries to avoid excessive resource usage
  attempt <- 1
  success <- FALSE
  
  cat(sprintf("Input set %d, simulation %d\n", param_set_idx, sim_index))
  
  while (attempt <= max_retries && !success) {
    eval <- tryCatch({
      # Create deterministic seed based on parameters and simulation index
      # This ensures reproducibility: same params + sim_index = same seed
      param_hash <- sum(c(
        match(params$dist, c("NO", "PO", "GA", "LO")) * 1000,
        ifelse(is.na(params$a), 0, params$a * 100),
        ifelse(is.na(params$b), 0, params$b * 10),
        ifelse(is.na(params$c), 0, params$c * 1),
        params$mu1 * 10000,
        params$mu2 * 1000,
        params$x1 * 100000,
        params$x2 * 1000000,
        params$n
      ))
      
      # Deterministic seed: parameter hash + param set + simulation index + attempt
      deterministic_seed <- (param_hash %% 100000) + param_set_idx * 10000 + sim_index * 100 + attempt
      
      # Use callr with built-in timeout
      callr::r(
        func = function(dist, a, b, c, mu1, mu2, x1, x2, n, sim_seed) {
          # Set deterministic seed for reproducibility
          set.seed(sim_seed)
          
          source("common_functions.R")
          dataset = generateBivDist_withCov(n=n, a=a, b=b, c=c, mu1=mu1, mu2=mu2, dist=dist, x1=x1, x2=x2)
          
          # Add some debugging info to verify different datasets (simple checksum)
          if (exists("dataset") && !is.null(dataset$Y1)) {
            dataset_sum <- sum(dataset$Y1, na.rm = TRUE)
            cat("Sim", sim_seed, "- Y1 sum:", round(dataset_sum, 4), "- First Y1:", round(dataset$Y1[1], 4), "\n")
          }
          
          # Set dataset as global as old S3 functions in gamlss and sometimes GJRM break without it
          assign("dataset", dataset, envir = .GlobalEnv)
          fits = fitBivModels_Bt_withCov(dataset=dataset, dist=dist, include="ALL",
                                         a=a, b=b, c=c, mu1=mu1, mu2=mu2,
                                         calc_actuals=FALSE, cv=FALSE)
          eval = evaluateModels(fits, vg_sims=100)
          return(eval)
        },
        args = c(params, list(sim_seed = deterministic_seed)),
        # Add process-level resource limits and timeout
        cmdargs = c("--max-ppsize=500000000"),  # 500MB per subprocess
        timeout = 300  # 5 minute timeout per simulation
      )
      
    }, error = function(e) {
      # Enhanced error reporting with better classification
      error_msg <- paste0("Error in simulation attempt ", attempt, " (set ", param_set_idx, ", sim ", sim_index, "):\n")
      
      # Check for specific error types
      if (grepl("FutureInterruptError|InterruptedError|interrupt", e$message, ignore.case = TRUE)) {
        error_msg <- paste0(error_msg, "INTERRUPTION ERROR - likely resource exhaustion\n")
      } else if (grepl("memory|cannot allocate", e$message, ignore.case = TRUE)) {
        error_msg <- paste0(error_msg, "MEMORY ERROR - insufficient memory\n")
      } else if (grepl("timeout", e$message, ignore.case = TRUE)) {
        error_msg <- paste0(error_msg, "TIMEOUT ERROR - process took too long\n")
      }
      
      error_msg <- paste0(error_msg, "Message: ", e$message, "\n")
      
      # Add call information if available
      if (!is.null(e$call)) {
        error_msg <- paste0(error_msg, "Call: ", deparse(e$call), "\n")
      }
      
      # Add subprocess error details for callr errors
      if ("parent" %in% names(e) && !is.null(e$parent)) {
        error_msg <- paste0(error_msg, "Subprocess error: ", e$parent$message, "\n")
      }
      
      cat(error_msg)
      
      # Force garbage collection and small delay for memory recovery
      gc()
      Sys.sleep(1)
      
      NULL
    })
    
    if (!is.null(eval)) {
      success <- TRUE
    } else {
      attempt <- attempt + 1
      if (attempt <= max_retries) {
        cat(sprintf("Retrying simulation %d (attempt %d/%d)...\n", sim_index, attempt, max_retries))
        # Progressive delay between retries
        Sys.sleep(attempt * 2)
      }
    }
  }
  
  if (!success) {
    cat(sprintf("All attempts failed for simulation %d in parameter set %d.\n", sim_index, param_set_idx))
    return(NA)
  } else {
    return(eval)
  }
}

# Enhanced function to process a single parameter set with chunking
process_parameter_set <- function(param_set_idx, params, outer_sims) {
  # Check if output file already exists (cache check)
  base_filename_pattern <- paste0(
    "Data/CoefSimData_", paste(
      params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims,
      sep="_"
    ), "_*.RData"
  )
  
  existing_files <- Sys.glob(base_filename_pattern)
  
  if (length(existing_files) > 0) {
    cat(sprintf("Input set %d: Output file already exists, skipping: %s\n", param_set_idx, existing_files[1]))
    return(NULL)
  }
  
  # Define the save filename for new file
  save_filename <- paste0(
    "Data/CoefSimData_", paste(
      params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims, Sys.Date(),
      sep="_"
    ), ".RData"
  )
  
  cat(sprintf("Input set %d: Processing parameters: %s\n", param_set_idx, paste(unlist(params), collapse=", ")))
  
  # Process simulations in smaller chunks to reduce memory pressure
  chunk_size <- min(20, outer_sims)  # Process 20 simulations at a time
  num_chunks <- ceiling(outer_sims / chunk_size)
  
  cat(sprintf("Running %d simulations in %d chunks of size %d for parameter set %d...\n", 
              outer_sims, num_chunks, chunk_size, param_set_idx))
  
  eval_outer <- list()
  
  for (chunk_idx in 1:num_chunks) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, outer_sims)
    
    cat(sprintf("Processing chunk %d/%d (simulations %d-%d)...\n", 
                chunk_idx, num_chunks, start_idx, end_idx))
    
    # Run chunk with retry mechanism for the entire chunk
    chunk_success <- FALSE
    chunk_retries <- 0
    max_chunk_retries <- 3
    
    while (!chunk_success && chunk_retries < max_chunk_retries) {
      tryCatch({
        # Generate deterministic seeds for each simulation in this chunk
        # Based on parameter set and simulation indices for reproducibility
        param_hash <- sum(c(
          match(params$dist, c("NO", "PO", "GA", "LO")) * 1000,
          ifelse(is.na(params$a), 0, params$a * 100),
          ifelse(is.na(params$b), 0, params$b * 10),
          ifelse(is.na(params$c), 0, params$c * 1),
          params$mu1 * 10000,
          params$mu2 * 1000,
          params$x1 * 100000,
          params$x2 * 1000000,
          params$n
        ))
        
        chunk_results <- future_lapply(seq_along(start_idx:end_idx), function(chunk_i) {
          actual_sim_idx <- start_idx + chunk_i - 1
          run_single_simulation(params, actual_sim_idx, param_set_idx)
        }, future.seed = NULL)  # Let our deterministic seeding handle this
        
        # Store results
        for (i in seq_along(chunk_results)) {
          eval_outer[[start_idx + i - 1]] <- chunk_results[[i]]
        }
        
        chunk_success <- TRUE
        
      }, error = function(e) {
        chunk_retries <<- chunk_retries + 1
        cat(sprintf("Chunk %d failed (attempt %d/%d): %s\n", 
                    chunk_idx, chunk_retries, max_chunk_retries, e$message))
        
        if (chunk_retries < max_chunk_retries) {
          cat("Forcing garbage collection and retrying chunk...\n")
          gc()
          Sys.sleep(5)  # Wait before retry
        }
      })
    }
    
    if (!chunk_success) {
      cat(sprintf("Chunk %d failed after %d attempts. Filling with NA values.\n", 
                  chunk_idx, max_chunk_retries))
      for (i in start_idx:end_idx) {
        eval_outer[[i]] <- NA
      }
    }
    
    # Force garbage collection between chunks
    gc()
    Sys.sleep(2)  # Brief pause between chunks
  }
  
  # Process results (same as original code but with safer indexing)
  score_items <- list()
  score_item_names <- c("vs2_wt", "vs2", "es", "vs1", "vs2_wt_coronly", "logliks")
  
  # Also extract df values separately for BIC calculation
  df_values <- list()
  
  successful_sims <- 0
  for (i in seq_along(eval_outer)) {
    if (!is.na(eval_outer[[i]])[1]) {  # Check if simulation was successful
      successful_sims <- successful_sims + 1
      for (item in score_item_names) {
        new_item <- if (item == "logliks") {
          eval_outer[[i]][[item]][, 1]  # Extract loglik values (column 1)
        } else {
          eval_outer[[i]][[item]]
        }
        if (successful_sims == 1) {
          score_items[[item]] <- new_item
        } else {
          score_items[[item]] <- rbind(score_items[[item]], new_item)
        }
      }
      
      # Extract df values (column 2) from logliks for BIC calculation
      df_item <- eval_outer[[i]][["logliks"]][, 2]
      if (successful_sims == 1) {
        df_values[["logliks"]] <- df_item
      } else {
        df_values[["logliks"]] <- rbind(df_values[["logliks"]], df_item)
      }
    }
  }
  
  cat(sprintf("Successful simulations: %d/%d\n", successful_sims, outer_sims))
  # Process parameter estimates
  par_estimates <- list()
  par_item_names <- c("coefficients", "ses", "sigmas", "correlations")
  for (i in 1:11) {
    par_estimates[[i]] <- matrix(nrow=length(eval_outer), ncol=16)
  }
  
  for (i in seq_along(eval_outer)) {
    if (!is.na(eval_outer[[i]])[1]) {  # Check if simulation was successful
      z <- 1
      for (item in par_item_names) {
        new_item <- data.frame(eval_outer[[i]][[item]])
        for (j in seq_len(ncol(new_item))) {
          par_estimates[[z]][i, ] <- new_item[, j]
          z <- z + 1
        }
      }
    }
  }
  
  names(par_estimates) <- c("mu1", "mu2", "x1", "x2", "mu1_se", "mu2_se", "x1_se", "x2_se", "sigma1", "sigma2", "corr")
  for (item in names(par_estimates)) {
    colnames(par_estimates[[item]]) <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
                                         "cop", "cop_n", "cop_j", "cop_g", "cop_f",
                                         "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t")
  }
  
  # Process timing and convergence data
  times <- matrix(NA, nrow=length(eval_outer), ncol=length(colnames(par_estimates[[1]])))
  conv <- matrix(NA, nrow=length(eval_outer), ncol=length(colnames(par_estimates[[1]])))
  colnames(times) <- colnames(par_estimates[[1]])
  colnames(conv) <- colnames(par_estimates[[1]])
  
  for (i in seq_along(eval_outer)) {
    if (!is.na(eval_outer[[i]])[1]) {  # Check if simulation was successful
      time2 <- eval_outer[[i]][["timer"]][,2]
      time1 <- eval_outer[[i]][["timer"]][,1]
      timediff <- difftime(time2, time1, units = "secs")
      times[i, ] <- c(timediff)
      conv[i, ] <- eval_outer[[i]][["conv"]]
    }
  }
  
  # Save the results for this parameter set
  save(list = c("par_estimates", "score_items", "df_values", "times", "conv"), file = save_filename)
  cat(sprintf("Parameter set %d completed and saved to: %s\n", param_set_idx, save_filename))
  
  return(save_filename)
}

# Main execution with progress tracking and recovery
cat("Starting robust parallel processing of parameter sets...\n")
cat(sprintf("Total parameter sets: %d\n", length(input_params_list)))
cat(sprintf("Simulations per set: %d\n", outer_sims))
cat(sprintf("Using %d cores for parallel processing\n", n_cores))

# Track progress and allow resuming
completed_sets <- c()
failed_sets <- c()

for (param_set_idx in seq_along(input_params_list)) {
  cat(sprintf("\n=== Processing Parameter Set %d/%d ===\n", param_set_idx, length(input_params_list)))
  
  params <- input_params_list[[param_set_idx]]
  
  tryCatch({
    result <- process_parameter_set(param_set_idx, params, outer_sims)
    completed_sets <- c(completed_sets, param_set_idx)
    cat(sprintf("✓ Parameter set %d completed successfully\n", param_set_idx))
    
  }, error = function(e) {
    failed_sets <- c(failed_sets, param_set_idx)
    cat(sprintf("✗ Parameter set %d failed: %s\n", param_set_idx, e$message))
    
    # Log the failure
    error_log <- paste0("Error in parameter set ", param_set_idx, " at ", Sys.time(), ":\n", e$message, "\n\n")
    write(error_log, file = "parameter_set_errors.log", append = TRUE)
  })
  
  # Progress summary
  cat(sprintf("Progress: %d completed, %d failed, %d remaining\n", 
              length(completed_sets), length(failed_sets), 
              length(input_params_list) - length(completed_sets) - length(failed_sets)))
}

# Final summary
cat("\n=== PROCESSING COMPLETE ===\n")
cat(sprintf("Completed parameter sets: %d\n", length(completed_sets)))
cat(sprintf("Failed parameter sets: %d\n", length(failed_sets)))

if (length(failed_sets) > 0) {
  cat("Failed parameter set indices:", paste(failed_sets, collapse = ", "), "\n")
  cat("Check parameter_set_errors.log for detailed error information.\n")
}

cat("All processing finished!\n")
