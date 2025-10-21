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

#For the broad mean simulations
load("input_params_list.RData")
outer_sims <- 1  # Number of simulations per input set
mode="broad_mean"
bt_mode=TRUE

# Testing for a single run
#params <- input_params_list[[1]]
#test_result <- run_single_simulation(params, 1, 1, bt_mode=TRUE)
#str(test_result)

# Define the list of input parameter sets as a list of lists
#input_params_list <- list(
#  list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
#  list(dist="NO", a=1, b=1, c=0.5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
#  list(dist="LO",a=NA,b=NA,c=.25,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
#  list(dist="LO",a=NA,b=NA,c=.5,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
#  list(dist="GA", a=.5, b=2, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
#  list(dist="GA", a=.75, b=1.5, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
#  list(dist="PO", a=NA, b=.2, c=5, mu1=1, mu2=1, x1=1, x2=0.01, n=1000),
#  list(dist="PO", a=NA, b=1, c=1, mu1=1, mu2=1, x1=1, x2=0.01, n=1000)
#)
#outer_sims <- 100  # Number of simulations per input set
#mode="specific_coef"
#bt_mode=TRUE

# Determine parallelization strategy based on outer_sims
if (outer_sims == 1) {
  cat("Using parameter-set-level parallelization (parallelizing multiple parameter sets)\n")
  parallel_strategy <- "parameter_sets"
  # Use fewer cores per parameter set, more parameter sets in parallel
  param_cores <- max(1, min(detectCores() - 1, 10))  # Cores for individual parameter sets
  param_sets_parallel <- max(1, min(n_cores, length(input_params_list)))
} else {
  cat("Using simulation-level parallelization (parallelizing multiple simulations within each parameter set)\n")
  parallel_strategy <- "simulations"
  param_cores <- n_cores
  param_sets_parallel <- 1
}

# Uncomment the next line to test with a smaller subset first
# input_params_list <- input_params_list[1:4]  # Test with first 4 parameter sets only

# Enhanced function to run a single simulation with better error handling
run_single_simulation <- function(params, sim_index, param_set_idx,bt_mode) {
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
        ifelse(is.na(params$a), 0, as.numeric(params$a) * 100),
        ifelse(is.na(params$b), 0, as.numeric(params$b) * 10),
        ifelse(is.na(params$c), 0, as.numeric(params$c) * 1),
        as.numeric(params$mu1) * 10000,
        as.numeric(params$mu2) * 1000,
        as.numeric(params$x1) * 100000,
        as.numeric(params$x2) * 1000000,
        as.numeric(params$n)
      ), na.rm = TRUE)
      
      # Deterministic seed: parameter hash + param set + simulation index + attempt
      # This ensures reproducibility: same inputs = same seed = same results
      # Use modulo and as.integer to ensure valid seed range
      deterministic_seed <- as.integer((param_hash %% 100000) + param_set_idx * 10000 + sim_index * 100 + attempt)
      
      # Ensure seed is within valid range (R requires seeds between 1 and 2^31-1)
      deterministic_seed <- max(1, min(deterministic_seed, 2147483647))
      
      # Use callr with built-in timeout
      callr::r(
        func = function(dist, a, b, c, mu1, mu2, x1, x2, n, sim_seed, bt_mode) {
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
          
          # Debug: Check bt_mode value in subprocess
          cat("DEBUG - bt_mode in subprocess:", bt_mode, "\n")
          
          fits = fitBivModels_Bt_withCov(dataset=dataset, dist=dist, include="ALL",
                                         a=a, b=b, c=c, mu1=mu1, mu2=mu2,
                                         calc_actuals=FALSE, cv=FALSE,bt_mode=bt_mode)
          
          # Debug: Check fits structure
          if(is.null(fits)) {
            cat("DEBUG - fits is NULL\n")
          } else {
            cat("DEBUG - fits has", length(fits), "components\n")
            if("coefficients" %in% names(fits)) {
              cat("DEBUG - coefficients dimensions:", dim(fits$coefficients), "\n")
            }
          }
          
          eval = evaluateModels(fits, vg_sims=100)
          
          # Debug: Check eval structure
          if(is.null(eval)) {
            cat("DEBUG - eval is NULL\n")
          } else {
            cat("DEBUG - eval has", length(eval), "components:", names(eval), "\n")
          }
          
          return(eval)
        },
        args = list(
          dist = params$dist,
          a = params$a, 
          b = params$b,
          c = params$c,
          mu1 = params$mu1,
          mu2 = params$mu2,
          x1 = params$x1,
          x2 = params$x2,
          n = params$n,
          sim_seed = deterministic_seed,
          bt_mode = bt_mode
        ),
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
      # Debug: Check what we got back
      cat("DEBUG - Simulation succeeded, eval has", length(eval), "components\n")
      if("coefficients" %in% names(eval)) {
        cat("DEBUG - coefficients in eval has dimensions:", dim(eval$coefficients), "\n")
        cat("DEBUG - First few coefficient values:", head(eval$coefficients[1,]), "\n")
      }
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
process_parameter_set <- function(param_set_idx, params, outer_sims, parallel_strategy = "simulations", param_cores = NULL, mode = "broad_mean",bt_mode) {
  # Adjust core usage based on parallel strategy
  cores_to_use <- if (!is.null(param_cores) && parallel_strategy == "parameter_sets") {
    param_cores
  } else {
    n_cores
  }
  
  cat(sprintf("Parameter set %d using %d cores (strategy: %s)\n", param_set_idx, cores_to_use, parallel_strategy))
  
  # Check if output file already exists (cache check)
  if(mode=="broad_mean") {
    if(bt_mode==TRUE) {
    base_filename_pattern <- paste0(
      "Data/broad_sims_bt/CoefSimData_", paste(
        params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims,
        sep="_"
      ), "_*.RData"
    )
  } else if(bt_mode==FALSE) {
    base_filename_pattern <- paste0(
      "Data/broad_sims/CoefSimData_", paste(
        params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims,
        sep="_"
      ), "_*.RData"
    )
    } else {
      stop("bt_mode must be TRUE or FALSE")
    }
  } else {
    base_filename_pattern <- paste0(
      "Data/CoefSimData_", paste(
        params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims,
        sep="_"
      ), "_*.RData"
    )
  }
  
  
  existing_files <- Sys.glob(base_filename_pattern)
  
  if (length(existing_files) > 0) {
    cat(sprintf("Input set %d: Output file already exists, skipping: %s\n", param_set_idx, existing_files[1]))
    return(NULL)
  }
  
  # Define the save filename for new file
  if(mode == "broad_mean") {
    if(bt_mode==TRUE) {
      save_filename <- paste0(
        "Data/broad_sims_bt/CoefSimData_", paste(
          params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims, Sys.Date(),
          sep="_"
        ), ".RData"
      )
    } else if(bt_mode==FALSE) {
      save_filename <- paste0(
        "Data/broad_sims/CoefSimData_", paste(
          params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims, Sys.Date(),
          sep="_"
        ), ".RData"
      )
    } else {
      stop("bt_mode must be TRUE or FALSE")
    }
  } else if(mode == "specific_coef") {
    save_filename <- paste0(
      "Data/CoefSimData_", paste(
        params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims, Sys.Date(),
        sep="_"
      ), ".RData"
    )
  } else {
    stop("Unknown mode specified")
  }
  
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
          ifelse(is.na(params$a), 0, as.numeric(params$a) * 100),
          ifelse(is.na(params$b), 0, as.numeric(params$b) * 10),
          ifelse(is.na(params$c), 0, as.numeric(params$c) * 1),
          as.numeric(params$mu1) * 10000,
          as.numeric(params$mu2) * 1000,
          as.numeric(params$x1) * 100000,
          as.numeric(params$x2) * 1000000,
          as.numeric(params$n)
        ), na.rm = TRUE)
        
        chunk_results <- future_lapply(seq_along(start_idx:end_idx), function(chunk_i) {
          actual_sim_idx <- start_idx + chunk_i - 1
          run_single_simulation(params, actual_sim_idx, param_set_idx, bt_mode)
        }, future.seed = TRUE, future.scheduling = min(cores_to_use, length(start_idx:end_idx)))  # Enable deterministic seeding
        
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
        if(item %in% names(eval_outer[[i]])) {
          new_item <- data.frame(eval_outer[[i]][[item]])
          cat("DEBUG - Processing", item, "for sim", i, "- dimensions:", dim(new_item), "\n")
          for (j in seq_len(ncol(new_item))) {
            par_estimates[[z]][i, ] <- new_item[, j]
            z <- z + 1
          }
        } else {
          cat("DEBUG - Item", item, "not found in eval_outer[[", i, "]], available:", names(eval_outer[[i]]), "\n")
          # Skip the expected number of z positions
          z <- z + ifelse(item == "coefficients", 4, ifelse(item == "ses", 4, ifelse(item == "sigmas", 2, 1)))
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

if (parallel_strategy == "parameter_sets" && outer_sims == 1) {
  # Parameter-set-level parallelization for outer_sims = 1
  cat(sprintf("Running %d parameter sets in parallel (batch size: %d)\n", 
              length(input_params_list), param_sets_parallel))
  
  # Process parameter sets in batches
  num_batches <- ceiling(length(input_params_list) / param_sets_parallel)
  
  for (batch_idx in 1:num_batches) {
    start_idx <- (batch_idx - 1) * param_sets_parallel + 1
    end_idx <- min(batch_idx * param_sets_parallel, length(input_params_list))
    batch_indices <- start_idx:end_idx
    
    cat(sprintf("\n=== Processing Batch %d/%d (Parameter sets %d-%d) ===\n", 
                batch_idx, num_batches, start_idx, end_idx))
    
    # Setup separate future plan for this batch
    old_plan <- plan(multisession, workers = length(batch_indices))
    
    tryCatch({
      # Run parameter sets in parallel within this batch
      batch_results <- future_lapply(batch_indices, function(param_set_idx) {
        
        # Set deterministic seed for this parameter set worker
        # This ensures reproducible results regardless of parallel execution order
        worker_seed <- 12345 + param_set_idx * 1000  # Fixed base + parameter set index
        set.seed(worker_seed)
        
        # Create a wrapper function that handles the parameter set processing
        tryCatch({
          params <- input_params_list[[param_set_idx]]
          cat(sprintf("Starting parameter set %d in parallel worker (seed: %d)\n", param_set_idx, worker_seed))
          
          # Call the original process_parameter_set function with explicit parameters
          result <- process_parameter_set(
            param_set_idx = param_set_idx, 
            params = params, 
            outer_sims = outer_sims, 
            parallel_strategy = "parameter_sets", 
            param_cores = param_cores,
            mode = mode,
            bt_mode = bt_mode
          )
          
          list(
            param_set_idx = param_set_idx,
            result = result,
            success = TRUE,
            error = NULL
          )
          
        }, error = function(e) {
          list(
            param_set_idx = param_set_idx,
            result = NULL,
            success = FALSE,
            error = e$message
          )
        })
        
      }, future.packages = c("callr", "parallel", "future", "future.apply"),
         future.globals = list(
           input_params_list = input_params_list,
           process_parameter_set = process_parameter_set,
           run_single_simulation = run_single_simulation,
           outer_sims = outer_sims,
           mode = mode,
           param_cores = param_cores,
           n_cores = n_cores,
           parallel_strategy = parallel_strategy,
           bt_mode = bt_mode
         ),
         future.seed = TRUE)  # Enable deterministic seeding for reproducibility
      
      # Process batch results
      for (batch_result in batch_results) {
        param_set_idx <- batch_result$param_set_idx
        
        if (batch_result$success) {
          completed_sets <- c(completed_sets, param_set_idx)
          cat(sprintf("✓ Parameter set %d completed successfully\n", param_set_idx))
        } else {
          failed_sets <- c(failed_sets, param_set_idx)
          cat(sprintf("✗ Parameter set %d failed: %s\n", param_set_idx, batch_result$error))
          
          # Log the failure
          error_log <- paste0("Error in parameter set ", param_set_idx, " at ", Sys.time(), 
                             ":\n", batch_result$error, "\n\n")
          write(error_log, file = "parameter_set_errors.log", append = TRUE)
        }
      }
      
    }, error = function(e) {
      cat(sprintf("Batch %d failed: %s\n", batch_idx, e$message))
      # Mark all sets in this batch as failed
      for (idx in batch_indices) {
        if (!idx %in% completed_sets) {
          failed_sets <- c(failed_sets, idx)
        }
      }
    }, finally = {
      # Restore the original plan
      plan(old_plan)
    })
    
    # Progress summary after each batch
    cat(sprintf("Batch %d complete. Progress: %d completed, %d failed, %d remaining\n", 
                batch_idx, length(completed_sets), length(failed_sets), 
                length(input_params_list) - length(completed_sets) - length(failed_sets)))
    
    # Brief pause between batches for system stability
    if (batch_idx < num_batches) {
      Sys.sleep(2)
      gc()  # Force garbage collection between batches
    }
  }
  
} else {
  # Original sequential processing for outer_sims > 1
  cat("Using sequential parameter set processing with simulation-level parallelization\n")
  
  for (param_set_idx in seq_along(input_params_list)) {
    cat(sprintf("\n=== Processing Parameter Set %d/%d ===\n", param_set_idx, length(input_params_list)))
    
    params <- input_params_list[[param_set_idx]]
    
    tryCatch({
      result <- process_parameter_set(param_set_idx, params, outer_sims, "simulations", n_cores, mode, bt_mode)
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
