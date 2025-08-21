# Load required library and source dependencies
library(callr)
source("common_functions.R")

# Define the list of input parameter sets as a list of lists
input_params_list <- list(
  list(dist="NO", a=1, b=1, c=0.25, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)#,
  #list(dist="NO", a=1, b=1, c=0.5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="NO", a=1, b=2, c=0.5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="LO",a=NA,b=NA,c=.25,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
  #list(dist="LO",a=NA,b=NA,c=.5,mu1=.25,mu2=.5,x1=1,x2=0.01,n=1000),
  #list(dist="NO", a=1, b=1, c=0.9, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="GA", a=.5, b=2, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="GA", a=1, b=1, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  #list(dist="GA", a=2, b=.5, c=NA, mu1=1, mu2=2, x1=1, x2=0.01, n=1000),
  
  #list(dist="PO", a=NA, b=1, c=5, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
)

outer_sims <- 10  # Number of simulations per input set

# Loop over each parameter set in the input_params_list
for (param_set_idx in seq_along(input_params_list)) {
  params <- input_params_list[[param_set_idx]]
  eval_outer <- list()
  max_retries <- 5
  for (i in 1:outer_sims) {
    cat(sprintf("Input set %d, simulation %d\n", param_set_idx, i))
    attempt <- 1
    success <- FALSE
    while (attempt <= max_retries && !success) {
      eval <- tryCatch({
        callr::r(
          func = function(dist, a, b, c, mu1, mu2, x1, x2, n) {
            source("common_functions.R")
            dataset = generateBivDist_withCov(n=n, a=a, b=b, c=c, mu1=mu1, mu2=mu2, dist=dist, x1=x1, x2=x2)
            # Set dataset as global as old S3 functions in gamlss and sometimes GJRM break without it
            assign("dataset", dataset, envir = .GlobalEnv)
            fits = fitBivModels_Bt_withCov(dataset=dataset, dist=dist, include="ALL",
                                           a=a, b=b, c=c, mu1=mu1, mu2=mu2,
                                           calc_actuals=FALSE, cv=FALSE)
            eval = evaluateModels(fits, vg_sims=100)
            return(eval)
          },
          args = params
        )
      }, error = function(e) {
        cat(sprintf("Error in callr::r attempt %d: %s\n", attempt, e$message))
        NULL
      })
      if (!is.null(eval)) {
        success <- TRUE
      } else {
        attempt <- attempt + 1
        if (attempt <= max_retries) {
          cat("Retrying...\n")
        }
      }
    }
    if (!success) {
      cat("All attempts failed for this simulation.\n")
      eval_outer[[i]] <- NA  # Or any other placeholder for failed runs
    } else {
      eval_outer[[i]] <- eval
    }
  }
  
  score_items <- list()
  score_item_names <- c("vs2_wt", "vs2", "es", "vs1", "vs2_wt_coronly", "logliks")
  for (i in 1:length(eval_outer)) {
    for (item in score_item_names) {
      new_item <- if (item == "logliks") {
        eval_outer[[i]][[item]][, 1]
      } else {
        eval_outer[[i]][[item]]
      }
      if (i == 1) {
        score_items[[item]] <- new_item
      } else {
        score_items[[item]] <- rbind(score_items[[item]], new_item)
      }
    }
  }
  
  par_estimates <- list()
  par_item_names <- c("coefficients", "ses", "sigmas", "correlations")
  for (i in 1:11) {
    par_estimates[[i]] <- matrix(nrow=length(eval_outer), ncol=16)
  }
  for (i in 1:length(eval_outer)) {
    z <- 1
    for (item in par_item_names) {
      new_item <- data.frame(eval_outer[[i]][[item]])
      for (j in 1:ncol(new_item)) {
        par_estimates[[z]][i, ] <- new_item[, j]
        z <- z + 1
      }
    }
  }
  names(par_estimates) <- c("mu1", "mu2", "x1", "x2", "mu1_se", "mu2_se", "x1_se", "x2_se", "sigma1", "sigma2", "corr")
  for (item in names(par_estimates)) {
    colnames(par_estimates[[item]]) <- c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
                                         "cop", "cop_n", "cop_j", "cop_g", "cop_f",
                                         "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t")
  }
  
  times=matrix(NA, nrow=length(eval_outer), ncol=length(colnames(par_estimates[[1]])))
  colnames(times) <- colnames(par_estimates[[1]])
  for (i in 1:length(eval_outer)) {
    time2=eval_outer[[i]][["timer"]][,2]
    time1=eval_outer[[i]][["timer"]][,1]
    timediff=difftime(time2, time1, units = "secs")
    times[i, ] <- c(timediff)
  }
  
  # Save the results for this parameter set
  save_filename <- paste0(
    "Data/CoefSimData_", paste(
      params$dist, params$a, params$b, params$c, params$mu1, params$mu2, params$x1, params$x2, params$n, outer_sims, Sys.Date(),
      sep="_"
    ), ".RData"
  )
  save(list = c("par_estimates", "score_items","times"), file = save_filename)
}
