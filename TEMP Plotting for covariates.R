library(ggplot2)
library(tidyr)
library(dplyr)
source("common_functions.R")

file_list=c("Data/results_combinedNO_1000_2025-07-08.RData"
            ,"Data/results_combinedGA_1000_2025-07-08.RData"
            ,"Data/results_combinedPO_1000_2025-07-08.RData"
            ,"Data/results_combinedLO_1000_2025-07-08.RData")

# Define the new order and new labels
model_order <- c("glm", "gee", "lme4", "re_nosig", "gamm", "re_np", "cop", "cop_n")
model_labels <- c("GLM", "GEE", "LME4", "GAMLSS", "GAMM", "GAMLSS NP", "GJRM (C)", "GJRM (N)")

# Plotting function (updated to use factor with new order and labels)
plot_results_gg_all_sims <- function(df, plot_title) {
  df$model <- factor(df$model, levels = model_order, labels = model_labels)
  ggplot(df, aes(x = model, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width=0.2) +
    geom_hline(aes(yintercept = true), color="red", linetype="dashed") +
    facet_grid(parameter ~ run, scales = "free_y") +
    labs(title = plot_title, x = "Model", y = "Estimate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Prepare and clean data for plotting all runs, now extract parameters to add to df
prepare_results_df <- function(results_combined, true_vals=rep(NA, 4)) {
  param_names <- c("t1","t2","x1","x2")
  dfs <- lapply(seq_along(results_combined), function(run) {
    results <- results_combined[[run]]
    est_df <- as.data.frame(results[[1]])
    se_df <- as.data.frame(results[[2]])*1.96 ### Making these into 95% confidence intervals
    colnames(est_df) <- param_names
    colnames(se_df) <- param_names
    est_df$model <- rownames(est_df)
    se_df$model <- rownames(se_df)
    
    # Extract parameters from results[[4]][[1]] according to order n,a,b,c,mu1,mu2,x1,x2
    # If result[[4]] exists and has at least one element
    if (!is.null(results[[4]]) && length(results[[4]]) >= 1) {
      params <- as.numeric(results[[4]][1,])
    } else {
      params <- rep(NA, 8)
    }
    names(params) <- c("n","a","b","c","mu1","mu2","x1","x2")
    
    # Add as columns to long df
    df_long <- est_df %>%
      pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "estimate") %>%
      left_join(se_df %>%
                  pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "se"),
                by = c("model", "parameter")
      ) %>%
      mutate(true = true_vals[match(parameter, param_names)],
             run = factor(run),
             n = params["n"], a = params["a"], b = params["b"], c = params["c"],
             mu1 = params["mu1"], mu2 = params["mu2"], x1 = params["x1"], x2 = params["x2"])
    df_long
  })
  do.call(rbind, dfs)
}

true_vals_list=list(c(1,2,1,1), c(log(1),log(2),1,1), c(log(1),log(2),1,1), c(logit_inv(1),logit_inv(2),1,1))
for (file in file_list) {
  # Extract the distribution from the filename  
  dist <- paste(regmatches(file, gregexpr("[A-Z]", file))[[1]][2],regmatches(file, gregexpr("[A-Z]", file))[[1]][3],sep="")
  load(file)  # loads results_combined
  df <- prepare_results_df(results_combined, true_vals = NA)
  # Exclude models that start with cop but are not cop or cop_n
  df <- df %>% filter(!(grepl("^cop", model) & !(model %in% c("cop", "cop_n"))))
  # Remove rows with NA estimates
  df$se[df$se < 0] <- NA
  
  df <- df %>%
    group_by(run, parameter) %>%
    mutate(glm_se = se[model == "glm"],
           se = ifelse(model %in% c("re_nosig", "re_np") & se > 10 * glm_se, NA, se)) %>%
    ungroup() %>%
    select(-glm_se)
  
  #Set the value of true for the four models for a given run based on a,b,c,mu1, mu2,x1,x2 columns
  if(dist=="NO") {
    df$true[df$parameter=="t1"] <- df$mu1[df$parameter=="t1"]
    df$true[df$parameter=="t2"] <- df$mu2[df$parameter=="t2"]
    df$true[df$parameter=="x1"] <- df$x1[df$parameter=="x1"]
    df$true[df$parameter=="x2"] <- df$x2[df$parameter=="x2"]
  } else if (dist=="GA") {
    df$true[df$parameter=="t1"] <- log(df$mu1[df$parameter=="t1"]*df$a[df$parameter=="t1"])
    df$true[df$parameter=="t2"] <- log(df$mu2[df$parameter=="t2"]*df$a[df$parameter=="t2"])
    df$true[df$parameter=="x1"] <- df$x1[df$parameter=="x1"]
    df$true[df$parameter=="x2"] <- df$x2[df$parameter=="x2"]
  } else if (dist=="PO") {
    df$true[df$parameter=="t1"] <- log(df$mu1[df$parameter=="t1"]*df$b[df$parameter=="t1"]*df$c[df$parameter=="t1"] )
    df$true[df$parameter=="t2"] <- log(df$mu2[df$parameter=="t2"]*df$b[df$parameter=="t1"]*df$c[df$parameter=="t1"])
    df$true[df$parameter=="x1"] <- df$x1[df$parameter=="x1"]
    df$true[df$parameter=="x2"] <- df$x2[df$parameter=="x2"]
  } else if (dist == "LO") {
    df$true[df$parameter == "t1"] <- logit(df$mu1[df$parameter=="t1"])
    df$true[df$parameter == "t2"] <- logit(df$mu2[df$parameter=="t2"])
    df$true[df$parameter == "x1"] <- df$x1[df$parameter=="x1"]
    df$true[df$parameter == "x2"] <- df$x2[df$parameter=="x2"]
  }
  
  print(plot_results_gg_all_sims(df, plot_title = basename(file)))
  readline(prompt = "Press [enter] to continue")
}