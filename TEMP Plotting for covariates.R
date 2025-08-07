library(ggplot2)
library(tidyr)
library(dplyr)
source("common_functions.R")

file_list = c(
  "Data/results_combinedNO_1000_2025-07-08.RData",
  "Data/results_combinedGA_1000_2025-07-08.RData",
  "Data/results_combinedPO_1000_2025-07-08.RData",
  "Data/results_combinedLO_1000_2025-07-08.RData"
)

model_order <- c("glm", "gee", "lme4", "re_nosig", "gamm", "re_np", "cop", "cop_n")
model_labels <- c("GLM", "GEE", "LME4", "GAMLSS", "GAMM", "GAMLSS NP", "GJRM (C)", "GJRM (N)")

# Chart for estimates (with error bars)
plot_estimates_gg <- function(df, plot_title) {
  df$model <- factor(df$model, levels = model_order, labels = model_labels)
  ggplot(df, aes(x = model, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.2) +
    geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
    facet_grid(parameter ~ run_label, scales = "free_y") +
    labs(title = plot_title, x = "Model", y = "Estimate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Chart for standard errors (against true value)
plot_se_gg <- function(df, plot_title) {
  df$model <- factor(df$model, levels = model_order, labels = model_labels)
  ggplot(df, aes(x = model, y = se)) +
    geom_point() +
    geom_hline(aes(yintercept = 0), color = "black", linetype = "dotted") +
    # ADD: true SE horizontal line
    geom_hline(aes(yintercept = true_se), color = "red", linetype = "dashed") +
    facet_grid(parameter ~ run_label, scales = "free_y") +
    labs(title = plot_title, x = "Model", y = "Standard Error (95% CI)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

prepare_results_df <- function(results_combined, true_vals = rep(NA, 4), true_ses = rep(NA, 4)) {
  param_names <- c("t1", "t2", "x1", "x2")
  dfs <- lapply(seq_along(results_combined), function(run) {
    results <- results_combined[[run]]
    est_df <- as.data.frame(results[[1]])
    se_df <- as.data.frame(results[[2]]) * 1.96
    colnames(est_df) <- param_names
    colnames(se_df) <- param_names
    est_df$model <- rownames(est_df)
    se_df$model <- rownames(se_df)
    if (!is.null(results[[4]]) && length(results[[4]]) >= 1) {
      params <- as.numeric(results[[4]][1, ])
    } else {
      params <- rep(NA, 8)
    }
    names(params) <- c("n", "a", "b", "c", "mu1", "mu2", "x1", "x2")
    run_label <- paste0("a=", signif(params["a"], 3),
                        ", b=", signif(params["b"], 3),
                        ", c=", signif(params["c"], 3),
                        ", mu1=", signif(params["mu1"], 3),
                        ", mu2=", signif(params["mu2"], 3))
    # true_se from argument, will be joined below
    df_long <- est_df %>%
      pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "estimate") %>%
      left_join(se_df %>%
                  pivot_longer(cols = all_of(param_names), names_to = "parameter", values_to = "se"),
                by = c("model", "parameter")
      ) %>%
      mutate(true = true_vals[match(parameter, param_names)],
             true_se = true_ses[match(parameter, param_names)],
             run = factor(run),
             n = params["n"], a = params["a"], b = params["b"], c = params["c"],
             mu1 = params["mu1"], mu2 = params["mu2"], x1 = params["x1"], x2 = params["x2"],
             run_label = run_label)
    df_long
  })
  do.call(rbind, dfs)
}

# NEW: Helper to get true SEs for one run
get_true_ses <- function(n, a, b, c, mu1, mu2, dist, x1, x2) {
  # sims=100 or 200 recommended for stable SE estimate
  sim_out <- simCovariateMLEs(sims = 200, n = n, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2,
                              dist = dist, x1 = x1, x2 = x2, trace = FALSE)
  # The returned ses vector is: c("mu1","mu2","x1","x2","s1","s2") -- only first 4 are for parameters of interest
  ses <- sim_out$ses[1:4]
  names(ses) <- c("t1", "t2", "x1", "x2")
  return(ses)
}

for (file in file_list) {
  dist <- paste(regmatches(file, gregexpr("[A-Z]", file))[[1]][2], regmatches(file, gregexpr("[A-Z]", file))[[1]][3], sep = "")
  load(file)
  # Get parameters from first run (usually same for all runs)
  if (!is.null(results_combined[[1]][[4]]) && length(results_combined[[1]][[4]]) >= 1) {
    params <- as.numeric(results_combined[[1]][[4]][1, ])
    names(params) <- c("n", "a", "b", "c", "mu1", "mu2", "x1", "x2")
    # Get true SEs for this scenario
    true_ses <- get_true_ses(params["n"], params["a"], params["b"], params["c"],
                             params["mu1"], params["mu2"], dist,
                             params["x1"], params["x2"])
  } else {
    true_ses <- rep(NA, 4)
  }
  df <- prepare_results_df(results_combined, true_vals = NA, true_ses = true_ses)
  df <- df %>% filter(!(grepl("^cop", model) & !(model %in% c("cop", "cop_n"))))
  df$se[df$se < 0] <- NA
  df <- df %>%
    group_by(run, parameter) %>%
    mutate(glm_se = se[model == "glm"],
           se = ifelse(model %in% c("re_nosig", "re_np") & se > 10 * glm_se, NA, se)) %>%
    ungroup() %>%
    dplyr::select(-glm_se)
  if (dist == "NO") {
    df$true[df$parameter == "t1"] <- df$mu1[df$parameter == "t1"]
    df$true[df$parameter == "t2"] <- df$mu2[df$parameter == "t2"]
    df$true[df$parameter == "x1"] <- df$x1[df$parameter == "x1"]
    df$true[df$parameter == "x2"] <- df$x2[df$parameter == "x2"]
  } else if (dist == "GA") {
    df$true[df$parameter == "t1"] <- log(df$mu1[df$parameter == "t1"] * df$a[df$parameter == "t1"])
    df$true[df$parameter == "t2"] <- log(df$mu2[df$parameter == "t2"] * df$a[df$parameter == "t2"])
    df$true[df$parameter == "x1"] <- df$x1[df$parameter == "x1"]
    df$true[df$parameter == "x2"] <- df$x2[df$parameter == "x2"]
  } else if (dist == "PO") {
    df$true[df$parameter == "t1"] <- log(df$mu1[df$parameter == "t1"] * df$b[df$parameter == "t1"] * df$c[df$parameter == "t1"])
    df$true[df$parameter == "t2"] <- log(df$mu2[df$parameter == "t2"] * df$b[df$parameter == "t1"] * df$c[df$parameter == "t1"])
    df$true[df$parameter == "x1"] <- df$x1[df$parameter == "x1"]
    df$true[df$parameter == "x2"] <- df$x2[df$parameter == "x2"]
  } else if (dist == "LO") {
    df$true[df$parameter == "t1"] <- logit(df$mu1[df$parameter == "t1"])
    df$true[df$parameter == "t2"] <- logit(df$mu2[df$parameter == "t2"])
    df$true[df$parameter == "x1"] <- df$x1[df$parameter == "x1"]
    df$true[df$parameter == "x2"] <- df$x2[df$parameter == "x2"]
  }
  
  # Chart 1: Estimates vs true value
  p_est <- plot_estimates_gg(df, plot_title = paste("Estimates for", dist, "distribution"))
  outfile_est <- file.path("Charts", paste0("Charts_", dist, "_EST_", tools::file_path_sans_ext(basename(file)), ".png"))
  ggsave(outfile_est, plot = p_est, width = 3 * length(results_combined), height = 9, dpi = 600)
  
  # Chart 2: SE vs true value (yintercept = 0, plus horizontal line for true SE)
  p_se <- plot_se_gg(df, plot_title = paste("Standard Errors for", dist, "distribution"))
  outfile_se <- file.path("Charts", paste0("Charts_", dist, "_SE_", tools::file_path_sans_ext(basename(file)), ".png"))
  ggsave(outfile_se, plot = p_se, width = 3 * length(results_combined), height = 9, dpi = 600)
}