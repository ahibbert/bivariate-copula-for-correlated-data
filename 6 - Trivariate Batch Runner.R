##############################################################################
# Batch runner for 5 - Trivariate Simulations
# Loops over multiple parameter sets, running the full simulation + plotting
# workflow for each one.
##############################################################################

source("common_functions.R"); source("link_functions.R")

library(ggplot2)
library(ggpubr)

##################### SHARED PARAMETERS ####################
true_sims   <- 10000
num_outer_sims <- 1000
n           <- 1000
copula_family <- 1

##################### PARAMETER GRID ####################
# Each row is one scenario: list(mu_intercept, mu_coefficients, cutoff, theta_intercept)

param_grid <- list(
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.9^2, .9, .9)),
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.75^2, .75, .75)),
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.667^2, .667, .667)),
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.5^2, .5, .5)),
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.25^2, .25, .25)),
  list(mu_intercept = c(-2,-1,0), mu_coefficients = c(1,0.01), cutoff = 0.5,  theta_intercept = c(.1^2, .1, .1))
)

##################### MAIN LOOP ####################

for (idx in seq_along(param_grid)) {

  p <- param_grid[[idx]]
  mu_intercept     <- p$mu_intercept
  mu_coefficients  <- p$mu_coefficients
  cutoff           <- p$cutoff
  theta_intercept  <- p$theta_intercept

  message(sprintf(
    "\n========== Scenario %d / %d ==========\n  mu_intercept=[%s]  mu_coefficients=[%s]  cutoff=%s  theta_intercept=[%s]\n",
    idx, length(param_grid),
    paste(mu_intercept, collapse = ","),
    paste(mu_coefficients, collapse = ","),
    cutoff,
    paste(theta_intercept, collapse = ",")
  ))

  params_file_tag <- paste0(
    "n", n,
    "_mu", paste(mu_intercept, collapse = "-"),
    "_cut", cutoff,
    "_coef", paste(mu_coefficients, collapse = "-"),
    "_theta", paste(theta_intercept, collapse = "-"),
    "_cop", copula_family,
    "_tsims", true_sims,
    "_nsims", num_outer_sims
  )

  simulation_params_title <- paste0(
    "Trivariate simulation (logit margin; R-vine copula): ",
    "n=", n,
    ", mu_intercept=[", paste(mu_intercept, collapse = ","), "]",
    ", cutoff=", cutoff,
    ", mu_coefficients=[", paste(mu_coefficients, collapse = ","), "]",
    ", theta_intercept=[", paste(theta_intercept, collapse = ","), "]",
    ", copula_family=", copula_family,
    " | true_sims=", true_sims,
    ", num_outer_sims=", num_outer_sims
  )

  #################### CALCULATE TRUE VALUES ####################

  true_vals <- calculate_trivariate_true_values(
    true_sims = true_sims,
    simulate_args = list(
      n = n,
      mu_intercept = mu_intercept,
      cutoff = cutoff,
      mu_coefficients = mu_coefficients,
      theta_intercept = theta_intercept,
      copula_family = copula_family,
      x1 = NULL,
      x2 = NULL
    ),
    seed_start = 1,
    verbose = TRUE
  )

  check <- rbind(
    round(true_vals$true_coef, 4),
    c(mu_intercept, mu_coefficients, logit(theta_intercept)),
    round(true_vals$true_coef / c(mu_intercept, mu_coefficients, logit(theta_intercept)), 2)
  )
  rownames(check) <- c("True Coef", "Simulated Params", "Ratio")
  print(check)

  #################### FIT MODELS ####################

  outer <- run_trivariate_outer_sims(
    num_outer_sims = num_outer_sims,
    simulate_args = list(
      n = n,
      mu_intercept = mu_intercept,
      cutoff = cutoff,
      mu_coefficients = mu_coefficients,
      theta_intercept = theta_intercept,
      copula_family = copula_family
    ),
    seed_start = 1000,
    verbose = TRUE
  )

  #################### PLOTS: COEFFICIENT BOXPLOTS VS TRUE ####################

  true_coef_vec <- true_vals$true_coef
  names(true_coef_vec) <- c("t1", "t2", "t3", "x1", "x2", "theta12", "theta23", "theta13")

  # Mathematical notation labels for coefficients
  coef_labels <- list(
    t1 = expression(beta[1]),
    t2 = expression(beta[2]),
    t3 = expression(beta[3]),
    x1 = expression(beta[x[1]]),
    x2 = expression(beta[x[2]])
  )

  coef_draws_long <- as.data.frame.table(outer$coefficients, responseName = "estimate")
  names(coef_draws_long) <- c("sim", "model", "coef", "estimate")
  coef_draws_long$sim <- as.integer(as.character(coef_draws_long$sim))
  coef_draws_long$estimate <- as.numeric(coef_draws_long$estimate)

  plot_one_coef <- function(coef_name) {
    df <- subset(coef_draws_long, coef == coef_name)
    ggplot(df, aes(x = model, y = estimate)) +
      geom_boxplot(outlier.alpha = 0.4) +
      geom_hline(yintercept = true_coef_vec[coef_name], linetype = "dashed") +
      labs(title = coef_labels[[coef_name]], x = NULL, y = "Estimated coefficient") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))
  }

  coef_plots <- lapply(c("t1", "t2", "t3", "x1", "x2"), plot_one_coef)
  names(coef_plots) <- c("t1", "t2", "t3", "x1", "x2")

  #################### PLOTS: SE BOXPLOTS VS TRUE ####################

  true_se_vec <- true_vals$true_coef_se
  names(true_se_vec) <- c("t1", "t2", "t3", "x1", "x2", "theta12", "theta23", "theta13")

  se_draws_long <- as.data.frame.table(outer$ses, responseName = "se")
  names(se_draws_long) <- c("sim", "model", "coef", "se")
  se_draws_long$sim <- as.integer(as.character(se_draws_long$sim))
  se_draws_long$se <- as.numeric(se_draws_long$se)

  plot_one_se <- function(coef_name) {
    df <- subset(se_draws_long, coef == coef_name)
    ggplot(df, aes(x = model, y = se)) +
      geom_boxplot(outlier.alpha = 0.4) +
      geom_hline(yintercept = as.numeric(true_se_vec[coef_name]), linetype = "dashed") +
      labs(title = coef_labels[[coef_name]], x = NULL, y = "Estimated SE") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))
  }

  se_plots <- lapply(c("t1", "t2", "t3", "x1", "x2"), plot_one_se)
  names(se_plots) <- c("t1", "t2", "t3", "x1", "x2")

  #################### PLOTS: LOGLIK / AIC / BIC / VS2 / VS2_WT ####################

  # Define consistent color palette across all model-selection plots
  all_models <- levels(coef_draws_long$model)
  model_colours <- setNames(scales::hue_pal()(length(all_models)), all_models)

  loglik_draws_long <- as.data.frame.table(outer$logliks, responseName = "value")
  names(loglik_draws_long) <- c("sim", "model", "stat", "value")
  loglik_draws_long$sim <- as.integer(as.character(loglik_draws_long$sim))
  loglik_draws_long$value <- as.numeric(loglik_draws_long$value)

  count_logliks <- sum(loglik_draws_long$stat == "logLiks")
  loglik_draws_long$stat_index <- c(
    rep(1, times = count_logliks),
    rep(2, times = count_logliks),
    rep(3, times = count_logliks),
    rep(4, times = count_logliks)
  )

  loglik_draws_long <- subset(loglik_draws_long, stat_index %in% c(1, 3, 4))
  loglik_draws_long <- subset(loglik_draws_long, model != "GAMM" & !is.na(model))
  loglik_draws_long$model <- droplevels(loglik_draws_long$model)
  loglik_draws_long$stat_label <- factor(
    ifelse(loglik_draws_long$stat_index == 1, "LogLik",
      ifelse(loglik_draws_long$stat_index == 3, "AIC", "BIC")),
    levels = c("LogLik", "AIC", "BIC")
  )

  plot_one_ll_stat <- function(stat_name) {
    df <- subset(loglik_draws_long, stat_label == stat_name)
    summ <- do.call(rbind, lapply(split(df, df$model), function(d) {
      data.frame(model = d$model[1],
                 median = median(d$value, na.rm = TRUE),
                 lo = quantile(d$value, 0.025, na.rm = TRUE),
                 hi = quantile(d$value, 0.975, na.rm = TRUE))
    }))
    ggplot(summ, aes(x = model, y = median, colour = model)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
      scale_colour_manual(values = model_colours) +
      labs(title = paste0(stat_name), x = NULL, y = stat_name) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5))
  }

  ll_plots <- list(
    LogLik = plot_one_ll_stat("LogLik"),
    AIC = plot_one_ll_stat("AIC"),
    BIC = plot_one_ll_stat("BIC")
  )

  # VS2 and VS2_wt boxplots
  vs2_long <- data.frame(
    sim = rep(seq_len(nrow(outer$vs2)), times = ncol(outer$vs2)),
    model = rep(colnames(outer$vs2), each = nrow(outer$vs2)),
    value = (as.vector(outer$vs2))
  )
  vs2_long$model <- factor(vs2_long$model, levels = colnames(outer$vs2))

  vs2_wt_long <- data.frame(
    sim = rep(seq_len(nrow(outer$vs2_wt)), times = ncol(outer$vs2_wt)),
    model = rep(colnames(outer$vs2_wt), each = nrow(outer$vs2_wt)),
    value = (as.vector(outer$vs2_wt))
  )
  vs2_wt_long$model <- factor(vs2_wt_long$model, levels = colnames(outer$vs2_wt))

  make_pointrange <- function(df_long, title_text, y_label) {
    summ <- do.call(rbind, lapply(split(df_long, df_long$model), function(d) {
      data.frame(model = d$model[1],
                 median = median(d$value, na.rm = TRUE),
                 lo = quantile(d$value, 0.025, na.rm = TRUE),
                 hi = quantile(d$value, 0.975, na.rm = TRUE))
    }))
    ggplot(summ, aes(x = model, y = median, colour = model)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
      scale_colour_manual(values = model_colours) +
      labs(title = title_text, x = NULL, y = y_label) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5))
  }

  ll_plots$VS2 <- make_pointrange(vs2_long, "Variogram Score", "VS2")
  ll_plots$VS2_wt <- make_pointrange(vs2_wt_long, "Variogram Score (Wt)", "VS2 Weighted")

  var_plot=ggarrange(plotlist = ll_plots, ncol = 5, nrow = 1)

  ggsave(
    filename = file.path("Charts", paste0("trivariate_model_selection_", params_file_tag, ".png")),
    plot = var_plot,
    width = 11, height = 4, dpi = 300
  )

  #################### PLOTS: MEDIAN ESTIMATE ± 95% CI VS TRUE ####################

  coef_names_for_median <- c("t1", "t2", "t3", "x1", "x2")
  model_levels <- levels(coef_draws_long$model)

  median_summary <- do.call(rbind, lapply(coef_names_for_median, function(cn) {
    coef_df <- subset(coef_draws_long, coef == cn)
    se_df   <- subset(se_draws_long, coef == cn)
    models  <- levels(coef_df$model)
    do.call(rbind, lapply(models, function(m) {
      med_est <- median(coef_df$estimate[coef_df$model == m], na.rm = TRUE)
      med_se  <- median(se_df$se[se_df$model == m], na.rm = TRUE)
      data.frame(coef = cn, model = m, med_est = med_est, med_se = med_se,
                 lo = med_est - qnorm(0.975) * med_se,
                 hi = med_est + qnorm(0.975) * med_se,
                 stringsAsFactors = FALSE)
    }))
  }))
  median_summary$model <- factor(median_summary$model, levels = model_levels)

  plot_one_median_ci <- function(coef_name) {
    df <- subset(median_summary, coef == coef_name)
    true_val <- true_coef_vec[coef_name]
    true_se  <- as.numeric(true_se_vec[coef_name])
    true_lo  <- true_val - qnorm(0.975) * true_se
    true_hi  <- true_val + qnorm(0.975) * true_se
    all_vals <- c(df$lo, df$hi, true_lo, true_hi)
    max_dev  <- max(abs(all_vals - true_val))
    ylims    <- true_val + c(-1, 1) * max_dev * 1.05

    ggplot(df, aes(x = model, y = med_est)) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = true_lo, ymax = true_hi),
                fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
      geom_hline(yintercept = true_val, linetype = "dashed", colour = "black") +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
      coord_cartesian(ylim = ylims) +
      labs(title = coef_labels[[coef_name]], x = NULL, y = "Estimate") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))
  }

  median_ci_plots <- lapply(coef_names_for_median, plot_one_median_ci)
  names(median_ci_plots) <- coef_names_for_median

  #################### PLOTS: MEAN ESTIMATE ± 95% CI VS TRUE ####################

  mean_summary <- do.call(rbind, lapply(coef_names_for_median, function(cn) {
    coef_df <- subset(coef_draws_long, coef == cn)
    se_df   <- subset(se_draws_long, coef == cn)
    models  <- levels(coef_df$model)
    do.call(rbind, lapply(models, function(m) {
      mean_est <- mean(coef_df$estimate[coef_df$model == m], na.rm = TRUE)
      mean_se  <- mean(se_df$se[se_df$model == m], na.rm = TRUE)
      data.frame(coef = cn, model = m, mean_est = mean_est, mean_se = mean_se,
                 lo = mean_est - qnorm(0.975) * mean_se,
                 hi = mean_est + qnorm(0.975) * mean_se,
                 stringsAsFactors = FALSE)
    }))
  }))
  mean_summary$model <- factor(mean_summary$model, levels = model_levels)

  plot_one_mean_ci <- function(coef_name) {
    df <- subset(mean_summary, coef == coef_name)
    true_val <- true_coef_vec[coef_name]
    true_se  <- as.numeric(true_se_vec[coef_name])
    true_lo  <- true_val - qnorm(0.975) * true_se
    true_hi  <- true_val + qnorm(0.975) * true_se
    all_vals <- c(df$lo, df$hi, true_lo, true_hi)
    max_dev  <- max(abs(all_vals - true_val))
    ylims    <- true_val + c(-1, 1) * max_dev * 1.05

    ggplot(df, aes(x = model, y = mean_est, colour = model)) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = true_lo, ymax = true_hi),
                fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
      geom_hline(yintercept = true_val, linetype = "dashed", colour = "black") +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
      scale_colour_manual(values = model_colours) +
      coord_cartesian(ylim = ylims) +
      labs(title = coef_labels[[coef_name]], x = NULL, y = "Estimate + 95% CI") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5))
  }

  mean_ci_plots <- lapply(coef_names_for_median, plot_one_mean_ci)
  names(mean_ci_plots) <- coef_names_for_median

  #################### SAVE MEAN CI PLOT ####################

  mean_ci_vs_true <- ggarrange(plotlist = mean_ci_plots, ncol = 5, nrow = 1)
  #mean_ci_vs_true <- ggpubr::annotate_figure(
  #  mean_ci_vs_true,
  #  top = ggpubr::text_grob(
  #    paste0("Mean estimate +/- 95% CI vs true (blue dashed / blue band)\n", simulation_params_title),
  #    size = 12, face = "bold"
  #  )
  #)


  ggsave(
    filename = file.path("Charts", paste0("trivariate_mean_ci_vs_true_", params_file_tag, ".png")),
    plot = mean_ci_vs_true,
    width = 11, height = 4, dpi = 300
  )

  #################### COMBINED PLOT ####################

  combined_plot <- ggarrange(
    ggarrange(plotlist = coef_plots, ncol = 5, nrow = 1),
    ggarrange(plotlist = se_plots, ncol = 5, nrow = 1),
    ggarrange(plotlist = median_ci_plots, ncol = 5, nrow = 1),
    ggarrange(plotlist = mean_ci_plots, ncol = 5, nrow = 1),
    ggarrange(plotlist = ll_plots, ncol = 5, nrow = 1),
    nrow = 5,
    heights = c(1, 1, 1, 1, 1),
    labels = c("Coefficients vs True", "SEs vs True", "Median Est +/- 95% CI vs True",
               "Mean Est +/- 95% CI vs True", "LogLik / AIC / BIC / VS2"),
    label.x = 0, label.y = 1,
    font.label = list(size = 12, face = "bold")
  )

  combined_plot <- ggpubr::annotate_figure(
    combined_plot,
    top = ggpubr::text_grob(simulation_params_title, size = 11, face = "bold")
  )

  ggsave(
    filename = file.path("Charts", paste0("trivariate_combined_", params_file_tag, ".png")),
    plot = combined_plot,
    width = 18, height = 22, dpi = 300
  )

  message(sprintf("Scenario %d / %d complete. Charts saved.", idx, length(param_grid)))
}

message("\n========== All scenarios complete ==========")
