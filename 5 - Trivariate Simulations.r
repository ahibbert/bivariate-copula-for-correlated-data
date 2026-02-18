source("common_functions.R"); source("link_functions.R");


##################### PARAMETERS ####################
true_sims=100
num_outer_sims=100
n=1000;
mu_intercept=c(-2,-1,0); mu_coefficients=c(1,0.01); cutoff=0.5; theta_intercept=c(.75^2,.75,.75); 

#TESTING
#sim=simulate_trivariate(
#  n = n,
#  mu_intercept = mu_intercept,
#  cutoff = cutoff,
#  mu_coefficients = mu_coefficients,
#  theta_intercept = theta_intercept,
#  copula_family = copula_family
#)
#data_long <- data.frame(
#      y = as.vector(sim$y),
#      x1 = rep(sim$x1, times = 3),
#      x2 = rep(sim$x2, times = 3),
#      t = factor(rep(c("1", "2", "3"), each = nrow(sim$y))),
#      id = rep(1:nrow(sim$y), times = 3)
#    )
#res=fit_trivariate_models(sim=sim, data_long=data_long, verbose=TRUE)

variogram_call_args=list(
  fits = res,
  sim = sim
)
var= do.call(calc_variogram_score, variogram_call_args)

#################### CALCULATE TRUE VALUES VIA SIMULATION ####################

#Calculate true coefficients and covariance matrix via simulation
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

check=rbind(round(true_vals$true_coef,4)
  ,c(mu_intercept,mu_coefficients,logit(theta_intercept))
  ,round(true_vals$true_coef
  /c(mu_intercept,mu_coefficients,logit(theta_intercept)),2) )
rownames(check)=c("True Coef", "Simulated Params", "Ratio")
print(check)

#################### FIT MODELS ####################

#Example: run many outer simulations and extract coefficient/SE/loglik/convergence draws
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

true_coef_vec = true_vals$true_coef

names(true_coef_vec)=c("t1", "t2", "t3", "x1", "x2", "theta12", "theta23", "theta13")

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for plotting")
if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required for ggarrange")

library(ggplot2)
library(ggpubr)

stopifnot(exists("outer"), exists("true_vals"))

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

coef_draws_long <- as.data.frame.table(outer$coefficients, responseName = "estimate")
names(coef_draws_long) <- c("sim", "model", "coef", "estimate")
coef_draws_long$sim <- as.integer(as.character(coef_draws_long$sim))
coef_draws_long$estimate <- as.numeric(coef_draws_long$estimate)

plot_one_coef <- function(coef_name) {
  df <- subset(coef_draws_long, coef == coef_name)
  ggplot(df, aes(x = model, y = estimate)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_hline(yintercept = true_coef_vec[coef_name], linetype = "dashed") +
    labs(
      title = paste0("Coefficient: ", coef_name),
      x = "Model",
      y = "Estimated coefficient"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

coef_plots <- lapply(c("t1", "t2", "t3", "x1", "x2"), plot_one_coef)
names(coef_plots) <- c("t1", "t2", "t3", "x1", "x2")

coef_boxplots_vs_true <- ggarrange(
  plotlist = coef_plots,
  ncol = 5,
  nrow = 1
)

coef_boxplots_vs_true <- ggpubr::annotate_figure(
  coef_boxplots_vs_true,
  top = ggpubr::text_grob(
    paste0("Coefficients vs true\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

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
    labs(
      title = paste0("SE: ", coef_name),
      x = "Model",
      y = "Estimated SE"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

se_plots <- lapply(c("t1", "t2", "t3", "x1", "x2"), plot_one_se)
names(se_plots) <- c("t1", "t2", "t3", "x1", "x2")

se_boxplots_vs_true <- ggarrange(
  plotlist = se_plots,
  ncol = 5,
  nrow = 1
)

se_boxplots_vs_true <- ggpubr::annotate_figure(
  se_boxplots_vs_true,
  top = ggpubr::text_grob(
    paste0("SEs vs true\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

#print(se_boxplots_vs_true)
#################### PLOTS: LOGLIK / AIC / BIC BOXPLOTS ####################

stopifnot(exists("outer"))

loglik_draws_long <- as.data.frame.table(outer$logliks, responseName = "value")
names(loglik_draws_long) <- c("sim", "model", "stat", "value")
loglik_draws_long$sim <- as.integer(as.character(loglik_draws_long$sim))
loglik_draws_long$value <- as.numeric(loglik_draws_long$value)

# stat is a factor with levels in the original array order

count_logliks=sum(loglik_draws_long$stat=="logLiks")

loglik_draws_long$stat_index <- c(rep(1, times = count_logliks),
                                 rep(2, times = count_logliks),
                                 rep(3, times = count_logliks),
                                 rep(4, times = count_logliks)
)

# In outer$logliks: stat 1 = logLiks, stat 3 = AIC, stat 4 = BIC
loglik_draws_long <- subset(loglik_draws_long, stat_index %in% c(1, 3, 4))
loglik_draws_long <- subset(loglik_draws_long, model != "GAMM")
loglik_draws_long$stat_label <- factor(
  ifelse(loglik_draws_long$stat_index == 1, "LogLik",
    ifelse(loglik_draws_long$stat_index == 3, "AIC", "BIC")
  ),
  levels = c("LogLik", "AIC", "BIC")
)

plot_one_ll_stat <- function(stat_name) {
  df <- subset(loglik_draws_long, stat_label == stat_name)
  ggplot(df, aes(x = model, y = value)) +
    geom_boxplot(outlier.alpha = 0.4) +
    labs(
      title = paste0(stat_name, " by model"),
      x = "Model",
      y = stat_name
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
  value = as.vector(outer$vs2)
)
vs2_long$model <- factor(vs2_long$model, levels = colnames(outer$vs2))

vs2_wt_long <- data.frame(
  sim = rep(seq_len(nrow(outer$vs2_wt)), times = ncol(outer$vs2_wt)),
  model = rep(colnames(outer$vs2_wt), each = nrow(outer$vs2_wt)),
  value = as.vector(outer$vs2_wt)
)
vs2_wt_long$model <- factor(vs2_wt_long$model, levels = colnames(outer$vs2_wt))

ll_plots$VS2 <- ggplot(vs2_long, aes(x = model, y = value)) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(title = "VS2 by model", x = "Model", y = "VS2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ll_plots$VS2_wt <- ggplot(vs2_wt_long, aes(x = model, y = value)) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(title = "VS2 Weighted by model", x = "Model", y = "VS2 Weighted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

loglik_aic_bic_boxplots <- ggarrange(
  plotlist = ll_plots,
  ncol = 5,
  nrow = 1
)

loglik_aic_bic_boxplots <- ggpubr::annotate_figure(
  loglik_aic_bic_boxplots,
  top = ggpubr::text_grob(
    paste0("LogLik / AIC / BIC / VS2 / VS2 Weighted\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

#################### PLOTS: MEDIAN ESTIMATE ± 95% CI VS TRUE ####################

# Build summary: median coefficient and median SE per model × coef
coef_names_for_median <- c("t1", "t2", "t3", "x1", "x2")

model_levels <- levels(coef_draws_long$model)

median_summary <- do.call(rbind, lapply(coef_names_for_median, function(cn) {
  coef_df <- subset(coef_draws_long, coef == cn)
  se_df   <- subset(se_draws_long, coef == cn)
  models  <- levels(coef_df$model)
  do.call(rbind, lapply(models, function(m) {
    med_est <- median(coef_df$estimate[coef_df$model == m], na.rm = TRUE)
    med_se  <- median(se_df$se[se_df$model == m], na.rm = TRUE)
    data.frame(
      coef  = cn,
      model = m,
      med_est = med_est,
      med_se  = med_se,
      lo = med_est - qnorm(0.975) * med_se,
      hi = med_est + qnorm(0.975) * med_se,
      stringsAsFactors = FALSE
    )
  }))
}))

median_summary$model <- factor(median_summary$model, levels = model_levels)

plot_one_median_ci <- function(coef_name) {
  df <- subset(median_summary, coef == coef_name)
  true_val <- true_coef_vec[coef_name]
  true_se  <- as.numeric(true_se_vec[coef_name])
  true_lo  <- true_val - qnorm(0.975) * true_se
  true_hi  <- true_val + qnorm(0.975) * true_se

  # Compute symmetric y-axis limits centred on the true value
  all_vals <- c(df$lo, df$hi, true_lo, true_hi)
  max_dev  <- max(abs(all_vals - true_val))
  ylims    <- true_val + c(-1, 1) * max_dev * 1.05

  ggplot(df, aes(x = model, y = med_est)) +
    # True 95% CI band
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = true_lo, ymax = true_hi),
              fill = "lightblue", alpha = 0.3, inherit.aes = FALSE) +
    # True coefficient line
    geom_hline(yintercept = true_val, linetype = "dashed", colour = "blue") +
    # Model median estimate ± 95% CI
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
    coord_cartesian(ylim = ylims) +
    labs(
      title = coef_name,
      x = "Model",
      y = "Estimate"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

median_ci_plots <- lapply(coef_names_for_median, plot_one_median_ci)
names(median_ci_plots) <- coef_names_for_median

median_ci_vs_true <- ggarrange(
  plotlist = median_ci_plots,
  ncol = 5,
  nrow = 1
)

median_ci_vs_true <- ggpubr::annotate_figure(
  median_ci_vs_true,
  top = ggpubr::text_grob(
    paste0("Median estimate +/- 95% CI vs true (blue dashed / blue band)\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

print(median_ci_vs_true)
#################### PLOTS: MEAN ESTIMATE ± 95% CI VS TRUE ####################

mean_summary <- do.call(rbind, lapply(coef_names_for_median, function(cn) {
  coef_df <- subset(coef_draws_long, coef == cn)
  se_df   <- subset(se_draws_long, coef == cn)
  models  <- levels(coef_df$model)
  do.call(rbind, lapply(models, function(m) {
    mean_est <- mean(coef_df$estimate[coef_df$model == m], na.rm = TRUE)
    mean_se  <- mean(se_df$se[se_df$model == m], na.rm = TRUE)
    data.frame(
      coef  = cn,
      model = m,
      mean_est = mean_est,
      mean_se  = mean_se,
      lo = mean_est - qnorm(0.975) * mean_se,
      hi = mean_est + qnorm(0.975) * mean_se,
      stringsAsFactors = FALSE
    )
  }))
}))

mean_summary$model <- factor(mean_summary$model, levels = model_levels)

plot_one_mean_ci <- function(coef_name) {
  df <- subset(mean_summary, coef == coef_name)
  true_val <- true_coef_vec[coef_name]
  true_se  <- as.numeric(true_se_vec[coef_name])
  true_lo  <- true_val - qnorm(0.975) * true_se
  true_hi  <- true_val + qnorm(0.975) * true_se

  # Compute symmetric y-axis limits centred on the true value
  all_vals <- c(df$lo, df$hi, true_lo, true_hi)
  max_dev  <- max(abs(all_vals - true_val))
  ylims    <- true_val + c(-1, 1) * max_dev * 1.05

  ggplot(df, aes(x = model, y = mean_est)) +
    # True 95% CI band
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = true_lo, ymax = true_hi),
              fill = "lightblue", alpha = 0.3, inherit.aes = FALSE) +
    # True coefficient line
    geom_hline(yintercept = true_val, linetype = "dashed", colour = "blue") +
    # Model mean estimate ± 95% CI
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
    coord_cartesian(ylim = ylims) +
    labs(
      title = coef_name,
      x = "Model",
      y = "Estimate"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

mean_ci_plots <- lapply(coef_names_for_median, plot_one_mean_ci)
names(mean_ci_plots) <- coef_names_for_median

mean_ci_vs_true <- ggarrange(
  plotlist = mean_ci_plots,
  ncol = 5,
  nrow = 1
)

mean_ci_vs_true <- ggpubr::annotate_figure(
  mean_ci_vs_true,
  top = ggpubr::text_grob(
    paste0("Mean estimate +/- 95% CI vs true (blue dashed / blue band)\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

print(mean_ci_vs_true)
#################### SAVE MEAN CI PLOT ####################

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

ggsave(
  filename = file.path("Charts", paste0("trivariate_mean_ci_vs_true_", params_file_tag, ".png")),
  plot = mean_ci_vs_true,
  width = 18, height = 5, dpi = 300
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
  labels = c("Coefficients vs True", "SEs vs True", "Median Est +/- 95% CI vs True", "Mean Est +/- 95% CI vs True", "LogLik / AIC / BIC / VS2"),
  label.x = 0,
  label.y = 1,
  font.label = list(size = 12, face = "bold")
)

combined_plot <- ggpubr::annotate_figure(
  combined_plot,
  top = ggpubr::text_grob(
    simulation_params_title,
    size = 11,
    face = "bold"
  )
)

print(combined_plot)
#################### SAVE COMBINED PLOT ####################

ggsave(
  filename = file.path("Charts", paste0("trivariate_combined_", params_file_tag, ".png")),
  plot = combined_plot,
  width = 18, height = 22, dpi = 300
)
