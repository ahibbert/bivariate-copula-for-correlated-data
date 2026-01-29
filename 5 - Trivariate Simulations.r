library(GJRM); source("common_functions.R"); source("link_functions.R");

simulate_trivariate <- function(
    n=1000, 
    mu_intercept=c(-1,0,1), 
    cutoff=0.5,
    mu_coefficients=c(1,0.01),
    theta_intercept=c(0.25,0.5,0.5),  #[2,3],[1,2],[1,3]|2
    copula_family=1,
    x1=NULL,
    x2=NULL,
    seed = NULL,
    plot_copula = FALSE
) {
    if (!is.null(seed)) set.seed(seed)

    if (length(mu_intercept) != 3) stop("mu_intercept must be length 3")
    if (length(mu_coefficients) != 2) stop("mu_coefficients must be length 2")
    if (length(theta_intercept) != 3) stop("theta_intercept must be length 3")

    # VineCopula provides RVineMatrix/RVineSim/contour methods
    if (!requireNamespace("VineCopula", quietly = TRUE)) {
        stop("Package 'VineCopula' is required for RVineMatrix/RVineSim. Install it via install.packages('VineCopula')")
    }

    Matrix <- matrix(
        c(3, 0, 0,
            1, 2, 0,
            2, 1, 1),
        nrow = 3, ncol = 3, byrow = TRUE
    )

    par <- matrix(
        c(0, 0, 0,
            theta_intercept[1], 0, 0,
            theta_intercept[3], theta_intercept[2], 0),
        nrow = 3, ncol = 3, byrow = TRUE
    )
    par2 <- 0 * par

    family <- matrix(
        c(0, 0, 0,
            copula_family, 0, 0,
            copula_family, copula_family, 0),
        nrow = 3, ncol = 3, byrow = TRUE
    )

    RVM <- VineCopula::RVineMatrix(Matrix = Matrix, family = family, par = par, par2 = par2)
    gaussian_copula_samples <- qnorm(VineCopula::RVineSim(n, RVM))

    if (is.null(x1)) x1 <- as.numeric(stats::runif(n, 0, 1) > 0.5)
    if (is.null(x2)) x2 <- sample(1:100, n, replace = TRUE)

    if (length(x1) != n) stop("x1 must be length n")
    if (length(x2) != n) stop("x2 must be length n")

    sim_true <- gaussian_copula_samples +
        matrix(rep(mu_intercept, n), ncol = 3, byrow = TRUE) +
        matrix(rep(x1, 3), byrow = FALSE, ncol = 3) * mu_coefficients[1] +
        matrix(rep(x2, 3), byrow = FALSE, ncol = 3) * mu_coefficients[2]

    y <- matrix(as.numeric(logit_inv(sim_true) > cutoff), ncol = 3)

    if (isTRUE(plot_copula)) {
        contour(RVM)
    }

    list(
        y = y,
        sim_true = sim_true,
        gaussian_copula_samples = gaussian_copula_samples,
        x1 = x1,
        x2 = x2,
        RVM = RVM,
        params = list(
            n = n,
            mu_intercept = mu_intercept,
            cutoff = cutoff,
            mu_coefficients = mu_coefficients,
            theta_intercept = theta_intercept,
            copula_family = copula_family
        )
    )
}

fit_trivariate_models <- function(
  sim,
  data_long,
  n = nrow(sim$y),
  d = ncol(sim$y),
  verbose = TRUE
) {
  if (!requireNamespace("gamlss", quietly = TRUE)) stop("Package 'gamlss' is required")
  if (!requireNamespace("glmtoolbox", quietly = TRUE)) stop("Package 'glmtoolbox' is required")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required")
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is required")
  if (!requireNamespace("GJRM", quietly = TRUE)) stop("Package 'GJRM' is required")

  library(gamlss)
  library(glmtoolbox)
  library(lme4)
  library(mgcv)
  library(GJRM)

  model_glm <- glm(y ~ -1 + t + x1 + x2,
    data = data_long,
    family = binomial(link = "logit")
  )

  model_gee <- tryCatch({
    glmgee(y ~ -1 + t + x1 + x2,
      data = data_long,
      family = binomial(link = "logit"),
      id = id,
      corstr = "Unstructured"
    )
  }, error = function(e) {
    warning(paste("glmgee failed:", conditionMessage(e)))
    NULL
  })

  model_re_nosig <- gamlss(
    formula = y ~ -1 + as.factor(t) + x1 + x2 + random(as.factor(id)),
    data = data_long,
    family = BI(),
    method = RS(100)
  )

  model_lme4 <- glmer(
    formula = y ~ -1 + t + x1 + x2 + (1 | id),
    data = data_long,
    family = binomial(link = "logit"),
    control = glmerControl(optCtrl = list(maxfun = 200000))
  )

  model_gamm <- tryCatch({
    gamm(
      formula = y ~ -1 + t + x1 + x2,
      random = list(id = ~1),
      data = data_long,
      family = binomial(link = "logit")
    )
  }, error = function(e) {
    warning(paste("gamm failed:", conditionMessage(e)))
    NULL
  })

  if (isTRUE(verbose)) print("Calculating effective degrees of freedom for LME4...")

  lme_EDF <- tryCatch({
    X <- getME(model_lme4, name = "X")[, 1:2]
    Z <- getME(model_lme4, name = "Z")
    U <- cbind(X, Z)
    W <- model_lme4@resp$sqrtrwt
    UWU <- (t(as.matrix(U)) %*% (diag(as.vector(W))) %*% as.matrix(U))
    D <- getME(model_lme4, name = "Lambda")

    if (sum(D) == 0) {
      NA_real_
    } else {
      D_inv <- solve(D)
      dinv_plus_00 <- c(0, 0, diag(D_inv))
      sum(diag(UWU %*% solve(UWU + diag(dinv_plus_00))))
    }
  }, error = function(e) {
    warning(paste("Failed to calculate LME4 EDF:", conditionMessage(e)))
    NA_real_
  })

  model_gjrm <- gjrm(
    formula = list(
      V1 ~ x1 + x2,
      V2 ~ x1 + x2,
      V3 ~ x1 + x2
    ),
    data = data.frame(
      V1 = sim$y[, 1],
      V2 = sim$y[, 2],
      V3 = sim$y[, 3],
      x1 = sim$x1,
      x2 = sim$x2
    ),
    margins = c("logit", "logit", "logit"),
    model = "T"
  )

  if (isTRUE(verbose)) print("Compiling results...")

  dfs <- c(
    (n * d) - df.residual(model_glm),
    if (is.null(model_gee)) NA else (n * d) - model_gee$df.residual,
    model_re_nosig$df.fit,
    lme_EDF,
    lme_EDF,
    model_gjrm$t.edf
  )

  logLiks <- c(
    logLik(model_glm),
    if (is.null(model_gee)) NA else model_gee$logLik,
    logLik(model_re_nosig),
    logLik(model_lme4),
    if (!is.null(model_gamm)) logLik(model_gamm$lme) else NA,
    logLik(model_gjrm)
  )

  results_table <- list()
  invisible(capture.output(results_table[[1]] <- if (!is.null(model_glm)) summary(model_glm)$coeff[, 1:2] else empty_coef_se()))
  invisible(capture.output(results_table[[2]] <- if (!is.null(model_gee)) summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients) - 2), 1:2] else empty_coef_se()))
  invisible(capture.output(results_table[[3]] <- if (!is.null(model_re_nosig)) cbind(summary(model_re_nosig)[1:5], summary(model_re_nosig)[6:10]) else empty_coef_se()))
  invisible(capture.output(results_table[[4]] <- if (!is.null(model_lme4)) summary(model_lme4)$coefficients[, c(1, 2)] else empty_coef_se()))
  invisible(capture.output(results_table[[5]] <- if (!is.null(model_gamm)) cbind(summary(model_gamm$lme)$coefficients[[1]], sqrt(diag(model_gamm$lme$varFix))) else empty_coef_se()))

  model_gjrm_vcov <- solve(model_gjrm$fit$hessian)

  model_gjrm_vcov_x1 <- model_gjrm_vcov[grep("x1", rownames(model_gjrm_vcov)), grep("x1", colnames(model_gjrm_vcov))]
  model_gjrm_vcov_x2 <- model_gjrm_vcov[grep("x2", rownames(model_gjrm_vcov)), grep("x2", colnames(model_gjrm_vcov))]
  model_gjrm_vcov_intercept <- model_gjrm_vcov[grep("(Intercept)", rownames(model_gjrm_vcov)), grep("(Intercept)", colnames(model_gjrm_vcov))]

  x1_var <- (sum(diag(model_gjrm_vcov_x1)) + 2 * (model_gjrm_vcov_x1[1, 2] + model_gjrm_vcov_x1[1, 3] + model_gjrm_vcov_x1[2, 3]))/(3^2)
  x2_var <- (sum(diag(model_gjrm_vcov_x2)) + 2 * (model_gjrm_vcov_x2[1, 2] + model_gjrm_vcov_x2[1, 3] + model_gjrm_vcov_x2[2, 3]))/(3^2)

  invisible(capture.output(results_table[[6]] <- if (!is.null(model_gjrm))
    cbind(
      c(model_gjrm$coefficients[c(1, 4, 7)], mean(model_gjrm$coefficients[c(2, 5, 8)]), mean(model_gjrm$coefficients[c(3, 6, 9)])),
      sqrt(c((diag(model_gjrm_vcov_intercept)), (x1_var), (x2_var)))
    )
  else empty_coef_se()))

  names(results_table) <- c("GLM", "GEE", "GAMLSS", "LME4", "GAMM", "GJRM")

  coefficients_table <- rbind(
    results_table[[1]][, 1],
    results_table[[2]][, 1],
    results_table[[3]][, 1],
    results_table[[4]][, 1],
    results_table[[5]][, 1],
    results_table[[6]][, 1]
  )
  ses_table <- rbind(
    results_table[[1]][, 2],
    results_table[[2]][, 2],
    results_table[[3]][, 2],
    results_table[[4]][, 2],
    results_table[[5]][, 2],
    results_table[[6]][, 2]
  )

  loglik_table <- cbind(
    logLiks,
    dfs,
    -2 * logLiks + 2 * dfs,
    -2 * logLiks + log(n * d) * dfs
  )

  conv_check <- c(
    model_glm$converged,
    if (is.null(model_gee)) NA else model_gee$converged,
    model_re_nosig$converged,
    if (!is.null(model_lme4)) !any(grepl("failed to converge", model_lme4@optinfo$conv$lme4$messages)) else NA,
    if (!is.null(model_gamm)) !any(grepl("converge", warnings(model_gamm))) else NA,
    model_gjrm$fit$converged
  )

  correlations <- list(
    0,
    if (is.null(model_gee)) NA else (model_gee$corr),
    getSmo(model_re_nosig)$sigb,
    summary(model_lme4)$varcor$id[1, 1],
    if (!is.null(model_gamm)) var(ranef(model_gamm$lme)[[1]]) else 0,
    c(model_gjrm$fit$theta12, model_gjrm$fit$theta23, model_gjrm$fit$theta13)
  )

  names(correlations)=rownames(loglik_table)=rownames(coefficients_table)=rownames(ses_table)=c("GLM","GEE","GAMLSS","LME4","GAMM","GJRM")

  list(
    coefficients = coefficients_table,
    ses = ses_table,
    logliks = loglik_table,
    convergences = conv_check,
    correlations = correlations
    #models = list(
    #  glm = model_glm,
    #  gee = model_gee,
    #  gamlss = model_re_nosig,
    #  lme4 = model_lme4,
    #  gamm = model_gamm,
    #  gjrm = model_gjrm
    #),
    #edf = list(lme4 = lme_EDF),
    #tables = results_table
  )
}

run_trivariate_outer_sims <- function(
  num_outer_sims,
  simulate_args = list(),
  fit_args = list(),
  seed_start = 1,
  verbose = TRUE
) {
  if (length(num_outer_sims) != 1 || is.na(num_outer_sims) || num_outer_sims < 1) {
    stop("num_outer_sims must be a single integer >= 1")
  }
  num_outer_sims <- as.integer(num_outer_sims)

  # We allocate after the first successful fit so we know dimensions.
  coeffs_arr <- NULL
  ses_arr <- NULL
  logliks_arr <- NULL
  convergences_mat <- NULL
  model_names <- NULL
  coef_names <- NULL
  loglik_colnames <- NULL

  seeds <- seed_start + (0:(num_outer_sims - 1))
  failures <- integer(0)

  for (s in seq_len(num_outer_sims)) {
    if (isTRUE(verbose)) message(sprintf("Outer sim %d / %d", s, num_outer_sims))

    sim_call_args <- simulate_args
    sim_call_args$seed <- seeds[s]
    if (is.null(sim_call_args$plot_copula)) sim_call_args$plot_copula <- FALSE

    sim <- do.call(simulate_trivariate, sim_call_args)

    data_long <- data.frame(
      y = as.vector(sim$y),
      x1 = rep(sim$x1, times = 3),
      x2 = rep(sim$x2, times = 3),
      t = factor(rep(c("1", "2", "3"), each = nrow(sim$y))),
      id = rep(1:nrow(sim$y), times = 3)
    )

    fit_call_args <- fit_args
    fit_call_args$sim <- sim
    fit_call_args$data_long <- data_long
    if (is.null(fit_call_args$verbose)) fit_call_args$verbose <- FALSE

    res <- tryCatch({
      do.call(fit_trivariate_models, fit_call_args)
    }, error = function(e) {
      warning(sprintf("fit_trivariate_models failed in outer sim %d: %s", s, conditionMessage(e)))
      NULL
    })

    if (is.null(res)) {
      failures <- c(failures, s)
      next
    }

    if (is.null(coeffs_arr)) {
      model_names <- rownames(res$coefficients)
      if (is.null(model_names)) model_names <- paste0("model", seq_len(nrow(res$coefficients)))

      coef_names <- colnames(res$coefficients)
      if (is.null(coef_names)) coef_names <- paste0("beta", seq_len(ncol(res$coefficients)))

      loglik_colnames <- colnames(res$logliks)
      if (is.null(loglik_colnames)) loglik_colnames <- paste0("stat", seq_len(ncol(res$logliks)))

      coeffs_arr <- array(NA_real_, dim = c(num_outer_sims, length(model_names), length(coef_names)))
      ses_arr <- array(NA_real_, dim = c(num_outer_sims, length(model_names), length(coef_names)))
      logliks_arr <- array(NA_real_, dim = c(num_outer_sims, length(model_names), ncol(res$logliks)))
      convergences_mat <- matrix(NA, nrow = num_outer_sims, ncol = length(model_names))

      dimnames(coeffs_arr) <- list(sim = seq_len(num_outer_sims), model = model_names, coef = coef_names)
      dimnames(ses_arr) <- list(sim = seq_len(num_outer_sims), model = model_names, coef = coef_names)
      dimnames(logliks_arr) <- list(sim = seq_len(num_outer_sims), model = model_names, stat = loglik_colnames)
      colnames(convergences_mat) <- model_names
      rownames(convergences_mat) <- seq_len(num_outer_sims)
    }

    # Store results (defensive: align by rownames where possible)
    coeffs_arr[s, , ] <- as.matrix(res$coefficients)
    ses_arr[s, , ] <- as.matrix(res$ses)
    logliks_arr[s, , ] <- as.matrix(res$logliks)
    convergences_mat[s, ] <- as.vector(res$convergences)
  }

  list(
    coefficients = coeffs_arr,
    ses = ses_arr,
    logliks = logliks_arr,
    convergences = convergences_mat,
    seeds = seeds,
    failures = failures
  )
}

calculate_trivariate_true_values <- function(
  true_sims,
  simulate_args = list(),
  seed_start = 1,
  verbose = TRUE
) {
  if (length(true_sims) != 1 || is.na(true_sims) || true_sims < 1) {
    stop("true_sims must be a single integer >= 1")
  }
  true_sims <- as.integer(true_sims)

  seeds <- seed_start + (0:(true_sims - 1))
  all_coefs <- matrix(NA_real_, ncol = 9, nrow = 0)

  for (i in seq_len(true_sims)) {
    if (isTRUE(verbose)) message(sprintf("True sim %d / %d", i, true_sims))

    sim_call_args <- simulate_args
    sim_call_args$seed <- seeds[i]
    if (is.null(sim_call_args$plot_copula)) sim_call_args$plot_copula <- FALSE

    sim <- do.call(simulate_trivariate, sim_call_args)


    ########NEED TO DO THIS WITH OPTIM INSTEAD

    coef <- c(
      stats::glm(sim$y[, 1] ~ sim$x1 + sim$x2, family = stats::binomial(link="logit"))$coefficients,
      stats::glm(sim$y[, 2] ~ sim$x1 + sim$x2, family = stats::binomial(link="logit"))$coefficients,
      stats::glm(sim$y[, 3] ~ sim$x1 + sim$x2, family = stats::binomial(link="logit"))$coefficients
    )

    all_coefs <- rbind(all_coefs, coef)
  }

  true_coef <- colMeans(all_coefs)

  true_covariance_matrix <- stats::cov(all_coefs)

  names(true_coef) <- c(
    "(Intercept)", "x1", "x2",
    "(Intercept)", "x1", "x2",
    "(Intercept)", "x1", "x2"
  )
  rownames(true_covariance_matrix) <- colnames(true_covariance_matrix) <- names(true_coef)

  model_true_vcov_x1 <- true_covariance_matrix[grep("x1", rownames(true_covariance_matrix)), grep("x1", colnames(true_covariance_matrix))]
  model_true_vcov_x2 <- true_covariance_matrix[grep("x2", rownames(true_covariance_matrix)), grep("x2", colnames(true_covariance_matrix))]
  model_true_vcov_intercept <- true_covariance_matrix[grep("(Intercept)", rownames(true_covariance_matrix)), grep("(Intercept)", colnames(true_covariance_matrix))]

  x1_var <- (sum(diag(model_true_vcov_x1)) + 2 * (model_true_vcov_x1[1, 2] + model_true_vcov_x1[1, 3] + model_true_vcov_x1[2, 3]))/(3^2)
  x2_var <- (sum(diag(model_true_vcov_x2)) + 2 * (model_true_vcov_x2[1, 2] + model_true_vcov_x2[1, 3] + model_true_vcov_x2[2, 3]))/(3^2)

  true_coef_se <- cbind(
    c(true_coef[c(1, 4, 7)], mean(true_coef[c(2, 5, 8)]), mean(true_coef[c(3, 6, 9)])),
    sqrt(c((diag(model_true_vcov_intercept)), (x1_var), (x2_var)))
  )
  colnames(true_coef_se) <- c("True Coef", "True SE")
  rownames(true_coef_se) <- c("t1", "t2", "t3", "x1", "x2")

  list(
    all_coefs = all_coefs,
    true_coef = true_coef,
    true_covariance_matrix = true_covariance_matrix,
    true_coef_se = true_coef_se,
    seeds = seeds
  )
}

##################### PARAMETERS ####################
true_sims=1000
n=1000;
mu_intercept=c(-1,0,1);
cutoff=0.5;
mu_coefficients=c(1,0.01);
theta_intercept=c(.025,.05,.05);  #[2,3],[1,2],[1,3]|2
copula_family=1
num_outer_sims=10

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

true_coef_vec <- true_vals$true_coef_se[, "True Coef"]
stopifnot(all(c("t1", "t2", "t3", "x1", "x2") %in% names(true_coef_vec)))

plot_one_coef <- function(coef_name) {
  df <- subset(coef_draws_long, coef == coef_name)
  ggplot(df, aes(x = model, y = estimate)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_hline(yintercept = as.numeric(true_coef_vec[coef_name]), linetype = "dashed") +
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

print(coef_boxplots_vs_true)


#################### PLOTS: SE BOXPLOTS VS TRUE ####################

se_draws_long <- as.data.frame.table(outer$ses, responseName = "se")
names(se_draws_long) <- c("sim", "model", "coef", "se")
se_draws_long$sim <- as.integer(as.character(se_draws_long$sim))
se_draws_long$se <- as.numeric(se_draws_long$se)

true_se_vec <- true_vals$true_coef_se[, "True SE"]
stopifnot(all(c("t1", "t2", "t3", "x1", "x2") %in% names(true_se_vec)))

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

print(se_boxplots_vs_true)


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

loglik_aic_bic_boxplots <- ggarrange(
  plotlist = ll_plots,
  ncol = 3,
  nrow = 1
)

loglik_aic_bic_boxplots <- ggpubr::annotate_figure(
  loglik_aic_bic_boxplots,
  top = ggpubr::text_grob(
    paste0("LogLik / AIC / BIC\n", simulation_params_title),
    size = 12,
    face = "bold"
  )
)

print(loglik_aic_bic_boxplots)

