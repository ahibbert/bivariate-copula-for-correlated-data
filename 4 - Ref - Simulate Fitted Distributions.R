run_sim_model_histograms <- function(dist = "NO", a = 1, b = 1, c = 0.5, mu1 = 1, mu2 = 2, x1 = 1, x2 = 0.01, n = 1000,log_in=FALSE) {
  source("common_functions.R")
  library(callr)
  
  dataset <- generateBivDist_withCov(n, a, b, c, mu1, mu2, dist, x1, x2)
  fits <- r(function(n, a, b, c, mu1, mu2, dist, x1, x2, dataset) {
    source("common_functions.R")
    assign("dataset", dataset, envir = .GlobalEnv)
    fitBivModels_Bt_withCov(data = dataset, dist = dist, include = "ALL", a = a, b = b, c = c, mu1 = mu1, mu2 = mu2)
  }, args = list(n = n, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2, dist = dist, x1 = x1, x2 = x2, dataset = dataset))
  
  out <- list()
  model_list <- c('glm', 'gee', 'lme4', 're_nosig', 'cop', 'cop_n')
  plot.new()
  par(mfrow = c(7, 2))

  for (model in model_list) {
    out[[model]] <- sim_model(model, dist, n, fits$coefficients, fits$sigmas, fits$correlations)
    if(log_in==TRUE){out[[model]]=log(out[[model]])}
  }

  print(  rbind(
          sapply(out, min)
          , sapply(out, max)
          , sapply(out, mean)
          , sapply(out, sd)
          ))
  
  min_val <- min(sapply(out, min))
  max_val <- max(sapply(out, max))
  time1=dataset$random_variable[dataset$time == 0]
  time2=dataset$random_variable[dataset$time == 1]
  if(log_in==TRUE){time1=log(time1);time2=log(time2)}

  
  hist(time1, breaks = 50,
    main = paste("Data |", dist, "| time==0"),
    xlab = 'y', xlim = c(min_val, max_val))
  hist(time2, breaks = 50,
    main = paste("Data |", dist, "| time==1"),
    xlab = 'y', xlim = c(min_val, max_val))
  for(model in model_list) {
    hist(out[[model]][dataset$time == 0], breaks = 50,
      main = paste(model, " | ", dist, " | time==0 | min=", round(min(out[[model]][dataset$time == 0]), 2), " max=", round(max(out[[model]][dataset$time == 0]), 2), sep = ""),
      xlab = 'y', xlim = c(min_val, max_val))
    hist(out[[model]][dataset$time == 1], breaks = 50,
      main = paste(model, " | ", dist, " | time==1 | min=", round(min(out[[model]][dataset$time == 1]), 2), " max=", round(max(out[[model]][dataset$time == 1]), 2), sep = ""),
      xlab = 'y', xlim = c(min_val, max_val))
  }
  
  dist_params <- paste0(dist, "_a", a, "_b", b, "_c", c, "_mu1", mu1, "_mu2", mu2)
  filename <- paste0("Charts/Additional/sim_model_histograms_", dist_params,log_in, ".png")
  dev.copy(png, filename = filename, width = 1200, height = 1800)
  dev.off()
}

# Example usage:
run_sim_model_histograms(dist="NO", a=1, b=1, c=0.9, mu1=1, mu2=2, x1=1, x2=0.01, n=1000)
run_sim_model_histograms(dist="GA", a=.5, b=2, c=NA, mu1=10, mu2=12, x1=1, x2=0.01, n=1000,log_in=FALSE)
run_sim_model_histograms(dist="GA", a=.5, b=2, c=NA, mu1=10, mu2=12, x1=1, x2=0.01, n=1000,log_in=TRUE)
run_sim_model_histograms(dist="GA", a=2, b=.5, c=NA, mu1=10, mu2=12, x1=1, x2=0.01, n=1000,log_in=FALSE)
run_sim_model_histograms(dist="GA", a=2, b=.5, c=NA, mu1=10, mu2=12, x1=1, x2=0.01, n=1000,log_in=TRUE)
run_sim_model_histograms(dist="PO", a=2, b=2, c=2, mu1=1, mu2=2, x1=1, x2=0.01, n=1000,log_in=FALSE)
run_sim_model_histograms(dist="PO", a=2, b=2, c=2, mu1=1, mu2=2, x1=1, x2=0.01, n=1000,log_in=TRUE)
