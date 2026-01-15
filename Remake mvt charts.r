# Create charts for all coefficient files that don't have corresponding chart files
create_missing_charts <- function() {
  chart_dir <- "Charts"
  data_dir <- "Charts/ChartData"
  
  # Get all coefficient files
  coef_files <- list.files(data_dir, pattern = "_Coef\\.csv$", full.names = TRUE)
  
  for (coef_file in coef_files) {
    # Extract base filename without _Coef.csv
    base_name <- sub("_Coef\\.csv$", "", basename(coef_file))
    chart_file <- file.path(chart_dir, paste0(base_name, ".png"))
    
    # Check if chart already exists
    if (!file.exists(chart_file)) {
      cat("Creating chart for:", base_name, "\n")
      
      # Load the coefficient and SE files
      coef_data <- read.csv(coef_file, row.names = 1)
      se_file <- sub("_Coef\\.csv$", "_SE.csv", coef_file)
      se_data <- read.csv(se_file, row.names = 1)
      
      # Load the true coefficient and SE files directly as true_sim
      true_sim <- read.csv(file.path(data_dir, paste0(base_name, "_TrueSim.csv")), row.names = 1)
      true_coef_data <- true_sim[, c("Mean")]
      true_se_data <- true_sim[, c("SD")]
      
      # Prepare data for plotting
      model_names <- rownames(coef_data)
      colnames(coef_data) <- c("Intercept", "X1", "X2")
        
        # Create data frames for each coefficient
        plot_data_list <- lapply(1:ncol(coef_data), function(i) {
        data.frame(
          Model = factor(model_names, levels = model_names),
          Coefficient = coef_data[, i],
          SE = se_data[, i],
          TrueCoefficient = true_coef_data[i],  # Add true coefficient
          TrueSE = true_se_data[i],  # Add true SE
          Lower = coef_data[,i] - 1.96 * se_data[, i],
          Upper = coef_data[,i] + 1.96 * se_data[, i]
        )
        })
        
        # Create plots for each coefficient using the same style as in Testing Multivariate Fits (GA)
        plots <- lapply(1:ncol(coef_data), function(i) {
        y_min <- max(min(plot_data_list[[i]]$Lower, true_sim[i,1] - 1.96 * true_sim[i,2]),true_sim[i,1] - 10 * true_sim[i,2])
        y_max <- min(max(plot_data_list[[i]]$Upper, true_sim[i,1] + 1.96 * true_sim[i,2]),true_sim[i,1] + 10 * true_sim[i,2])
        y_range <- y_max - y_min
        
        ggplot(plot_data_list[[i]], aes(x = Model, y = Coefficient)) +
          geom_point(size = 3) +
          geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
          geom_hline(yintercept = true_sim[i,1], color = "red", linetype = "dashed") +
          geom_hline(yintercept = true_sim[i,1] - 1.96 * true_sim[i,2], color = "blue", linetype = "dashed") +
          geom_hline(yintercept = true_sim[i,1] + 1.96 * true_sim[i,2], color = "blue", linetype = "dashed") +
          labs(
          title = paste("Coefficient Estimates for", colnames(coef_data)[i]),
          x = NULL,
          y = "Coefficient Estimate"
          ) +
          theme_minimal() +
          theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none"
          ) +
          ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
        })

      # Extract distribution name from base_name
      #install.packages("readr")
      library(readr)

      len_base=length(strsplit(base_name, "_")[[1]])
      dist_name <- strsplit(base_name, "_")[[1]][2]
      skew_out <- parse_number(strsplit(base_name, "_")[[1]][len_base])

      rho_out <- parse_number(strsplit(base_name, "_")[[1]][len_base-2])
      sims_out <- parse_number(strsplit(base_name, "_")[[1]][len_base-1])
      if(dist_name=="Gamma") {dist_name_out="Gamma"} else if (dist_name=="Negative Binomial") {dist_name_out="Negative Binomial"} else if (dist_name=="Binomial") {dist_name_out="Binomial"} else {dist_name_out="Normal"}
      
      # Display all three plots together
      p <- grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3,
            top = paste("Multivariate ", dist_name_out, " (T=5) Coefficient Estimates with 95% CI | rho=",
                rho_out, " | skew=", skew_out, 
                "| sims=", sims_out, sep=""))

      ggsave(p, filename = chart_file, width = 9, height = 4)
      cat("Saved chart to:", chart_file, "\n")
    }
  }
  
  cat("Chart creation complete.\n")
}

create_missing_charts()


# Load and consolidate all LogLik files
consolidate_loglik_data <- function() {
  library(readr)
  library(tidyr)
  
  data_dir <- "Charts/ChartData"
  
  # Get all LogLik files
  loglik_files <- list.files(data_dir, pattern = "_LogLik\\.csv$", full.names = TRUE)
  
  cat("Found", length(loglik_files), "LogLik files\n")
  
  # Initialize lists to store data for each metric
  loglik_list <- list()
  df_list <- list()
  aic_list <- list()
  bic_list <- list()
  
  # Process each file
  for (file in loglik_files) {
    # Read the CSV file
    data <- read.csv(file, row.names = 1, check.names = FALSE)
    
    # Extract base filename
    base_name <- sub("_LogLik\\.csv$", "", basename(file))
    
    # Extract parameters from filename
    parts <- strsplit(base_name, "_")[[1]]
    len_parts <- length(parts)
    
    # parts[2] contains the distribution name (may have a space like "Negative Binomial")
    dist_name <- parts[2]
    if (dist_name == "Gamma") {
      dist_name_out <- "Gamma"
    } else if (dist_name == "Negative Binomial") {
      dist_name_out <- "Negative Binomial"
    } else if (dist_name == "Binomial") {
      dist_name_out <- "Binomial"
    } else if (dist_name == "Normal") {
      dist_name_out <- "Normal"
    } else {
      dist_name_out <- dist_name  # Fallback to original name
    }
    
    # Extract T, rho, sims, skew from the filename (always at the end)
    t_value <- parse_number(parts[len_parts - 3])  # T5
    rho <- parse_number(parts[len_parts - 2])       # rho0.75
    sims <- parse_number(parts[len_parts - 1])      # sims50
    skew <- parse_number(parts[len_parts])          # skew5.04
    
    # Create a row with metadata
    metadata <- data.frame(
      Distribution = dist_name_out,
      T = t_value,
      rho = rho,
      sims = sims,
      skew = skew
    )
    
    # Extract values for each model and add to lists
    for (model in rownames(data)) {
      row_data <- cbind(metadata, Model = model, 
                       LogLik = data[model, "LogLik"],
                       DF = data[model, "DF"],
                       AIC = data[model, "AIC"],
                       BIC = data[model, "BIC"])
      
      loglik_list[[length(loglik_list) + 1]] <- row_data[, c("Distribution", "T", "rho", "sims", "skew", "Model", "LogLik")]
      df_list[[length(df_list) + 1]] <- row_data[, c("Distribution", "T", "rho", "sims", "skew", "Model", "DF")]
      aic_list[[length(aic_list) + 1]] <- row_data[, c("Distribution", "T", "rho", "sims", "skew", "Model", "AIC")]
      bic_list[[length(bic_list) + 1]] <- row_data[, c("Distribution", "T", "rho", "sims", "skew", "Model", "BIC")]
    }
  }
  
  # Combine all rows into data frames
  loglik_long <- do.call(rbind, loglik_list)
  df_long <- do.call(rbind, df_list)
  aic_long <- do.call(rbind, aic_list)
  bic_long <- do.call(rbind, bic_list)
  
  # Reset row names
  rownames(loglik_long) <- NULL
  rownames(df_long) <- NULL
  rownames(aic_long) <- NULL
  rownames(bic_long) <- NULL
  
  # Pivot wider - make each model a separate column
  loglik_table <- pivot_wider(loglik_long, 
                              names_from = Model, 
                              values_from = LogLik)
  
  df_table <- pivot_wider(df_long, 
                          names_from = Model, 
                          values_from = DF)
  
  aic_table <- pivot_wider(aic_long, 
                           names_from = Model, 
                           values_from = AIC)
  
  bic_table <- pivot_wider(bic_long, 
                           names_from = Model, 
                           values_from = BIC)
  
  # Convert tibbles to data frames if needed
  loglik_table <- as.data.frame(loglik_table)
  df_table <- as.data.frame(df_table)
  aic_table <- as.data.frame(aic_table)
  bic_table <- as.data.frame(bic_table)
  
  # Round all model columns to 0 decimal places
  model_cols <- c("GLM", "GEE", "GAMLSS", "LME4", "GAMM", "VineCopula")
  
  for (col in model_cols) {
    if (col %in% colnames(loglik_table)) {
      loglik_table[[col]] <- round(loglik_table[[col]], 0)
    }
    if (col %in% colnames(df_table)) {
      df_table[[col]] <- round(df_table[[col]], 0)
    }
    if (col %in% colnames(aic_table)) {
      aic_table[[col]] <- round(aic_table[[col]], 0)
    }
    if (col %in% colnames(bic_table)) {
      bic_table[[col]] <- round(bic_table[[col]], 0)
    }
  }
  
  # Return a list with all tables
  result <- list(
    LogLik = loglik_table,
    DF = df_table,
    AIC = aic_table,
    BIC = bic_table
  )
  
  cat("\nConsolidated tables created:\n")
  cat("- LogLik table:", nrow(loglik_table), "rows\n")
  cat("- DF table:", nrow(df_table), "rows\n")
  cat("- AIC table:", nrow(aic_table), "rows\n")
  cat("- BIC table:", nrow(bic_table), "rows\n")
  
  return(result)
}

# Run the consolidation
consolidated_data <- consolidate_loglik_data()

# Access the tables
loglik_table <- consolidated_data$LogLik
df_table <- consolidated_data$DF
aic_table <- consolidated_data$AIC
bic_table <- consolidated_data$BIC

# Preview the LogLik table
cat("\nPreview of LogLik table:\n")
print(head(loglik_table, 20))

# Optional: Save the tables to CSV files
write.csv(loglik_table, "Charts/ChartData/mvt_Consolidated_LogLik.csv", row.names = FALSE)
write.csv(df_table, "Charts/ChartData/mvt_Consolidated_DF.csv", row.names = FALSE)
write.csv(aic_table, "Charts/ChartData/mvt_Consolidated_AIC.csv", row.names = FALSE)
write.csv(bic_table, "Charts/ChartData/mvt_Consolidated_BIC.csv", row.names = FALSE)

cat("\nConsolidated tables saved to Charts/ directory\n")


# Generate LaTeX table with best value bolded in each row
generate_model_comparison_table <- function(data_table, metric_name, select_max = TRUE, 
                                            output_file = NULL, caption = NULL, label = NULL) {
  # Round values for display
  display_table <- data_table
  model_cols <- c("GLM", "GAMLSS", "LME4", "VineCopula")  # Removed GAMM and GEE
  
  # Round the numeric columns
  for (col in model_cols) {
    if (col %in% colnames(display_table)) {
      display_table[[col]] <- round(display_table[[col]], 0)
    }
  }
  
  # Set defaults for output file, caption, and label
  if (is.null(output_file)) {
    output_file <- paste0("Charts/", metric_name, "_Table.tex")
  }
  if (is.null(caption)) {
    caption <- paste(metric_name, "Comparison Across Models and Distributions")
  }
  if (is.null(label)) {
    label <- paste0("tab:", tolower(metric_name))
  }
  
  # Start LaTeX table
  latex <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lrrrrcccc}",
    "\\toprule",
    "Distribution & T & $\\rho$ & Sims & Skew & GLM & GAMLSS & LME4 & VineCopula \\\\",
    "\\midrule"
  )
  
  # Add each row
  for (i in 1:nrow(display_table)) {
    # Get the model values for this row
    model_values <- as.numeric(display_table[i, model_cols])
    
    # Find the best value (maximum or minimum)
    if (select_max) {
      best_idx <- which.max(model_values)
    } else {
      best_idx <- which.min(model_values)
    }
    
    # Build the row
    row_parts <- c(
      display_table$Distribution[i],
      display_table$T[i],
      display_table$rho[i],
      display_table$sims[i],
      round(display_table$skew[i], 2)
    )
    
    # Add model values, bolding the best one
    for (j in 1:length(model_cols)) {
      if (j == best_idx) {
        row_parts <- c(row_parts, paste0("\\textbf{", model_values[j], "}"))
      } else {
        row_parts <- c(row_parts, as.character(model_values[j]))
      }
    }
    
    latex <- c(latex, paste(paste(row_parts, collapse = " & "), "\\\\"))
  }
  
  # Close the table
  latex <- c(
    latex,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  
  # Write to file
  writeLines(latex, output_file)
  cat("\nLaTeX table saved to", output_file, "\n")
  
  # Also print to console
  cat("\n--- LaTeX Table Code ---\n")
  cat(paste(latex, collapse = "\n"))
  cat("\n--- End LaTeX Table Code ---\n")
  
  return(latex)
}

# Generate the LaTeX tables
latex_loglik <- generate_model_comparison_table(loglik_table, "LogLik", select_max = TRUE,
                                                caption = "Log-Likelihood Comparison Across Models and Distributions")
latex_aic <- generate_model_comparison_table(aic_table, "AIC", select_max = FALSE,
                                             caption = "AIC Comparison Across Models and Distributions")
latex_bic <- generate_model_comparison_table(bic_table, "BIC", select_max = FALSE,
                                             caption = "BIC Comparison Across Models and Distributions")
