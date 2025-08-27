# Interactive Data File Selector
# This script creates a menu to select data files and automatically extracts parameters

library(tools)

# Function to extract parameters from filename
extract_params_from_filename <- function(filename) {
  # Remove extension and path
  basename_no_ext <- file_path_sans_ext(basename(filename))

  # Split by underscore
  parts <- strsplit(basename_no_ext, "_")[[1]]

  # Expected pattern: CoefSimData_DIST_a_b_c_mu1_mu2_x1_x2_n_runs_date
  if (length(parts) >= 11 && parts[1] == "CoefSimData") {
    return(list(
      dist = parts[2],
      a = as.numeric(parts[3]),
      b = as.numeric(parts[4]),
      c = as.numeric(parts[5]),
      mu1 = as.numeric(parts[6]),
      mu2 = as.numeric(parts[7]),
      x1 = as.numeric(parts[8]),
      x2 = as.numeric(parts[9]),
      n = as.numeric(parts[10]),
      runs = as.numeric(parts[11]),
      date = parts[12]
    ))
  } else {
    warning(paste("Could not parse filename:", filename))
    return(NULL)
  }
}

# Function to create interactive file selector
select_data_files <- function(data_dir = "Data") {
  # Find all .RData files in the Data directory
  rdata_files <- list.files(data_dir, pattern = "CoefSimData.*\\.RData$", full.names = TRUE)

  if (length(rdata_files) == 0) {
    stop("No CoefSimData .RData files found in ", data_dir)
  }

  cat("Available data files:\n")
  cat("====================\n\n")

  # Display files with extracted parameters
  for (i in seq_along(rdata_files)) {
    params <- extract_params_from_filename(rdata_files[i])
    if (!is.null(params)) {
      cat(sprintf("%d. %s\n", i, basename(rdata_files[i])))
      cat(sprintf("   Parameters: dist=%s, a=%g, b=%g, c=%g, mu1=%g, mu2=%g, x1=%g, x2=%g, n=%g\n\n",
                  params$dist, params$a, params$b, params$c,
                  params$mu1, params$mu2, params$x1, params$x2, params$n))
    }
  }

  # Interactive selection
  cat("Select files to process:\n")
  cat("- Enter numbers separated by commas (e.g., 1,3,5)\n")
  cat("- Enter 'all' to select all files\n")
  cat("- Enter 'q' to quit\n\n")

  selection <- readline("Your selection: ")

  if (tolower(selection) == "q") {
    return(NULL)
  }

  if (tolower(selection) == "all") {
    selected_indices <- seq_along(rdata_files)
  } else {
    # Parse comma-separated numbers
    selected_indices <- as.numeric(strsplit(gsub(" ", "", selection), ",")[[1]])
    selected_indices <- selected_indices[!is.na(selected_indices)]
    selected_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(rdata_files)]
  }

  if (length(selected_indices) == 0) {
    cat("No valid selection made.\n")
    return(NULL)
  }

  # Create input_datasets list
  input_datasets <- list()
  for (i in selected_indices) {
    params <- extract_params_from_filename(rdata_files[i])
    if (!is.null(params)) {
      input_datasets[[length(input_datasets) + 1]] <- list(
        rdata = rdata_files[i],
        meta = list(
          dist = params$dist,
          a = params$a,
          b = params$b,
          c = params$c,
          mu1 = params$mu1,
          mu2 = params$mu2,
          x1 = params$x1,
          x2 = params$x2,
          n = params$n
        )
      )
    }
  }

  cat(sprintf("\nSelected %d file(s) for processing.\n", length(input_datasets)))
  return(input_datasets)
}

# Function to auto-generate input_datasets from all files
auto_generate_input_datasets <- function(data_dir = "Data") {
  rdata_files <- list.files(data_dir, pattern = "CoefSimData.*\\.RData$", full.names = TRUE)

  input_datasets <- list()
  for (file in rdata_files) {
    params <- extract_params_from_filename(file)
    if (!is.null(params)) {
      input_datasets[[length(input_datasets) + 1]] <- list(
        rdata = file,
        meta = list(
          dist = params$dist,
          a = params$a,
          b = params$b,
          c = params$c,
          mu1 = params$mu1,
          mu2 = params$mu2,
          x1 = params$x1,
          x2 = params$x2,
          n = params$n
        )
      )
    }
  }

  cat(sprintf("Auto-generated input_datasets with %d file(s).\n", length(input_datasets)))
  return(input_datasets)
}

# Example usage:
cat("Data File Selector Loaded!\n")
cat("========================\n\n")
cat("Usage:\n")
cat("1. Interactive selection: input_datasets <- select_data_files()\n")
cat("2. Auto-select all files: input_datasets <- auto_generate_input_datasets()\n")
cat("3. Custom directory: input_datasets <- select_data_files('path/to/data')\n\n")
