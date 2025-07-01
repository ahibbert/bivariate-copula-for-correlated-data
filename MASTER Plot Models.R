# Master script to run "Simulations - 2 - Plot Models.R" for different files_in and bt_mode values,
# showing live output from each run using callr::r_bg() with polling.

library(callr)

# Define the combinations you want to run
input_list <- list(
  list(files_in=c("Data/results_combined_B1_B2_NO_1000_2024-11-27.RData"
                  ,"Data/results_combined_B1_B2_PO_1000_2024-12-04.RData"
                  ,"Data/results_combined_B1_B2_GA_1000_2024-11-28.RData"
                  ,"Data/results_combined_B1_B2_LO_1000_2024-11-27.RData")
       , bt_mode = TRUE),
  list(files_in=c("Data/results_combined_B1_Bt_NO_1000_2024-12-05.RData"
                  ,"Data/results_combined_B1_Bt_PO_1000_2024-12-05.RData"
                  ,"Data/results_combined_B1_Bt_GA_1000_2024-11-28.RData"
                  ,"Data/results_combined_B1_Bt_LO_1000_2024-11-26.RData")
       , bt_mode = FALSE)
)

# Read the plotting script
script_lines <- readLines("Simulations - 2 - Plot Models.R")
# If files_in and bt_mode are defined in the first N lines, skip them:
start_line <- 3 # Adjust as needed!
script_to_run <- script_lines[start_line:length(script_lines)]

for (inputs in input_list) {
  # Start a background R process with the specified inputs
  p <- callr::r_bg(
    function(inputs, script) {
      list2env(inputs, envir = .GlobalEnv)
      eval(parse(text = script))
    },
    args = list(inputs = inputs, script = script_to_run),
    supervise = TRUE,
    stdout = "|",
    stderr = "|"
  )
  cat("---- Running with files_in =", inputs$files_in, ", bt_mode =", inputs$bt_mode, "----\n")
  # Poll for live output/errors
  while (p$is_alive()) {
    out <- p$read_output_lines()
    if (length(out)) cat(out, sep = "\n")
    err <- p$read_error_lines()
    if (length(err)) cat(err, sep = "\n")
    Sys.sleep(0.2)
  }
  # Print remaining output after process ends
  out <- p$read_all_output()
  if (length(out)) cat(out, sep = "\n")
  err <- p$read_all_error()
  if (length(err)) cat(err, sep = "\n")
  cat("---- Finished run ----\n\n")
}