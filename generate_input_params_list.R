# Generate input_params_list from parameter specifications

#a=.5*1:5; b=.5*1:5;c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9); mu1=1; mu2=2; n=1000;dist="NO"
#a=NA; b=c(.2,.5,1,2,5);c=c(.2,.5,1,2,5); mu1=c(.5,1,2,5); mu2=c(.5,1,2,5); n=1000;dist="PO"
#a=.1+.1*1:20; b=.1+.1*1:20; c=NA; mu1=10; mu2=12; n=1000; dist="GA"
#a=NA; b=NA;c=c(.1,.25,.5,.75,.9); mu1=c(.1,.25,.5,.75,.9); mu2=c(.1,.25,.5,.75,.9); n=1000;dist="LO"

# Define parameter sets for each distribution
param_sets <- list(
  NO = list(
    a = 0.5 * 1:5,
    b = 0.5 * 1:5,
    c = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    mu1 = 1,
    mu2 = 2,
    n = 1000,
    dist = "NO"
  ),
  PO = list(
    a = NA,
    b = c(0.2, 0.5, 1, 2, 5),
    c = c(0.2, 0.5, 1, 2, 5),
    mu1 = c(0.5, 1, 2, 5),
    mu2 = c(0.5, 1, 2, 5),
    n = 1000,
    dist = "PO"
  ),
  GA = list(
    a = 0.1 + 0.1 * 1:20,
    b = 0.1 + 0.1 * 1:20,
    c = NA,
    mu1 = 10,
    mu2 = 12,
    n = 1000,
    dist = "GA"
  ),
  LO = list(
    a = NA,
    b = NA,
    c = c(0.1, 0.25, 0.5, 0.75, 0.9),
    mu1 = c(0.1, 0.25, 0.5, 0.75, 0.9),
    mu2 = c(0.1, 0.25, 0.5, 0.75, 0.9),
    n = 1000,
    dist = "LO"
  )
)

# Function to generate all combinations for a given distribution
generate_combinations <- function(params) {
  # Get all parameter names except 'dist'
  param_names <- names(params)[names(params) != "dist"]
  
  # Create a list of all parameters, handling NA values
  param_list <- list()
  for(name in param_names) {
    if(length(params[[name]]) == 1 && is.na(params[[name]])) {
      param_list[[name]] <- NA
    } else {
      param_list[[name]] <- params[[name]]
    }
  }
  
  # Generate all combinations using expand.grid
  combinations <- expand.grid(param_list, stringsAsFactors = FALSE)
  
  # Convert to list of lists and add x1, x2, and dist
  result <- list()
  for(i in seq_len(nrow(combinations))) {
    entry <- list(
      dist = params$dist,
      a = combinations$a[i],
      b = combinations$b[i],
      c = combinations$c[i],
      mu1 = combinations$mu1[i],
      mu2 = combinations$mu2[i],
      x1 = 1,
      x2 = .01,
      n = combinations$n[i]
    )
    result[[length(result) + 1]] <- entry
  }
  
  return(result)
}

# Generate combinations for all distributions
input_params_list <- list()

for(dist_name in names(param_sets)) {
  combinations <- generate_combinations(param_sets[[dist_name]])
  input_params_list <- c(input_params_list, combinations)
}

# Print summary of generated list
cat("Generated", length(input_params_list), "parameter combinations:\n")
for(dist_name in names(param_sets)) {
  count <- sum(sapply(input_params_list, function(x) x$dist == dist_name))
  cat("  ", dist_name, ":", count, "combinations\n")
}

# Display first few entries as example
cat("\nFirst few entries:\n")
for(i in seq_len(min(5, length(input_params_list)))) {
  cat("Entry", i, ":\n")
  print(input_params_list[[i]])
  cat("\n")
}

# Save the list
save(input_params_list, file = "input_params_list.RData")
cat("Saved input_params_list to input_params_list.RData\n")
