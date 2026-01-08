# Master file for testing multivariate fits with different parameter values
# This file runs the testing scripts for different families (NO, GA, NB, LO)

# ============================================================================
# NORMAL COPULA (NO) - Testing Multivariate Fits
# ============================================================================

n=1000;d=5; sim_mean=1; sim_sigma=1; num_outer_sims=1; true_sims=100

cat("\n=== Running Testing Multivariate Fits (NO) with rho=0.5 ===\n")
rho <- 0.5
source("Testing Multivariate Fits (NO).R")

cat("\n=== Running Testing Multivariate Fits (NO) with rho=0.75 ===\n")
rho <- 0.75
source("Testing Multivariate Fits (NO).R")

# ============================================================================
# GAMMA COPULA (GA) - Testing Multivariate Fits
# ============================================================================
# Uncomment when ready to add GA testing
# cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.5 ===\n")
# rho <- 0.5
# source("Testing Multivariate Fits (GA).R")
#
# cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.75 ===\n")
# rho <- 0.75
# source("Testing Multivariate Fits (GA).R")

# ============================================================================
# NEGATIVE BINOMIAL COPULA (NB) - Testing Multivariate Fits
# ============================================================================
# Uncomment when ready to add NB testing
# cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.5 ===\n")
# rho <- 0.5
# source("Testing Multivariate Fits (NB).R")
#
# cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.75 ===\n")
# rho <- 0.75
# source("Testing Multivariate Fits (NB).R")

# ============================================================================
# LOGISTIC COPULA (LO) - Testing Multivariate Fits
# ============================================================================
# Uncomment when ready to add LO testing
# cat("\n=== Running Testing Multivariate Fits (LO) with rho=0.5 ===\n")
# rho <- 0.5
# source("Testing Multivariate Fits (LO).R")
#
# cat("\n=== Running Testing Multivariate Fits (LO) with rho=0.75 ===\n")
# rho <- 0.75
# source("Testing Multivariate Fits (LO).R")

cat("\n=== All tests completed ===\n")
