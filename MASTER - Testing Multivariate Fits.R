# Master file for testing multivariate fits with different parameter values
# This file runs the testing scripts for different families (NO, GA, NB, LO)

set.seed(12345)
n=1000;d=5;
num_outer_sims=50; true_sims=10000;
# num_outer_sims=1; true_sims=5;


# ============================================================================
# NORMAL COPULA (NO) - Testing Multivariate Fits
# ============================================================================

sim_mean=1; sim_sigma=1; 

cat("\n=== Running Testing Multivariate Fits (NO) with rho=0.5 ===\n")
rho <- 0.5
source("Testing Multivariate Fits (NO).R")

cat("\n=== Running Testing Multivariate Fits (NO) with rho=0.75 ===\n")
rho <- 0.75
source("Testing Multivariate Fits (NO).R")

# ============================================================================
# GAMMA COPULA (GA) - Testing Multivariate Fits
# ============================================================================

sim_mean=1;dist_name="GA"

cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.5, low skew ===\n")
rho <- 0.5; sim_sigma=0.5
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.5, high skew ===\n")
rho <- 0.5; sim_sigma=2
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.75, low skew ===\n")
rho <- 0.75; sim_sigma=0.5
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (GA) with rho=0.75, high skew ===\n")
rho <- 0.75; sim_sigma=2
source("Testing Multivariate Fits (GA).R")

# ============================================================================
# NEGATIVE BINOMIAL COPULA (NB) - Testing Multivariate Fits
# ============================================================================

sim_mean=1;dist_name="NB"

cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.5, low skew ===\n")
rho <- 0.5; sim_sigma=.2; 
source("Testing Multivariate Fits (GA).R") # This is intentional, as GA script handles NB too

cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.5, high skew ===\n")
rho <- 0.5; sim_sigma=5; 
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.75, low skew  ===\n")
rho <- 0.75; sim_sigma=.2; 
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.75, high skew ===\n")
rho <- 0.75; sim_sigma=5; 
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (NB) with rho=0.75, extreme skew ===\n")
rho <- 0.75; sim_sigma=20; 
source("Testing Multivariate Fits (GA).R")

# ============================================================================
# LOGISTIC COPULA (LO) - Testing Multivariate Fits
# ============================================================================
sim_mean=1; dist_name="LO"; sim_sigma=NA;

cat("\n=== Running Testing Multivariate Fits (LO) with rho=0.5 ===\n")
rho <- 0.25
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (LO) with rho=0.5 ===\n")
rho <- 0.5
source("Testing Multivariate Fits (GA).R")

cat("\n=== Running Testing Multivariate Fits (LO) with rho=0.75 ===\n")
rho <- 0.75
source("Testing Multivariate Fits (GA).R")

cat("\n=== All tests completed ===\n")
