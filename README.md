# Bivariate Copula for Correlated Data

This repository includes all the code required to replicate results included in the paper 'A comparison between copula-based, mixed model, and estimating equation methods for regression of bivariate correlated data' by Aydin Sareff-Hibbert and Gillian Z. Heller. 

It includes R code for simulating, fitting, and comparing models for correlated outcomes for bivariate and trivariate data.

## Overview

The following components of the paper map to the different sets of code:

| Component | Script |
|---|---|
| Bivariate simulations (sections 3.1-3.3) | 1, 2 (run sequentially) |
| Application (section 4) | 3 |
| Trivariate simulations (section 3.4) | 6 |

This repository contains research scripts for:

- generating correlated bivariate datasets under several outcome families,
- fitting multiple competing model classes,
- evaluating model performance across simulation settings,
- producing comparison plots, and
- running an applied example on longitudinal doctor-visit data.

The main workflow compares:

- `GLM`
- `GEE`
- `GAMLSS`
- `GAMLSS` mixture / random-effects variants
- `lme4` mixed models
- `mgcv::gamm`
- `GJRM`

## Repository layout

### Main scripts

- [1 - Simulations - Generate and Fit Models.R](1%20-%20Simulations%20-%20Generate%20and%20Fit%20Models.R)  
  Runs the main simulation-and-fitting pipeline over parameter grids stored in `input_params_list.RData`.

- [2 - Simulations - Plot Models.R](2%20-%20Simulations%20-%20Plot%20Models.R)  
  Aggregates saved simulation outputs and produces comparison plots.

- [2-1 Plotting - Evaluation_Metric_Plots.R](2-1%20Plotting%20-%20Evaluation_Metric_Plots.R)  
  Focused plotting for evaluation metrics.

- [2-2 Plotting - MU1_MU2_Coefficient_Plots.R](2-2%20Plotting%20-%20MU1_MU2_Coefficient_Plots.R)  
  Plots intercept / mean coefficient results.

- [2-3 Plotting - X1_Coefficient_Plots.R](2-3%20Plotting%20-%20X1_Coefficient_Plots.R)  
  Plots `x1` coefficient results.

- [2-4 Plotting - X2_Coefficient_Plots.R](2-4%20Plotting%20-%20X2_Coefficient_Plots.R)  
  Plots `x2` coefficient results.

- [3 - Example - Fits and Plots.R](3%20-%20Example%20-%20Fits%20and%20Plots.R)  
  Applied example using either real data or simulated data, followed by model fitting and plotting.

- [4 - Ref - Simulate Fitted Distributions.R](4%20-%20Ref%20-%20Simulate%20Fitted%20Distributions.R)  
  Reference script for simulation from fitted distributions.

- [6 - Trivariate Batch Runner.R](6%20-%20Trivariate%20Batch%20Runner.R)  
  Batch runner for the trivariate simulation workflow.

### Core support files

- [common_functions.R](common_functions.R)  
  Main function library for simulation, fitting, evaluation, and helper utilities.

- [link_functions.R](link_functions.R)  
  Link-function helpers used by the modelling code.

- [generate_input_params_list.R](generate_input_params_list.R)  
  Builds the parameter grid and saves it to `input_params_list.RData`.

- [start_httpgd.R](start_httpgd.R)  
  Starts an `httpgd` graphics device for browser-based plot viewing.

### Data and outputs

- `Data/` — input data and saved simulation results
- `Charts/` — generated figures
- `Cache/` — cached intermediate objects, including true SE calculations
- `Submission/` — output material prepared for submission
- `Archive/` — older or exploratory scripts retained for reference

## Statistical scope

The simulation code includes several outcome families identified by short labels in the scripts:

- `NO` — Gaussian / normal
- `PO` — count-data workflow based on Poisson-style settings
- `GA` — gamma
- `LO` — Bernoulli / logistic-style binary setting

Across these settings, the code studies estimation accuracy, uncertainty, fit statistics, and predictive scoring under correlated outcomes.

## Requirements

This project is written in R. The codebase uses a broad set of packages; the most important ones referenced in the main scripts include:

- `callr`
- `parallel`
- `future`
- `future.apply`
- `haven`
- `e1071`
- `gee`
- `geeM`
- `gamlss`
- `gamlss.mx`
- `lme4`
- `mgcv`
- `MASS`
- `GJRM`
- `VineCopula`
- `scoringRules`
- `ggplot2`
- `dplyr`
- `tidyr`
- `gridExtra`
- `cowplot`
- `ggpubr`
- `latex2exp`
- `glmtoolbox`
- `foreach`
- `doParallel`
- `httpgd`

A simple setup command is:

```r
install.packages(c(
  "callr", "future", "future.apply", "haven", "e1071", "gee", "geeM",
  "gamlss", "gamlss.mx", "lme4", "mgcv", "MASS", "VineCopula",
  "scoringRules", "ggplot2", "dplyr", "tidyr", "gridExtra", "cowplot",
  "ggpubr", "latex2exp", "glmtoolbox", "foreach", "doParallel", "httpgd"
))
```

Some packages may need to be installed from alternative repositories depending on the R version and platform. `GJRM` in particular should be checked separately if installation fails.

## Typical workflow

### 1. Generate the parameter grid

Run [generate_input_params_list.R](generate_input_params_list.R) to rebuild `input_params_list.RData`. Note this script is very long running and can take days to execute even with the implmented parallelisation.

### 2. Run simulation fits

Run [1 - Simulations - Generate and Fit Models.R](1%20-%20Simulations%20-%20Generate%20and%20Fit%20Models.R).

This script:

- loads the parameter grid,
- generates datasets,
- fits the candidate models in parallel,
- evaluates the fitted models, and
- saves the outputs for later plotting.

### 3. Plot simulation results

Run [2 - Simulations - Plot Models.R](2%20-%20Simulations%20-%20Plot%20Models.R), which calls other files as needed. Uses the outputs from [1 - Simulations - Generate and Fit Models.R](1%20-%20Simulations%20-%20Generate%20and%20Fit%20Models.R). Plots are placed in Charts/ directory.

### 4. Run the applied example

Run [3 - Example - Fits and Plots.R](3%20-%20Example%20-%20Fits%20and%20Plots.R).

This script supports:

- an applied longitudinal example based on RAND HRS doctor-visit counts, and
- a simulation example using the same fitting framework.

Unfortunately the RAND HRS data cannot be freely shared, but is available by request at the official source (details below).

### 5. Run trivariate experiments

Run [6 - Trivariate Batch Runner.R](6%20-%20Trivariate%20Batch%20Runner.R) for the trivariate batch simulations and plotting. Charts are placed in Charts/ directory. Note this script is very long running and can take days to execute even with the implmented parallelisation.

## Data notes

The applied example expects external data files that are not guaranteed to be distributable through this repository.

In particular, [3 - Example - Fits and Plots.R](3%20-%20Example%20-%20Fits%20and%20Plots.R) reads:
- `Data/randhrs1992_2020v2.sas7bdat`

The RAND HRS data should be obtained directly from the official source:
- https://www.rand.org/well-being/social-and-behavioral-policy/portfolios/aging-longevity/dataprod/hrs-data.html

## Notes on execution

- Several scripts assume the working directory is the repository root.
- Some fitting steps are isolated with `callr` because package interactions can affect model fitting in a shared R session. particularly if GJRM is ever loaded it often breaks functionality for GAMLSS model fits and functions so must either be loaded after or in an isolated session.
- Parallel settings are configured inside the simulation scripts and may need adjustment for the local machine, though sensible defaults are provided.
- Plot viewing can be improved by using `httpgd` via [start_httpgd.R](start_httpgd.R).

## Reproducibility

The simulation scripts use fixed seeds and deterministic seed construction in key sections to improve reproducibility across parameter sets and repeated runs.