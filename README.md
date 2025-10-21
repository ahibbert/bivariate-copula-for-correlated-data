# bivariate-copula-for-correlated-data

# Introduction

This project provides all code used for the paper 'A comparison between copula-based, mixed model, and estimating equation methods for analysis of bivariate correlated data' written by Aydin Sareff-Hibbert and Gillian Z. Heller.

# Overview

This project contains the following components:
1. The code used to generate datasets, fit models and plot results for simlations presented in the paper, 'Simulations - 1 - Generate and Fit Models.R' and 'Simulations - 2 - Plot Models.R'. All fitted models data from (1) are saved into the Data/ folder while data on true standard errors from simulations are saved into the Cache/ folder and are generated as a part of plotting calculations in (2). All charts are saved in the Charts/ folder. These two scripts can be run sequentially to reproduce all simulation results in the paper.
2. The code used to fit a single model, used for a single applications case and for spot checking by the reader, '3 - Example - Fits and Plots.R'. This code can be run manually to reproduce the results for the applications case in the paper.
3. Miscellaneous code used for figures and various calculations denoted 'Ref - XX' 
4. Most main pieces of code rely on common functions provided by 'common_functions.R' and 'link_functions.R'

