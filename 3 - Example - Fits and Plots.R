############# Overview of file #####################

# CHOOSE EITHER FROM 1) Applications data or 2) Simulation data before moving to 3) to fit all models. 
# These create an object 'dataset' that can then be fit with functions from common_functions.

# Section 3) Plots the given dataset and fits all specified models to the given data for comparison in a simple table.
# Section 4) Provides a function to fit all models multiple times at different sample sizes then time it
# Section 5) Function for fitting a GJRM model with manual gamlss ZIS margins and a VineCopula fit for the copula

#############Required functions #################
source("common_functions.R"); source("link_functions.R")
set.seed(1000);options(scipen=999);

####### 1) Applications Data #############
#### RAND HLS - Data is available on application from https://www.rand.org/well-being/social-and-behavioral-policy/portfolios/aging-longevity/dataprod/hrs-data.html

library(haven)

#Load dataset with only required columns from base randHRS file available at link above
rand <- read_sas("Data/randhrs1992_2020v2.sas7bdat",col_select=c("HHIDPN","R14IWENDY","R15IWENDY","RAGENDER","RABYEAR"
                                                                 ,"R14DOCTIM" #Doctor visits last 2 years
                                                                 ,"R15DOCTIM" #Doctor visits last 2 years
))

#load("Data/rand_mvt.rds")

rand_doc_visits <-as.data.frame(rand[!(is.na(rand$R14DOCTIM))&!(is.na(rand$R15DOCTIM))&rand$R14DOCTIM>=0&rand$R15DOCTIM>=0, ])
rand_doc_visits_SAMPLED<-rand_doc_visits[runif(nrow(rand_doc_visits))<0.2,]
rand_doc_visits=rand_doc_visits_SAMPLED
rm(rand)

#save(rand_doc_visits_SAMPLED,file="Data/rand_doc_visits_SAMPLED")
#load("Data/rand_doc_visits_SAMPLED") ; dist="PO" #Used once the data is sampled once to ensure the same sample across models. 
dist="PO"
gamma_c_mu1<-as.vector(rand_doc_visits[,"R14DOCTIM"])
gamma_c_mu2<-as.vector(rand_doc_visits[,"R15DOCTIM"])
gender=as.vector(rand_doc_visits[,"RAGENDER"])-1
age=as.vector(2020-rand_doc_visits[,"RABYEAR"])

# Set up as longitudinal structured data - MUST DO THIS STEP TO SETUP APPLICATIONS DATA CORRECTLY FOR MODEL FITTING 
patient<-as.factor(seq(1:length(gamma_c_mu1)))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0,gender,age),cbind(patient,gamma_c_mu2,1,gender,age)))
colnames(dataset)<-c("patient","random_variable","time","sex","age")
dataset<-dataset[order(dataset$patient),]

a=NA; b=NA; c=NA; mu1=NA; mu2=NA; n=NA #Dummy values to pass to function


#### 2) Simulation Data #############

#n=1000;a=NA;b=NA;c=0.75;mu1=.25;mu2=.5;dist="LO";x1=1;x2=0.01
#dataset=generateBivDist_withCov(n=1000,a=NA,b=NA,c=0.75,mu1=.25,mu2=.5,dist="LO",x1=1,x2=0.01)
#results<-fitBivModels_Bt_withCov(data=dataset,dist=dist,include="ALL",a,b,c,mu1,mu2)
#eval_out=evaluateModels(results,rownames(results$correlations),vg_sims = 100)


####### 3) Plot and fit data ###########################################################

# Plot margins and correlation structure
plotDist(dataset,dist)

# Print basic dataset info: mean, skewness, correlation
library(e1071)
mean(dataset$random_variable[dataset$time==0])
mean(dataset$random_variable[dataset$time==1])
log(mean(dataset$random_variable[dataset$time==1])/mean(dataset$random_variable[dataset$time==0]))
skewness(dataset$random_variable[dataset$time==0])
skewness(dataset$random_variable[dataset$time==1])
cor(dataset$random_variable[dataset$time==0],dataset$random_variable[dataset$time==1],method="kendall")
cor(dataset$random_variable[dataset$time==0],dataset$random_variable[dataset$time==1],method="pearson")

# Fit models to the dataset and cleanup results 
 

# Running in isolated callr session to avoid issues with GJRM package breaking gamlss functionality
library(callr)

results <- r(function(data_input, dist, a, b, c, mu1, mu2) {
  # Source required functions in the isolated session
  source("common_functions.R")
  source("link_functions.R")
  
  # Make dataset available globally in the isolated session
  assign("dataset", data_input, envir = .GlobalEnv)
  
  # Call the function
  fitBivModels_Bt_withCov(data=data_input, dist=dist, include="ALL", a, b, c, mu1, mu2)
}, args = list(data_input = dataset, dist = dist, a = a, b = b, c = c, mu1 = mu1, mu2 = mu2))

coef=cbind(results$coefficients,results$ses,results$logliks)
colnames(coef)=c("mu1","mu2","x1","x2","mu1_se","mu2_se","x1_se","x2_se","loglik","edf")
coef

# Save workspace after fitting models
#save.image(file="workspace_after_fitting.RData")
#load("workspace_after_fitting.RData")

####### 5) Manual fitting of ZIS #############

set.seed(1000)

library(gamlss); source("common_functions.R"); source("link_functions.R")
#GAMLSS MODEL
#model_re_nosig_re
model_re_nosig <- gamlss(formula=random_variable~ -1+as.factor(time==1)+(sex)+age+random(as.factor(patient))
        , tau.formula = ~1
        , nu.formula = ~1
        , sigma.formula = ~1
        , data=dataset, family=ZISICHEL(), method=CG(100))

links<-ZISICHEL()

sigma_try=0 # Sigma ft is not significant compared to zero in model so we fix it at zero.
# Fit margins separately
fit1<-gamlss(gamma_c_mu1~(gender)+age
  , tau.formula = ~1
  , nu.formula = ~1
  , sigma.formula = ~1
  , family=ZISICHEL(), method=CG(1000)
  , tau.start = links$tau.linkinv( model_re_nosig$tau.coefficients )
  , nu.start = links$nu.linkinv(model_re_nosig$nu.coefficients )
  #, sigma.start = links$sigma.linkinv(model_re_nosig$sigma.coefficients)
  , sigma.start = exp(sigma_try)
  , sigma.fix=TRUE
  )

fit2<-gamlss(gamma_c_mu2~(gender)+age
  , tau.formula = ~1
  , nu.formula = ~1
  , sigma.formula = ~1
  , family=ZISICHEL()
  , method=CG(1000)
  , tau.start = links$tau.linkinv( model_re_nosig$tau.coefficients )
  , nu.start = links$nu.linkinv(model_re_nosig$nu.coefficients )
  #, sigma.start = links$sigma.linkinv(model_re_nosig$sigma.coefficients)
  , sigma.start = exp(sigma_try)
  , sigma.fix=TRUE
  ) 

fit2$sigma.coefficients=sigma_try; fit1$sigma.coefficients=sigma_try

mu1_all <- links$mu.linkinv(fit1$mu.coefficients[1] + fit1$mu.coefficients[2]*gender + fit1$mu.coefficients[3]*age)
mu2_all <- links$mu.linkinv(fit2$mu.coefficients[1] + fit2$mu.coefficients[2]*gender + fit2$mu.coefficients[3]*age)
sigma1_all <- links$sigma.linkinv(fit1$sigma.coefficients[1])
sigma2_all <- links$sigma.linkinv(fit2$sigma.coefficients[1])
nu1_all <- links$nu.linkinv(fit1$nu.coefficients[1])
nu2_all <- links$nu.linkinv(fit2$nu.coefficients[1])
tau1_all <- links$tau.linkinv(fit1$tau.coefficients[1])
tau2_all <- links$tau.linkinv(fit2$tau.coefficients[1])

# Randomized PIT for discrete margins
dMargin1 <- runif(length(gamma_c_mu1),
  pZISICHEL(gamma_c_mu1 - 1, mu1_all, sigma1_all, nu1_all, tau1_all),
  pZISICHEL(gamma_c_mu1, mu1_all, sigma1_all, nu1_all, tau1_all)
)
dMargin2 <- runif(length(gamma_c_mu2),
  pZISICHEL(gamma_c_mu2 - 1, mu2_all, sigma2_all, nu2_all, tau2_all),
  pZISICHEL(gamma_c_mu2, mu2_all, sigma2_all, nu2_all, tau2_all)
)

library(VineCopula)
copFits<-BiCopSelect(dMargin1,dMargin2,familyset=c(1,3,4,5,6),rotations=FALSE)

-2*logLik(fit1)[1]+(-2*logLik(fit2)[1]) + -2*copFits$logLik 

copula_model<-BiCopEst(dMargin1,dMargin2,family=copFits$family)

ll_m1<--1*fit1$G.deviance/2
ll_m2<--1*fit2$G.deviance/2
ll_cop<-copula_model$logLik

ll_combined<-ll_m1+ll_m2+ll_cop
df_fit<-fit1$df.fit+fit2$df.fit+1

c(ll_combined,-2*ll_combined+2*df_fit,-2*ll_combined+4*df_fit,-2*ll_combined+(log(nrow(dataset))*df_fit))

optim_cop_like <- function(fit1,fit2,copula_model) {
  # Simulate from the fitted copula and margins to get new data
  bicop_sim<-BiCopSim(length(gamma_c_mu1),copula_model)

  qMargin1<-qZISICHEL(bicop_sim[,1],links$mu.linkinv(fit1$mu.coefficients[1]+fit1$mu.coefficients[2]*gender+fit1$mu.coefficients[3]*age)
  ,links$sigma.linkinv(fit1$sigma.coefficients[1])
  ,links$nu.linkinv(fit1$nu.coefficients[1])
  ,links$tau.linkinv(fit1$tau.coefficients[1]))
  qMargin2<-qZISICHEL(bicop_sim[,2],links$mu.linkinv(fit2$mu.coefficients[1]+fit2$mu.coefficients[2]*gender+fit2$mu.coefficients[3]*age)
  ,links$sigma.linkinv(fit2$sigma.coefficients[1])
  ,links$nu.linkinv(fit2$nu.coefficients[1])
  ,links$tau.linkinv(fit2$tau.coefficients[1]))
  
  # Fit models to the SIMULATED data instead of original data
  fit1<-gamlss(qMargin1~(gender)+age,family=ZISICHEL(),method=RS(100))
  fit2<-gamlss(qMargin2~(gender)+age,family=ZISICHEL(),method=RS(100))

  # Return negative log-likelihood for minimization
  return(list(c(coef(fit1),coef(fit2)),
    -(logLik(fit1)[1]+logLik(fit2)[1])))
}

# Estimate standard error for beta_t via simulation of fitted model

mean_margins=matrix(NA,ncol=6,nrow=0)

for (i in 1:100) {
  #i=1
  set.seed(1000+i) # Set different seed for each run to get different simulations
  print(i)
  optim_fit <- optim_cop_like(fit1,fit2,copula_model)
  mean_margins <- rbind(mean_margins, c(optim_fit[[1]],optim_fit[[2]]))
}

par(mfrow=c(2,3))
for (i in 1:ncol(mean_margins)) {
  hist(mean_margins[,i],main=colnames(mean_margins)[i])
}

#colMeans(mean_margins)
cov_sims=cov(mean_margins)

coef_results_se=rbind(model_re_nosig$mu.coefficients[1:4]
      , summary(model_re_nosig)[1:4,2]
      , c(fit1$mu.coefficients[1],fit2$mu.coefficients[1],(fit1$mu.coefficients[2:3]+fit2$mu.coefficients[2:3])/2)
      , c(sqrt(diag(cov_sims))[c(1,4)],sqrt(cov_sims[2,2]+cov_sims[5,5]+2*cov_sims[2,5])
      ,sqrt(cov_sims[3,3]+cov_sims[6,6]+2*cov_sims[3,6])*10)
)

rownames(coef_results_se)=c("ZIS GAMLSS", "ZIS GAMLSS SE", "ZIS GJRM", "ZIS GJRM SE")
coef_results_se

########## DIAGNOSTIC: Check ZIS GAMLSS (model_re_nosig) simulation setup #########
cat("\n=== DIAGNOSTIC: ZIS GAMLSS (model_re_nosig) Simulation Parameters ===\n")
cat("\nCoefficients being passed to sim_model_zis:\n")
cat("  mu.coefficients: c(", paste(round(model_re_nosig$mu.coefficients, 4), collapse=", "), ")\n")
cat("    Length:", length(model_re_nosig$mu.coefficients), "\n")
cat("  sigma.coefficients: c(", paste(round(model_re_nosig$sigma.coefficients, 4), collapse=", "), ")\n")
cat("    Length:", length(model_re_nosig$sigma.coefficients), "\n")
cat("  nu.coefficients: c(", paste(round(model_re_nosig$nu.coefficients, 4), collapse=", "), ")\n")
cat("    Length:", length(model_re_nosig$nu.coefficients), "\n")
cat("  tau.coefficients: c(", paste(round(model_re_nosig$tau.coefficients, 4), collapse=", "), ")\n")
cat("    Length:", length(model_re_nosig$tau.coefficients), "\n")
cat("  Random effect variance:", round(as.numeric(VarCorr(getSmo(model_re_nosig))[[1]]), 6), "\n")

cat("\n=== DIAGNOSTIC: ZIS GJRM Simulation Parameters ===\n")
cat("\nCoefficients being passed to sim_model_zis:\n")
cat("  fit1 mu: c(", paste(round(fit1$mu.coefficients, 4), collapse=", "), ")\n")
cat("  fit2 mu: c(", paste(round(fit2$mu.coefficients, 4), collapse=", "), ")\n")
cat("  fit1 sigma: ", round(fit1$sigma.coefficients[1], 4), " [intercept only]\n")
cat("  fit2 sigma: ", round(fit2$sigma.coefficients[1], 4), " [intercept only]\n")
cat("  fit1 nu: ", round(fit1$nu.coefficients[1], 4), " [intercept only]\n")
cat("  fit2 nu: ", round(fit2$nu.coefficients[1], 4), " [intercept only]\n")
cat("  fit1 tau: ", round(fit1$tau.coefficients[1], 4), " [intercept only]\n")
cat("  fit2 tau: ", round(fit2$tau.coefficients[1], 4), " [intercept only]\n")
# Test manual simulation for first few observations
cat("\n=== Testing manual simulation for first 3 observations ===\n")
set.seed(1000)
test_simCop <- BiCopSim(N=nrow(dataset)/2, family=copula_model$family, par=copula_model$par)
for(i in 1:3) {
  lp_mu1 <- fit1$mu.coefficients[1] + fit1$mu.coefficients[2]*gender[i] + fit1$mu.coefficients[3]*age[i]
  lp_mu2 <- fit2$mu.coefficients[1] + fit2$mu.coefficients[2]*gender[i] + fit2$mu.coefficients[3]*age[i]
  
  lp_sigma1 <- fit1$sigma.coefficients[1]
  lp_sigma2 <- fit2$sigma.coefficients[1]
  
  lp_nu1 <- fit1$nu.coefficients[1]
  lp_nu2 <- fit2$nu.coefficients[1]
  
  lp_tau1 <- fit1$tau.coefficients[1]
  lp_tau2 <- fit2$tau.coefficients[1]
  
  cat(paste0("\nObs ", i, " (gender=", gender[i], ", age=", age[i], "):\n"))
  cat(paste0("  Copula values: u1=", round(test_simCop[i,1], 4), ", u2=", round(test_simCop[i,2], 4), "\n"))
  cat(paste0("  Transformed mu: mu1=", round(links$mu.linkinv(lp_mu1), 4), ", mu2=", round(links$mu.linkinv(lp_mu2), 4), "\n"))
  cat(paste0("  Transformed sigma: sigma1=", round(links$sigma.linkinv(lp_sigma1), 4), ", sigma2=", round(links$sigma.linkinv(lp_sigma2), 4), "\n"))
}

########## CALC VARIOGRAM SCORES MANUALLY #########

set.seed(1000)
source("common_functions.R")
library(scoringRules)

coefficients=results$coefficients
sigmas=results$sigmas
correlations=results$correlations
dist="PO"
n=nrow(dataset)/2

rownames(correlations)=rownames(sigmas)

model_list=rownames(correlations)
model_list=c(model_list,"ZIS GAMLSS", "ZIS GJRM")
model_list_complete=model_list
vg_sims=100; n_runs=100

# Pre-compute weight matrices (constant across runs)
w_vs=matrix(1,ncol=nrow(dataset),nrow=nrow(dataset))
w_vs_0=matrix(1,ncol=nrow(dataset),nrow=nrow(dataset))
for (i in 1:nrow(w_vs)) {
        for (j in 1:ncol(w_vs)) {
          if (abs(i-j)==n) {
            w_vs[i,j]=4*nrow(dataset)
            w_vs_0[i,j]=1
          }
        }
      }      

# Storage: one row per run, one column per model
vs2_all=matrix(NA, nrow=n_runs, ncol=length(model_list))
vs2_wt_all=matrix(NA, nrow=n_runs, ncol=length(model_list))
colnames(vs2_all)=colnames(vs2_wt_all)=model_list

for(run in 1:n_runs) {
  cat("=== Run", run, "of", n_runs, "===\n")
  
  # Simulate from each model
  sim_model_out=list()
  model_list_run=model_list
  for(model in model_list) {
    print(model)
    sim_model_out[[model]]=matrix(NA,ncol=vg_sims,nrow=nrow(dataset))
    for(i in 1:vg_sims) {
      print(i)
      set.seed(i+1000*(run)) # Ensure different seed for each sim and run, but reproducible
      if(model=="ZIS GAMLSS") {
        ## GAMLSS SIM
        sim_model_out[[model]][,i]=
        sim_model_zis(model="re_nosig",dist="ZISICHEL",n=n
          ,coefficients=model_re_nosig$mu.coefficients,sigmas=model_re_nosig$sigma.coefficients
          ,nus=model_re_nosig$nu.coefficients,taus=model_re_nosig$tau.coefficients
          ,correlations=(getSmo(model_re_nosig)$sigb) #getSmo(model_re_nosig)$sigb
          ,age=age,sex=gender)
      } else if (model=="ZIS GJRM") {
        ## ZIS GJRM SIM
        sim_model_out[[model]][,i]=sim_model_zis(model="cop_f",dist="ZISICHEL",n=n
          , coefficients=c(fit1$mu.coefficients, fit2$mu.coefficients)
          , sigmas=c(fit1$sigma.coefficients[1], fit2$sigma.coefficients[1])
          , nus=c(fit1$nu.coefficients[1], fit2$nu.coefficients[1])
          , taus=c(fit1$tau.coefficients[1], fit2$tau.coefficients[1])
          , correlations= c(copula_model$par),age=age,sex=gender)
      } else {
        sim_model_out[[model]][,i]=sim_model(model=model,dist=dist,n=n,coefficients=coefficients,sigmas=sigmas,correlations=correlations,
        age=age, sex=gender)
      }
    }
    if(all(is.na(sim_model_out[[model]]))) {
      model_list_run=model_list_run[model_list_run!=model]
    }
  }
  
  # Calculate Variogram scores for this run (with p=2)
  for (model in model_list_run) {
    print(model)
    vs2_wt_all[run, model]= vs_sample(y=c(gamma_c_mu1,gamma_c_mu2),dat=sim_model_out[[model]],p=2,w_vs=w_vs)
    vs2_all[run, model]=    vs_sample(y=c(gamma_c_mu1,gamma_c_mu2),dat=sim_model_out[[model]],p=2)
    }
}

# Keep only models that had valid results in all runs
#valid_models=colnames(vs2_all)[colSums(is.na(vs2_all))==0&&!c("ZIS GAMLSS", "ZIS GJRM")]
valid_models=colnames(vs2_all)[colSums(is.na(vs2_all))==0]

# Plot boxplots with ggplot2
library(ggplot2)
library(tidyr)
library(gridExtra)

rename_model <- function(x) {
  main_map <- c(
    glm = "GLM",
    gee = "GEE",
    re_nosig = "GAMLSS",
    re_np = "GAMLSS NP",
    lme4 = "LME4",
    gamm = "GAMM",
    "ZIS GJRM" = "GJRM (ZIS-F)",
    "ZIS GAMLSS" = "GAMLSS (ZIS)"
  )
  cop_map <- c(
    cop = "GJRM (C)",
    cop_n = "GJRM (N)",
    cop_j = "GJRM (J)",
    cop_g = "GJRM (G)",
    cop_f = "GJRM (F)",
    cop_amh = "GJRM (AMH)",
    cop_fgm = "GJRM (FGM)",
    cop_pl = "GJRM (PL)",
    cop_h = "GJRM (H)",
    cop_t = "GJRM (T)"
  )
  if (x %in% names(main_map)) return(main_map[x])
  if (x %in% names(cop_map)) return(cop_map[x])
  return(x)
}

valid_model_labels <- vapply(valid_models, rename_model, character(1))

# Prepare data for all VS scores
vs2_df=as.data.frame((vs2_all))
vs2_df$run=1:nrow(vs2_df)
vs2_long=pivot_longer(vs2_df, cols=-run, names_to="Model", values_to="Score")
vs2_long=vs2_long[vs2_long$Model %in% valid_models,]
vs2_long$Model=factor(vs2_long$Model, levels=valid_models)
vs2_long$ModelLabel=factor(vapply(as.character(vs2_long$Model), rename_model, character(1)), levels=valid_model_labels)

vs2_wt_df=as.data.frame((vs2_wt_all))
vs2_wt_df$run=1:nrow(vs2_wt_df)
vs2_wt_long=pivot_longer(vs2_wt_df, cols=-run, names_to="Model", values_to="Score")
vs2_wt_long=vs2_wt_long[vs2_wt_long$Model %in% valid_models,]
vs2_wt_long$Model=factor(vs2_wt_long$Model, levels=valid_models)
vs2_wt_long$ModelLabel=factor(vapply(as.character(vs2_wt_long$Model), rename_model, character(1)), levels=valid_model_labels)

# VS p=2 plots
min_median_vs2 <- min(tapply(log(vs2_long$Score), vs2_long$Model, median, na.rm=TRUE))
p1=ggplot(vs2_long, aes(x=ModelLabel, y=log(Score))) +
  geom_hline(yintercept = min_median_vs2, linetype="dashed", color="red", linewidth=1) +
  stat_summary(fun = median, geom = "point", size = 3, color = "darkblue") +
  stat_summary(fun.data = function(x) {
    q = quantile(x, c(0.025, 0.975), na.rm = TRUE)
    data.frame(y = median(x), ymin = q[1], ymax = q[2])
  }, geom = "linerange", width = 0.2, color = "darkblue") +
  labs(title="VS (p=2, Unweighted)", x=NULL, y="Variogram Score") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

min_median_vs2_wt <- min(tapply(log(vs2_wt_long$Score), vs2_wt_long$Model, median, na.rm=TRUE))
p2=ggplot(vs2_wt_long, aes(x=ModelLabel, y=log(Score))) +
  geom_hline(yintercept = min_median_vs2_wt, linetype="dashed", color="red", linewidth=1) +
  stat_summary(fun = median, geom = "point", size = 3, color = "darkred") +
  stat_summary(fun.data = function(x) {
    q = quantile(x, c(0.025, 0.975), na.rm = TRUE)
    data.frame(y = median(x), ymin = q[1], ymax = q[2])
  }, geom = "linerange", width = 0.2, color = "darkred") +
  labs(title="VS (p=2, Weighted)", x=NULL, y="Variogram Score") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
variogram_plot_width <- 8.5
variogram_plot_height <- 6.6



########## ADDITIONAL FILTERED PLOTS (excluding ZIS GJRM and lme4) #########
# Filter data to exclude problematic models
library(dplyr)

exclude_models <- c("lme4", "gamm", "ZIS GAMLSS", "re_nosig")
filtered_model_labels <- vapply(setdiff(valid_models, exclude_models), rename_model, character(1))

vs2_long_filt <- vs2_long %>% filter(!Model %in% exclude_models)
vs2_long_filt$Model <- factor(vs2_long_filt$Model, levels=setdiff(levels(vs2_long_filt$Model), exclude_models))
vs2_long_filt$ModelLabel <- factor(vapply(as.character(vs2_long_filt$Model), rename_model, character(1)), levels=filtered_model_labels)

vs2_wt_long_filt <- vs2_wt_long %>% filter(!Model %in% exclude_models)
vs2_wt_long_filt$Model <- factor(vs2_wt_long_filt$Model, levels=setdiff(levels(vs2_wt_long_filt$Model), exclude_models))
vs2_wt_long_filt$ModelLabel <- factor(vapply(as.character(vs2_wt_long_filt$Model), rename_model, character(1)), levels=filtered_model_labels)

# VS p=2 plots - filtered
min_median_vs2_filt <- min(tapply(log(vs2_long_filt$Score), vs2_long_filt$Model, median, na.rm=TRUE))
p1_filt=ggplot(vs2_long_filt, aes(x=ModelLabel, y=log(Score))) +
  geom_hline(yintercept = min_median_vs2_filt, linetype="dashed", color="red", linewidth=1) +
  stat_summary(fun = median, geom = "point", size = 3, color = "darkblue") +
  stat_summary(fun.data = function(x) {
    q = quantile(x, c(0.025, 0.975), na.rm = TRUE)
    data.frame(y = median(x), ymin = q[1], ymax = q[2])
  }, geom = "linerange", width = 0.2, color = "darkblue") +
  labs(title="VS (p=2, Unweighted) - Filtered", x=NULL, y="Variogram Score") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

min_median_vs2_wt_filt <- min(tapply(log(vs2_wt_long_filt$Score), vs2_wt_long_filt$Model, median, na.rm=TRUE))
p2_filt=ggplot(vs2_wt_long_filt, aes(x=ModelLabel, y=log(Score))) +
  geom_hline(yintercept = min_median_vs2_wt_filt, linetype="dashed", color="red", linewidth=1) +
  stat_summary(fun = median, geom = "point", size = 3, color = "darkred") +
  stat_summary(fun.data = function(x) {
    q = quantile(x, c(0.025, 0.975), na.rm = TRUE)
    data.frame(y = median(x), ymin = q[1], ymax = q[2])
  }, geom = "linerange", width = 0.2, color = "darkred") +
  labs(title="VS (p=2, Weighted) - Filtered", x=NULL, y="Variogram Score") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

vs_main_grid <- arrangeGrob(p1, p2, p1_filt, p2_filt, ncol=2)
grid.arrange(vs_main_grid)

dir.create("Charts", showWarnings = FALSE)
ggsave(
  filename = "Charts/variogram_scores_all.png",
  plot = vs_main_grid,
  width = variogram_plot_width,
  height = variogram_plot_height,
  units = "in",
  dpi = 300
)

########## DIAGNOSTIC: Compare distributions #########
cat("\n=== DIAGNOSTIC: Distribution Comparison ===\n")
cat("\nOriginal data (all time points):\n")
cat("  Mean:", round(mean(dataset$random_variable), 4), "\n")
cat("  SD:", round(sd(dataset$random_variable), 4), "\n")
cat("  Skewness:", round(skewness(dataset$random_variable), 4), "\n")
cat("  N zeros:", sum(dataset$random_variable==0), "/", length(dataset$random_variable), "\n")
cat("  Prop zeros:", round(mean(dataset$random_variable==0), 4), "\n")

cat("\nOriginal time 1 (gamma_c_mu1):\n")
cat("  Mean:", round(mean(gamma_c_mu1), 4), "\n")
cat("  SD:", round(sd(gamma_c_mu1), 4), "\n")
cat("  N zeros:", sum(gamma_c_mu1==0), "/", length(gamma_c_mu1), "\n")

cat("\nOriginal time 2 (gamma_c_mu2):\n")
cat("  Mean:", round(mean(gamma_c_mu2), 4), "\n")
cat("  SD:", round(sd(gamma_c_mu2), 4), "\n")
cat("  N zeros:", sum(gamma_c_mu2==0), "/", length(gamma_c_mu2), "\n")

cat("\nZIS GJRM simulated (all, first sim):\n")
cat("  Mean:", round(mean(sim_model_out[["ZIS GJRM"]][,1]), 4), "\n")
cat("  SD:", round(sd(sim_model_out[["ZIS GJRM"]][,1]), 4), "\n")
cat("  Skewness:", round(skewness(sim_model_out[["ZIS GJRM"]][,1]), 4), "\n")
cat("  N zeros:", sum(sim_model_out[["ZIS GJRM"]][,1]==0), "/", length(sim_model_out[["ZIS GJRM"]][,1]), "\n")
cat("  Prop zeros:", round(mean(sim_model_out[["ZIS GJRM"]][,1]==0), 4), "\n")

cat("\nZIS GJRM simulated time 1 (first n obs):\n")
cat("  Mean:", round(mean(sim_model_out[["ZIS GJRM"]][1:n,1]), 4), "\n")
cat("  SD:", round(sd(sim_model_out[["ZIS GJRM"]][1:n,1]), 4), "\n")
cat("  N zeros:", sum(sim_model_out[["ZIS GJRM"]][1:n,1]==0), "/", n, "\n")

cat("\nZIS GJRM simulated time 2 (last n obs):\n")
cat("  Mean:", round(mean(sim_model_out[["ZIS GJRM"]][(n+1):(2*n),1]), 4), "\n")
cat("  SD:", round(sd(sim_model_out[["ZIS GJRM"]][(n+1):(2*n),1]), 4), "\n")
cat("  N zeros:", sum(sim_model_out[["ZIS GJRM"]][(n+1):(2*n),1]==0), "/", n, "\n")

# EXTREME VALUE DIAGNOSTIC
cat("\n=== EXTREME VALUE CHECK (ZIS GJRM) ===\n")
zis_sim <- sim_model_out[["ZIS GJRM"]][,1]
cat("Max value simulated:", max(zis_sim), "\n")
cat("95th percentile:", quantile(zis_sim, 0.95), "\n")
cat("99th percentile:", quantile(zis_sim, 0.99), "\n")
cat("N values > 50:", sum(zis_sim > 50), "/", length(zis_sim), "\n")
cat("N values > 100:", sum(zis_sim > 100), "/", length(zis_sim), "\n")

cat("\nOriginal data extremes:\n")
orig_data <- c(gamma_c_mu1, gamma_c_mu2)
cat("Max value:", max(orig_data), "\n")
cat("95th percentile:", quantile(orig_data, 0.95), "\n")
cat("99th percentile:", quantile(orig_data, 0.99), "\n")
cat("N values > 50:", sum(orig_data > 50), "/", length(orig_data), "\n")
cat("N values > 100:", sum(orig_data > 100), "/", length(orig_data), "\n")

cat("\n** If simulated data has many more extreme values, that explains the overdispersion **\n")

cat("\n=== CORRELATION STRUCTURE CHECK ===\n")
cat("Original correlation (Kendall):", round(cor(gamma_c_mu1, gamma_c_mu2, method="kendall"), 4), "\n")
cat("Original correlation (Pearson):", round(cor(gamma_c_mu1, gamma_c_mu2, method="pearson"), 4), "\n")

n=length(gamma_c_mu1)
sim_t1 <- sim_model_out[["ZIS GJRM"]][1:n, 1]
sim_t2 <- sim_model_out[["ZIS GJRM"]][(n+1):(2*n), 1]
cat("\nZIS GJRM simulated correlation (Kendall):", round(cor(sim_t1, sim_t2, method="kendall"), 4), "\n")
cat("ZIS GJRM simulated correlation (Pearson):", round(cor(sim_t1, sim_t2, method="pearson"), 4), "\n")

glm_t1 <- sim_model_out[["glm"]][1:n, 1]
glm_t2 <- sim_model_out[["glm"]][(n+1):(2*n), 1]
cat("\nGLM simulated correlation (Kendall):", round(cor(glm_t1, glm_t2, method="kendall"), 4), "\n")
cat("GLM simulated correlation (Pearson):", round(cor(glm_t1, glm_t2, method="pearson"), 4), "\n")

cat("\n** If ZIS GJRM correlation is similar to original, the issue is ONLY variance **\n")
cat("** If correlation is also wrong, both variance AND copula need fixing **\n")

cat("\n=== Checking if fitted model parameters match expectations ===\n")
test_indices <- 1:min(10, length(gender))
cat("For first", length(test_indices), "observations (time 1):\n")
cat("  Expected mean mu1 range: [", round(min(links$mu.linkinv(fit1$mu.coefficients[1] + fit1$mu.coefficients[2]*gender[test_indices] + fit1$mu.coefficients[3]*age[test_indices])), 4), 
  ",", round(max(links$mu.linkinv(fit1$mu.coefficients[1] + fit1$mu.coefficients[2]*gender[test_indices] + fit1$mu.coefficients[3]*age[test_indices])), 4), "]\n")
cat("For time 2 observations:\n")
cat("  Expected mean mu2 range: [", round(min(links$mu.linkinv(fit2$mu.coefficients[1] + fit2$mu.coefficients[2]*gender[test_indices] + fit2$mu.coefficients[3]*age[test_indices])), 4),
  ",", round(max(links$mu.linkinv(fit2$mu.coefficients[1] + fit2$mu.coefficients[2]*gender[test_indices] + fit2$mu.coefficients[3]*age[test_indices])), 4), "]\n")

cat("\n=== OVERDISPERSION DIAGNOSTIC ===\n")
cat("\nChecking ZISICHEL parameters (all observations):\n")
all_mu1 <- links$mu.linkinv(fit1$mu.coefficients[1] + fit1$mu.coefficients[2]*gender + fit1$mu.coefficients[3]*age)
all_mu2 <- links$mu.linkinv(fit2$mu.coefficients[1] + fit2$mu.coefficients[2]*gender + fit2$mu.coefficients[3]*age)
all_sigma1 <- rep(links$sigma.linkinv(fit1$sigma.coefficients[1]), length(gender))
all_sigma2 <- rep(links$sigma.linkinv(fit2$sigma.coefficients[1]), length(gender))
all_nu1 <- rep(links$nu.linkinv(fit1$nu.coefficients[1]), length(gender))
all_nu2 <- rep(links$nu.linkinv(fit2$nu.coefficients[1]), length(gender))
all_tau1 <- rep(links$tau.linkinv(fit1$tau.coefficients[1]), length(gender))
all_tau2 <- rep(links$tau.linkinv(fit2$tau.coefficients[1]), length(gender))

cat("  Time 1 mu range: [", round(min(all_mu1), 4), ",", round(max(all_mu1), 4), "], mean:", round(mean(all_mu1), 4), "\n")
cat("  Time 2 mu range: [", round(min(all_mu2), 4), ",", round(max(all_mu2), 4), "], mean:", round(mean(all_mu2), 4), "\n")
cat("  Time 1 sigma: ", round(all_sigma1[1], 6), " (constant)\n")
cat("  Time 2 sigma: ", round(all_sigma2[1], 6), " (constant)\n")
cat("  Time 1 nu: ", round(all_nu1[1], 6), " (constant)\n")
cat("  Time 2 nu: ", round(all_nu2[1], 6), " (constant)\n")
cat("  Time 1 tau: ", round(all_tau1[1], 4), " (constant)\n")
cat("  Time 2 tau: ", round(all_tau2[1], 4), " (constant)\n")

cat("\n** Check for extreme parameter values that could cause overdispersion **\n")
cat("Sigma values close to 0 or >1 can cause issues\n")
cat("Nu values far from -1 can cause heavy tails\n")
cat("Tau (zero-inflation) should match observed proportion of zeros\n")

# Check correlation strength
cat("\nCopula parameter:", copula_model$par, "(family:", copula_model$family, ")\n")
cat("Empirical Kendall's tau:", round(cor(gamma_c_mu1, gamma_c_mu2, method="kendall"), 4), "\n")
implied_tau <- BiCopPar2Tau(copula_model$family, copula_model$par)
cat("Copula implied Kendall's tau:", round(implied_tau, 4), "\n")

cat("\nFit1 (time 1) coefficients:\n")
cat("  Mu: intercept=", round(fit1$mu.coefficients[1], 4), ", sex=", round(fit1$mu.coefficients[2], 4), ", age=", round(fit1$mu.coefficients[3], 4), "\n")
cat("  Sigma: ", round(fit1$sigma.coefficients[1], 4), " (intercept only)\n")
cat("  Nu: ", round(fit1$nu.coefficients[1], 4), " (intercept only)\n")
cat("  Tau: ", round(fit1$tau.coefficients[1], 4), " (intercept only)\n")
cat("\nFit2 (time 2) coefficients:\n")
cat("  Mu: intercept=", round(fit2$mu.coefficients[1], 4), ", sex=", round(fit2$mu.coefficients[2], 4), ", age=", round(fit2$mu.coefficients[3], 4), "\n")
cat("  Sigma: ", round(fit2$sigma.coefficients[1], 4), " (intercept only)\n")
cat("  Nu: ", round(fit2$nu.coefficients[1], 4), " (intercept only)\n")
cat("  Tau: ", round(fit2$tau.coefficients[1], 4), " (intercept only)\n")

rbind(
c(mean(dataset$random_variable),sd(dataset$random_variable),skewness(dataset$random_variable),sum(as.numeric(dataset$random_variable==0)))
,c(mean(sim_model_out[["ZIS GJRM"]][,1]),sd(sim_model_out[["ZIS GJRM"]][,1]),skewness(sim_model_out[["ZIS GJRM"]][,1]),sum(as.numeric(sim_model_out[["ZIS GJRM"]][,1]==0)))
,c(mean(sim_model_out[["ZIS GAMLSS"]][,1]),sd(sim_model_out[["ZIS GAMLSS"]][,1]),skewness(sim_model_out[["ZIS GAMLSS"]][,1]),sum(as.numeric(sim_model_out[["ZIS GAMLSS"]][,1]==0)))
,c(mean(sim_model_out[["glm"]][,1]),sd(sim_model_out[["glm"]][,1]),skewness(sim_model_out[["glm"]][,1]),sum(as.numeric(sim_model_out[["glm"]][,1]==0)))
)
########## VS Summary #########

cat("\n=== VARIOGRAM SCORES SUMMARY ===\n")
cat("\nVS (p=2) - Unweighted:\n")
#print(vs2_all)
cat("\nVS (p=2) - Weighted:\n")
#print(vs2_wt_all)

cat("\n=== LOG VARIOGRAM SCORES ===\n")
log_vs_matrix <- log(t(rbind((vs2_all),(vs2_wt_all))))
colnames(log_vs_matrix) <- c("VS(p=2)", "VS(p=2,wt)")

#print(log_vs_matrix)

# Create summary table with mean, median, and 95% range for each model
cat("\n=== VARIOGRAM SCORES SUMMARY TABLE ===\n")

# Function to compute summary statistics
compute_summary <- function(scores) {
  c(
    Mean = mean(scores, na.rm=TRUE),
    Median = median(scores, na.rm=TRUE),
    Q2.5 = quantile(scores, 0.025, na.rm=TRUE),
    Q97.5 = quantile(scores, 0.975, na.rm=TRUE)
  )
}

# Compute summaries for unweighted scores
log_vs2 <- log(t(vs2_all))
summary_unweighted <- t(apply(log_vs2, 1, compute_summary))
colnames(summary_unweighted) <- c("Mean_VS2", "Median_VS2", "Q2.5_VS2", "Q97.5_VS2")

# Compute summaries for weighted scores
log_vs2_wt <- log(t(vs2_wt_all))
summary_weighted <- t(apply(log_vs2_wt, 1, compute_summary))
colnames(summary_weighted) <- c("Mean_VS2wt", "Median_VS2wt", "Q2.5_VS2wt", "Q97.5_VS2wt")

# Combine into one table and add readable model names
summary_table <- cbind(summary_unweighted, summary_weighted)
rownames(summary_table) <- vapply(rownames(summary_table), rename_model, character(1))

# Round for readability
summary_table_rounded <- round(summary_table, 4)

print(summary_table_rounded)
