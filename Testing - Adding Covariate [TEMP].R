# OVERVIEW: Run any single simulation and plot the results of the fitted models against the estimated true values.
# WARNING: You must restart R before each run because loading the GJRM function breaks GAMLSS

######### INPUT PARAMETERS ##########

dist="NO" # Distribution to use for the simulation
sims=100 # Number of simulations to run for the true value calculation
n=1000 # Number of observations in the dataset (note there are two timepoints so this will create n x 2 data points)
a=1 # Parameter a for the distribution
b=1 # Parameter b for the distribution
c=.5 # Parameter c for the distribution (not used for GA)
mu1=1 # Mean for time 1
mu2=2 # Mean for time 2
x1=1 # Covariate 1 value
x2=1 # Covariate 2 value
# Set the seed for reproducibility
set.seed(1000)

############## Generate, Fit Models, and Calculate True Values ##############

# Load the necessary libraries
source("common_functions.R")
# Generate a dataset with covariates
dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
# Fit all models to the dataset
results=fitBivModels_Bt_withCov(dataset,dist=dist,include="ALL",a=a,b=b,c=c,mu1=mu1,mu2=mu2,calc_actuals=FALSE)
# Calculate the MLE estimate to find true value for the coefficients and standard errors based on simulation
true=simCovariateMLEs(sims=sims,n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2,trace=TRUE)

############## Plot Models ##############
df=results$coefficients
plot.new(); par(mfrow=c(ncol(df),2),las=3)
true_vals= true$coefficients
factor_names=c("time 1 intercept","time 2 intercept","binary","linear")
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", ylab = "Est",xlab="",main= paste("Est |",factor_names[i]),ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}

df=results$ses
true_vals = true$ses
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", ylab = "SE", xlab="",main= paste("SE |",factor_names[i]),ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}



################################## VARIGORAM METHOD TESTING ####
source("common_functions.R")
#sims=100
#dist="GA" ;sims=100 ;n=1000 ;a=.5 ;b=2 ;c=NA ;mu1=1 ;mu2=2 ;x1=1 ;x2=1 # Covariate 2 value
#dist="PO"; a=NA;b=1;c=2;mu1=1;mu2=2;x1=1;x2=1;n=1000
#dist="NO"; a=1;b=1;c=.75;mu1=.25;mu2=.75;x1=1;x2=1;n=1000
outer_runs=100
vs_outer=vs_outer_wt=es_outer=vs1_outer=matrix(NA, ncol=4,nrow=outer_runs)

for (outer_run in 1:outer_runs) {
  
  print(outer_run)
  dist="NO";a=1;b=1;c=.75;mu1=1;mu2=2;x1=1;x2=1;n=100
  # Generate a dataset with covariates
  dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
  # Fit all models to the dataset
  #results=fitBivModels_Bt_withCov(dataset,dist=dist,include="ALL",a=a,b=b,c=c,mu1=mu1,mu2=mu2,calc_actuals=FALSE)
  # Calculate the MLE estimate to find true value for the coefficients and standard errors based on simulation
  #true=simCovariateMLEs(sims=sims,n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2,trace=TRUE)
  
  dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
  #dataset_test=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2); dataset_test$patient=dataset_test$patient+1000
  #true=simCovariateMLEs(sims=sims,n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2,trace=TRUE)
  
  library("gamlss2")
  require(gamlss)
  require(gee)
  require(lme4)
  require(MASS)
  require(gamlss.mx)
  library(mgcv)
  library(glmtoolbox)
  
  if(dist=="GA") {
    invisible(capture.output(model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=GA)))
    invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=Gamma(link = "log"), maxit=1000)))
    invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
    invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                     , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
    invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset, family=Gamma(link="log"))))
    invisible(capture.output(model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=Gamma(link="log"))))
  } else if(dist=="NO") {
    invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=gaussian, maxit=1000)))
    invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
    #invisible(capture.output(model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NO)))
    invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NO)))
    
    invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                     , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
    invisible(capture.output(model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)))
    invisible(capture.output(model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=gaussian)))
  } else if(dist=="PO") {
    invisible(capture.output(model_glm <- glm.nb(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, maxit=1000)))
    #invisible(capture.output(model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
    #model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, init.beta=model_glm$coefficients,
    #                  family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")
    invisible(capture.output(model_gee<-overglm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family="nb1(log)", maxiter=25, corstr = "exchangeable")))
    invisible(capture.output(model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI())))
    #invisible(capture.output(model_re <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
    invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                                                     , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
    model_lme4 <- glmer.nb(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)
    model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=nb(link="log"))
  } else if (dist=="LO") {
    invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=binomial, maxit=1000)))
    invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))
    #model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=BI)
    invisible(capture.output(model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI)))
    #invisible(capture.output(model_re <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
    invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                                                     , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
    
    invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) +as.factor(sex)+age+ (1|patient), data=dataset,family=binomial)))
    
    invisible(capture.output(model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=binomial)))
  }

  ################### VARIOGRAM ######################
  # OK to calculate the varigram score for p=2 we can start with the values for the realisations
  #Take pairwise difference for all components squared
  df=dataset
  
  #For each model generate the data X times based on the fitted parameters
  results_table=list()
  results_table[["glm"]]=summary(model_glm)$coeff[,1:2]
  results_table[["gee"]]=summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2]
  results_table[["re_nosig"]]=cbind(summary(model_re_nosig)[1:4],summary(model_re_nosig)[6:9])
  #results_table[["re_np"]]=cbind(summary(model_re_np)[1:4],summary(model_re_np)[8:11])
  results_table[["lme4"]]=summary(model_lme4)$coefficients[,c(1,2)]
  #results_table[["gamm"]]=cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix)))
  
  sigmas=list()
  sigmas[["glm"]]=sqrt(rep(summary(model_glm)$dispersion,2))
  sigmas[["gee"]]=sqrt(rep(model_gee$phi,2))
  sigmas[["re_nosig"]]=sqrt(rep(exp(model_re_nosig$sigma.coefficients),2)^2+as.numeric(VarCorr(getSmo(model_re_nosig))[[1]]))
  #sigmas[["re_np"]]=exp(model_re_np$sigma.coefficients)
  sigmas[["lme4"]]=sqrt(rep(summary(model_lme4)$sigma,2)^2+summary(model_lme4)$varcor$patient[1,1])
  #sigmas[["gamm"]]=rep(model_gamm$lme$sigma,2)
  
  #Extract correlations
  correlations=list()
  correlations[["glm"]]=0
  correlations[["gee"]]=(model_gee$corr[2,1])
  correlations[["re_nosig"]]=as.numeric(VarCorr(getSmo(model_re_nosig))[[1]]) / (sigmas[["re_nosig"]][1]^2)
  #correlations[["re_np"]]=0
  correlations[["lme4"]]=summary(model_lme4)$varcor$patient[1,1] / (sigmas[["lme4"]][1]^2)
  #correlations[["gamm"]]=model_gamm$lme$apVar
  
  model_list=c("glm","gee","re_nosig","lme4")
  
  ###ORIGINAL
  vg_sims=100; vg_x=vg_x_mean=list(); vg=list(); datasets=list()
  
  for(model in model_list) {
    for (run in 1:vg_sims) {
      mu1=results_table[[model]][1,1]
      mu2=results_table[[model]][2,1]
      x1=results_table[[model]][3,1]
      x2=results_table[[model]][4,1]
      
      a=sigmas[[model]][1]
      b=sigmas[[model]][2]
      c=correlations[[model]]
      
      datasets[[model]][[run]]=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)$random_variable
    }
  }
  
  library(scoringRules)
  
  w_vs=matrix(1,ncol=length(df$random_variable),nrow=length(df$random_variable))
  # Set values where row = col + n or col = row + n to 10
  # For every pair of adjacent observations
  for (i in seq(1, 2*n, by=2)) {
    w_vs[i, i+1] <- n
    w_vs[i+1, i] <- n
  }
  
  
  vs_outer_wt[outer_run,1]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["glm"]])),p=2,w_vs=w_vs)
  vs_outer_wt[outer_run,2]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["gee"]])),p=2,w_vs=w_vs)
  vs_outer_wt[outer_run,3]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["re_nosig"]])),p=2,w_vs=w_vs)
  vs_outer_wt[outer_run,4]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["lme4"]])),p=2,w_vs=w_vs)
  
  vs_outer[outer_run,1]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["glm"]])),p=2)
  vs_outer[outer_run,2]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["gee"]])),p=2)
  vs_outer[outer_run,3]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["re_nosig"]])),p=2)
  vs_outer[outer_run,4]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["lme4"]])),p=2)
  
  es_outer[outer_run,1]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["glm"]])))
  es_outer[outer_run,2]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["gee"]])))
  es_outer[outer_run,3]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["re_nosig"]])))
  es_outer[outer_run,4]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["lme4"]])))
  
  vs1_outer[outer_run,1]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["glm"]])),p=1)
  vs1_outer[outer_run,2]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["gee"]])),p=1)
  vs1_outer[outer_run,3]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["re_nosig"]])),p=1)
  vs1_outer[outer_run,4]=vs_sample(y=df$random_variable,dat=as.matrix(do.call(cbind, datasets[["lme4"]])),p=1)
  
}

# Assign model names to columns
colnames(vs_outer)  <- model_list
colnames(es_outer)  <- model_list
colnames(vs1_outer) <- model_list
colnames(vs_outer_wt) <- model_list

# Convert matrices to data frames
vs_df   <- as.data.frame(vs_outer)
es_df   <- as.data.frame(es_outer)
vs1_df  <- as.data.frame(vs1_outer)
vs_wt_df <- as.data.frame(vs_outer_wt)

# Calculate ranks for each row (1 = lowest value, 4 = highest value)
vs_ranks   <- t(apply(vs_df,   1, rank))
es_ranks   <- t(apply(es_df,   1, rank))
vs1_ranks  <- t(apply(vs1_df,  1, rank))
vs_wt_ranks <- t(apply(vs_wt_df, 1, rank))

# Add rank columns to the data frames
rank_colnames <- paste0(model_list, "_rank")
vs_df[rank_colnames]   <- vs_ranks
es_df[rank_colnames]   <- es_ranks
vs1_df[rank_colnames]  <- vs1_ranks
vs_wt_df[rank_colnames] <- vs_wt_ranks

out=rbind(
colMeans(vs_df)
,colMeans(es_df)
,colMeans(vs1_df)
,colMeans(vs_wt_df)
)
rownames(out)=c("VS","ES","VS1","VS_WT")
round(out,0)[,c(1:4)]



