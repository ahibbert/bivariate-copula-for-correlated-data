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
rand <- read_sas("Data/randhrs1992_2020v2.sas7bdat",col_select=c("HHIDPN","R14IWENDY","R15IWENDY"
                                                                 ,"R14DOCTIM" #Doctor visits last 2 years
                                                                 ,"R15DOCTIM" #Doctor visits last 2 years
))

rand_doc_visits <-as.data.frame(rand[!(is.na(rand$R14DOCTIM))&!(is.na(rand$R15DOCTIM))&rand$R14DOCTIM>=0&rand$R15DOCTIM>=0, c("HHIDPN","R14IWENDY","R15IWENDY","R14DOCTIM","R15DOCTIM")])
rand_doc_visits_SAMPLED<-rand_doc_visits[runif(nrow(rand_doc_visits))<0.2,]

#save(rand_doc_visits_SAMPLED,file="rand_doc_visits_SAMPLED")
#load("rand_doc_visits_SAMPLED") ; dist="PO" #Used once the data is sampled once to ensure the same sample across models. 
dist="PO"
gamma_c_mu1<-as.vector(as.data.frame(rand_doc_visits_SAMPLED[,4])$`rand_doc_visits_SAMPLED[, 4]`)
gamma_c_mu2<-as.vector(as.data.frame(rand_doc_visits_SAMPLED[,5])$`rand_doc_visits_SAMPLED[, 5]`)

# Set up as longitudinal structured data - MUST DO THIS STEP TO SETUP APPLICATIONS DATA CORRECTLY FOR MODEL FITTING 
patient<-as.factor(seq(1:length(gamma_c_mu1)))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")
dataset<-dataset[order(dataset$patient),]

a=NA; b=NA; c=NA; mu1=NA; mu2=NA; n=NA #Dummy values to pass to function

####### 2) Simulation ##################

#Bivariate distribution parameters to simulate

#dist="NO";a=1; b=2; c=0.75; mu1=1; mu2=2; n=1000
#dist="GA";a=.25; b=1.75; c=NA; mu1=10; mu2=12; n=1000
dist="GA";a=.2; b=.2; c=NA; mu1=10; mu2=12; n=1000
#dist="PO";a=NA; b=1; c=.1; mu1=5; mu2=5; n=1000 ## Highly skewed
#dist="PO";a=NA; b=.5; c=9; mu1=5; mu2=5; n=1000 ## Not highly skewed
#dist="LO";a=NA; b=NA; c=.5; mu1=.25; mu2=.75; n=1000

# Generate a realisation of the bivariate distribution
dataset <- generateBivDist(n,a,b,c,mu1,mu2,dist)

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
results<-fitBivModels_Bt(data=dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=FALSE)
if(dist=="NO"){clean_results<-results} else 
  if(dist=="LO"){clean_results<-cbind(results,round(logit_inv(results[,c(1,2)]),4));} else
    {clean_results<-cbind(results,round(exp(results[,c(1,2)]),4));}
clean_results<-cbind(clean_results[,c(9,10)],clean_results[,1:8]);colnames(clean_results)<-c("mu_1","mu_2",colnames(clean_results[,c(3:10)]))
rownames(clean_results)<-c("GLM","GEE","GAMLSS (4)","GAMLSS NP (5)","LME4","GAMM","GJRM (Clayton)","GJRM (Normal)","GJRM (Joe)"    ,"GJRM (Gumbel)","GJRM (Frank)" ,"GJRM (AMH)"   ,"GJRM (FGM)"    ,"GJRM (Plackett)","GJRM (Hougaard)","GJRM (T)","Actual")

clean_results_2<-cbind(clean_results,clean_results[,"EDF"])
clean_results_2[,8]<--2*clean_results[,"LogLik"]+2*clean_results[,"EDF"]
clean_results_2[,9]<--2*clean_results[,"LogLik"]+4*clean_results[,"EDF"]
clean_results_2[,10]<-round(-2*clean_results[,"LogLik"]+log(nrow(dataset))*clean_results[,"EDF"])

colnames(clean_results_2)<-c(colnames(clean_results)[1:7],"AIC (2)","AIC (4)","BIC","EDF")
clean_results_2

####### 4) Time to run#############
source("common_functions.R")

#dist="GA";a=.2; b=.2; c=NA; mu1=10; mu2=12; n=1000
#dist="PO";a=NA; b=1; c=.1; mu1=5; mu2=5; n=1000 ## Highly skewed
#dist="PO";a=NA; b=.5; c=9; mu1=5; mu2=5; n=1000 ## Not highly skewed

dist="GA";a=.5; b=1; c=NA; mu1=10; mu2=12; n=1000

# Functions for fitting models and timing runs
timer <- function(dataset) {
  
  require(gamlss)
  require(gee)
  require(gamlss.mx)
  require(lme4)
  require(mgcv)
  times=rep(0,6);i=0
  start=Sys.time()
  
  tryCatch({model_glm <- glm(formula=random_variable~-1+as.factor(time==1), data=dataset, family=Gamma(link="log")) 
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  tryCatch({
  model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  tryCatch({
  model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  tryCatch({
  model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  tryCatch({
    model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=Gamma(link="log"))
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  tryCatch({
  model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                          , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)
  }, error=function(error) {})
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  
  #time[2:length(time)]=time[2:length(time)]-time[1:(length(time)-1)]########ADDED NOT TESTED
  return(times)
}

timerGJRM <- function(dataset) {
  
  require(GJRM)
  
  gamma_c_mu1<-dataset[dataset$time==0,]
  gamma_c_mu2<-dataset[dataset$time==1,]
  
  #Setting up GJRM equations
  eq.mu.1 <- formula(random_variable~1)
  eq.mu.2 <- formula(random_variable.1~1)
  fl <- list(eq.mu.1, eq.mu.2)
  
  margin_dist="GA"
  
  times=rep(0,2)
  start=Sys.time();i=0
  
  model_copula<-    gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
  i=i+1;times[i]=difftime(Sys.time(), start, units = "secs")[[1]]
  
  #time[2:length(time)]=time[2:length(time)]-time[1:(length(time)-1)]########ADDED NOT TESTED

  return(times)
}

# Fitting LM models for each sample size
samplesizes=c(100,500,1000,5000)
sumtime_lms<-matrix(nrow=length(samplesizes)*10,ncol=7)
z=1
for (n in 1:length(samplesizes)) {
  for (i in 1:10) {
    dataset <- generateBivDist(samplesizes[n],a,b,c,mu1,mu2,dist)
    sumtime_lms[z,] = c(n,timer(dataset))
    z=z+1
    print(c(n,i))
  }
}

#Fitting GJRM models for each sample size (this has to be separate as loading GJRM breaks GAMLSS functionality)
sumtimeGJRM<-matrix(nrow=length(samplesizes)*10,ncol=3)
z=1
for (n in 1:length(samplesizes)) {
  for (i in 1:10) {
    tryCatch({
      dataset <- generateBivDist(samplesizes[n],a,b,c,mu1,mu2,dist)
      sumtimeGJRM[z,] = c(n,timerGJRM(dataset))}, error=function(error) {sumtimeGJRM[i,]=c(n,rep(NA,2))})
    z=z+1
    print(c(n,i))
  }
}

# Calculate summary of runtimes
sumtime_lm_diff<-cbind(sumtime_lms[,c(1:2)],sumtime_lms[,3:ncol(sumtime_lms)]-sumtime_lms[,2:(ncol(sumtime_lms)-1)])
sumtime_gjrm_diff<-cbind(sumtimeGJRM[,c(1:2)],sumtimeGJRM[,3:ncol(sumtimeGJRM)]-sumtimeGJRM[,2:(ncol(sumtimeGJRM)-1)])
sumtime_summary<-cbind(sumtime_lm_diff,sumtime_gjrm_diff[,c(2,3)])

sumtime_summary_avg<-rbind(colMeans(sumtime_summary[sumtime_summary[,1]==1,])
                             ,colMeans(sumtime_summary[sumtime_summary[,1]==2,])
                             ,colMeans(sumtime_summary[sumtime_summary[,1]==3,])
                             ,colMeans(sumtime_summary[sumtime_summary[,1]==4,])
                             ,colMeans(sumtime_summary[sumtime_summary[,1]==5,]))

colnames(sumtime_summary_avg)<-c("n","GLM","GEE","GAMLSS (4)","LME4","GAMM","GAMLSS NP (5)","GJRM (C)","GJRM (N)")
sumtime_summary_avg

####### 5) Manual fitting of ZIS for GJRM #############

library(gamlss)

links<-ZISICHEL(mu.link = "log", sigma.link = "log", nu.link = "identity", 
                tau.link = "logit")

# Fit margins
fit1<-gamlss(gamma_c_mu1~1,family=ZISICHEL(),method=RS(100))
fit2<-gamlss(gamma_c_mu2~1,family=ZISICHEL(),method=RS(100))

dMargin1<-pZISICHEL(gamma_c_mu1,links$mu.linkinv(fit1$mu.coefficients),links$sigma.linkinv(fit1$sigma.coefficients),links$nu.linkinv(fit1$nu.coefficients),links$tau.linkinv(fit1$tau.coefficients))
dMargin2<-pZISICHEL(gamma_c_mu2,links$mu.linkinv(fit2$mu.coefficients),links$sigma.linkinv(fit2$sigma.coefficients),links$nu.linkinv(fit2$nu.coefficients),links$tau.linkinv(fit2$tau.coefficients))

par(mfrow=c(1,2))
hist(dMargin1)
hist(dMargin2)

plot(dMargin1,dMargin2)

# Choose and fit best copula
library(VineCopula)
copFits<-BiCopSelect(dMargin1,dMargin2)

copula_model<-BiCopEst(dMargin1,dMargin2,family=copFits$family)

ll_m1<--1*fit1$G.deviance/2
ll_m2<--1*fit2$G.deviance/2
ll_cop<-copula_model$logLik

ll_combined<-ll_m1+ll_m2+ll_cop
df_fit<-fit1$df.fit+fit2$df.fit+2

c(ll_combined,-2*ll_combined+2*df_fit,-2*ll_combined+4*df_fit,-2*ll_combined+(log(nrow(dataset))*df_fit))

# Estimate standard error for beta_t via simulation of fitted model
mean_margins=matrix(NA,ncol=2,nrow=0)
set.seed(100)
for (i in 1:1000) {
  bicop_sim<-BiCopSim(length(gamma_c_mu1),copula_model)
  qMargin1<-qZISICHEL(bicop_sim[,1],links$mu.linkinv(fit1$mu.coefficients),links$sigma.linkinv(fit1$sigma.coefficients),links$nu.linkinv(fit1$nu.coefficients),links$tau.linkinv(fit1$tau.coefficients))
  qMargin2<-qZISICHEL(bicop_sim[,2],links$mu.linkinv(fit2$mu.coefficients),links$sigma.linkinv(fit2$sigma.coefficients),links$nu.linkinv(fit2$nu.coefficients),links$tau.linkinv(fit2$tau.coefficients))
  mean_margins=rbind(mean_margins,c(log(mean(qMargin1)),log(mean(qMargin2))))
}

summary(fit2)[1]-summary(fit1)[1]
sqrt(var(mean_margins[,1])+var(mean_margins[,2])-2*cov(mean_margins)[1,2])
sqrt(vcov(fit1)[1,1]+vcov(fit2)[1,1]-2*cov(mean_margins)[1,2])
