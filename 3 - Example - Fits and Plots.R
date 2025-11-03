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
# results<-fitBivModels_Bt_withCov(data=dataset,dist=dist,include="ALL",a,b,c,mu1,mu2)

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

eval_out=evaluateModels(results,rownames(results$correlations),vg_sims = 100)
eval_out

####### 5) Manual fitting of ZIS #############

library(gamlss)
#GAMLSS MODEL
model_re_nosig <- gamlss(formula=random_variable~ -1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=ZISICHEL(),method=RS(6))

#GJRM MODEL
library(gamlss)
#Test we are working on the same data and model as other fits
#model_re_nosig <- gamlss(formula=random_variable~ as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI(),method=RS(1000))
#model_re_nosig <- gamlss(formula=random_variable~ as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=ZISICHEL(),method=CG(1000))


links<-ZISICHEL(mu.link = "log", sigma.link = "log", nu.link = "identity", 
                tau.link = "logit")

# Fit margins
fit1<-gamlss(gamma_c_mu1~gender+age,family=ZISICHEL(),method=RS(100))
fit2<-gamlss(gamma_c_mu2~gender+age,family=ZISICHEL(),method=RS(100))

dMargin1<-pZISICHEL(gamma_c_mu1,links$mu.linkinv(fit1$mu.coefficients[1]+fit1$mu.coefficients[2]*gender+fit1$mu.coefficients[3]*age),links$sigma.linkinv(fit1$sigma.coefficients),links$nu.linkinv(fit1$nu.coefficients),links$tau.linkinv(fit1$tau.coefficients))
dMargin2<-pZISICHEL(gamma_c_mu2,links$mu.linkinv(fit2$mu.coefficients[1]+fit2$mu.coefficients[2]*gender+fit2$mu.coefficients[3]*age),links$sigma.linkinv(fit2$sigma.coefficients),links$nu.linkinv(fit2$nu.coefficients),links$tau.linkinv(fit2$tau.coefficients))

# Choose and fit best copula
library(VineCopula)
#copFits<-BiCopSelect(dMargin1,dMargin2)
copFits<-BiCopSelect(dMargin1,dMargin2,familyset=c(1,3,4,5,6,13,14))

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

  qMargin1<-qZISICHEL(bicop_sim[,1],links$mu.linkinv(fit1$mu.coefficients[1]+fit1$mu.coefficients[2]*gender+fit1$mu.coefficients[3]*(age))
  ,links$sigma.linkinv(fit1$sigma.coefficients),links$nu.linkinv(fit1$nu.coefficients),links$tau.linkinv(fit1$tau.coefficients))
  qMargin2<-qZISICHEL(bicop_sim[,2],links$mu.linkinv(fit2$mu.coefficients[1]+fit2$mu.coefficients[2]*gender+fit2$mu.coefficients[3]*(age))
  ,links$sigma.linkinv(fit2$sigma.coefficients),links$nu.linkinv(fit2$nu.coefficients),links$tau.linkinv(fit2$tau.coefficients))
  
  fit1<-gamlss(qMargin1~gender+age,family=ZISICHEL(),method=RS(100))
  fit2<-gamlss(qMargin2~gender+age,family=ZISICHEL(),method=RS(100))

  # Return negative log-likelihood for minimization
  return(list(coef(fit1),coef(fit2),
    -(logLik(fit1)[1]+logLik(fit2)[1])))
}

# Estimate standard error for beta_t via simulation of fitted model

mean_margins=matrix(NA,ncol=6,nrow=0)
set.seed(100)

for (i in 1:100) {
  #i=1
  print(i)
  optim_fit <- optim_cop_like(fit1,fit2,copula_model)
  mean_margins <- rbind(mean_margins, c(optim_fit[[1]],optim_fit[[2]]))
}

par(mfrow=c(2,3))
for (i in 1:ncol(mean_margins)) {
  hist(mean_margins[,i],main=colnames(mean_margins)[i])
}

colMeans(mean_margins)
cov_sims=cov(mean_margins)

sqrt(diag(cov_sims))

sqrt(cov_sims[2,2]+cov_sims[5,5]+2*cov_sims[2,5])
sqrt(cov_sims[3,3]+cov_sims[6,6]+2*cov_sims[3,6])*10

summary(fit2)[1]-summary(fit1)[1]
sqrt(var(mean_margins[,1])+var(mean_margins[,2])-2*cov(mean_margins)[1,2])
sqrt(vcov(fit1)[1,1]+vcov(fit2)[1,1]-2*cov(mean_margins)[1,2])
