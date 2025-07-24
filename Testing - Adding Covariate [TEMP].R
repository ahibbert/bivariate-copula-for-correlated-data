# OVERVIEW: Run any single simulation and plot the results of the fitted models against the estimated true values.
# WARNING: You must restart R before each run because loading the GJRM function breaks GAMLSS

######### INPUT PARAMETERS ##########

dist="GA" # Distribution to use for the simulation
sims=100 # Number of simulations to run for the true value calculation
n=1000 # Number of observations in the dataset (note there are two timepoints so this will create n x 2 data points)
a=.5 # Parameter a for the distribution
b=2 # Parameter b for the distribution
c=NA # Parameter c for the distribution (not used for GA)
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



################################## TEST CODE TO ARCHIVE
#dist="NO";a=1;b=1;c=.5;mu1=1;mu2=2;x1=1;x2=1;n=1000
#dist="GA" ;sims=100 ;n=1000 ;a=.5 ;b=2 ;c=NA ;mu1=1 ;mu2=2 ;x1=1 ;x2=1 # Covariate 2 value
#dist="PO"; a=NA;b=1;c=2;mu1=1;mu2=2;x1=1;x2=1;n=1000
dist="LO"; a=NA;b=NA;c=.5;mu1=.25;mu2=.75;x1=1;x2=1;n=1000

dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
dataset_test=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2); dataset_test$patient=dataset_test$patient+1000
true=simCovariateMLEs(sims=sims,n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2,trace=TRUE)

library("gamlss2")
require(gamlss)
require(gee)
require(lme4)
require(MASS)
require(gamlss.mx)
library(mgcv)
library(glmtoolbox)

if(dist=="GA") {
  model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=GA)
  invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=Gamma(link = "log"), maxit=1000)))
  invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
  invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                   , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
  invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset, family=Gamma(link="log"))))
  model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=Gamma(link="log"))
} else if(dist=="NO") {
  invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=gaussian, maxit=1000)))
  invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
  invisible(capture.output(model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NO)))
  invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                   , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
  model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)
  model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=gaussian)
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
  model_re_nosig <- gamlss2(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI)
  #invisible(capture.output(model_re <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
  invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                                                   , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
  
  model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) +as.factor(sex)+age+ (1|patient), data=dataset,family=binomial)
  
  model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=binomial)
}

model_list= list(
  glm=model_glm,
  gee=model_gee,
  re_nosig=model_re_nosig,
  lme4=model_lme4,
  gamm=model_gamm$gam
)

dfun=if(dist=="NO"){dNO}else if(dist=="GA"){dGA}else if(dist=="PO"){dNBI}else if(dist=="LO"){dNBI} else {stop("Unsupported distribution")}

pred_list=list()
evaluation=matrix(NA,ncol=length(model_list),nrow=3,dimnames=list(c("MSEP","LS","logLik"), names(model_list)))
for (model_name in names(model_list)) {
  model=model_list[[model_name]]
  pred_list[[model_name]]= predict(model, newdata=dataset_test, type="response", allow.new.levels=TRUE)
  if(dist=="LO") {
    ls_temp=log(dfun(dataset_test$random_variable, mu=pred_list[[model_name]])) #Evaluated at the TRUE sigma. Does this make sense?
    ls_temp[is.infinite(ls_temp)]=0; # Replace infinite values with 0 until we have a fix
    evaluation["LS",model_name]=sum(ls_temp)
  } else {
    ls_temp=log(dfun(dataset_test$random_variable, mu=pred_list[[model_name]],sigma=exp(sigma_temp))) #Evaluated at the TRUE sigma. Does this make sense?
    ls_temp[is.infinite(ls_temp)]=0; # Replace infinite values with 0 until we have a fix
    evaluation["LS",model_name]=sum(ls_temp)
  }
  evaluation["MSEP",model_name]= mean((dataset_test$random_variable - pred_list[[model_name]])^2)
  sigma_temp=dataset_test$time* true$coefficients["s2"] - (dataset_test$time-1) * true$coefficients["s2"]
  
  if(model_name=="gamm") {
    model_temp=model_gamm$lme
    evaluation["logLik",model_name]=logLik(model_temp)
  } else {
    evaluation["logLik",model_name]=logLik(model)  
  }
}
evaluation

#############Proving we can combine multiple gammas########
shape1=6
shape2=6
scale1=10
scale2=30

gamma1=rgamma(10000,shape=shape1,scale=scale1)
gamma2=rgamma(10000,shape=shape2,scale=scale2)

meanguess=(shape1*scale1 + shape2*scale2)/2
#meanweighted
#varguess=(shape2*scale2*shape2*(scale2^2)+(shape1*scale1*shape1*(scale1^2)))/(shape2*scale2+shape1*scale1)
#varweighted
varguess=(shape2*(scale2^2)*shape2*(scale2^2)+(shape1*(scale1^2)*shape1*(scale1^2)))/(shape1*(scale1^2)+shape2*(scale2^2))

scaleguess=varguess/meanguess
shapeguess=(meanguess^2)/varguess

out=matrix(data=c(
  shape1,scale1,
  shape1*scale1 #shape*scale=mean of 6
  ,shape1*(scale1^2) #variance of 12 3 * 2^2
  ,shape2,scale2
  ,shape2*scale2 #shape*scale=mean of 12
  ,shape2*(scale2^2) #variance of 3 * 4^2 = 48
  ,(mean(c(gamma1,gamma2))^2)/var(c(gamma1,gamma2))
  ,var(c(gamma1,gamma2))/mean(c(gamma1,gamma2))
  ,mean(c(gamma1,gamma2)) #shape*scale=mean of 9 - so we add the scales and divide by 2
  ,var(c(gamma1,gamma2)) #This is the scale parameter (4.3????)
  ,shapeguess
  ,scaleguess
  ,meanguess
  ,varguess
  ,(((mean(c(gamma1,gamma2))^2)/var(c(gamma1,gamma2)))/shapeguess) -1
  ,((var(c(gamma1,gamma2))/mean(c(gamma1,gamma2)))/scaleguess) -1
  ,(mean(c(gamma1,gamma2))/meanguess) -1 #shape*scale=mean of 9 - so we add the scales and divide by 2
  ,(var(c(gamma1,gamma2))/varguess) -1#This is the scale parameter (4.3????)
  
),ncol=4,byrow=TRUE,dimnames=list(c(paste("gamma 1:"),paste("gamma 2:"),paste("gamma (1,2) calc"),"Par method","Diff"),c("shape","scale","mean","var")))

round(out,3)

#fits=fitDist(gamma1)
#fits$fits
#fits=fitDist(c(gamma1,gamma2))
#fits$fits
###VARIANCE GUESS IS CORRECT OK SO WE CAN ESTIMATE SCALE


##########OK What about negbin 
