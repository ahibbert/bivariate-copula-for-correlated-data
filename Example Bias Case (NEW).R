############# Overview of file
#############
#############

########### 0. DATA SETUP ##############

############### 0.2 Applications #############
####CHOOSE A DATASET

#Lipids Data
#require(sas7bdat)
#lipid <- read.sas7bdat("LipidsData.sas7bdat")
#lipids_merged<-(merge(lipid[lipid$MONTH==0,],lipid[lipid$MONTH==24,],by="PATIENT"))
lipids_merged<-readRDS("lipids_merged.rds")
gamma_c_mu1<-lipids_merged$TRG.x
gamma_c_mu2<-lipids_merged$TRG.y

lipids_merged<-readRDS("lipids_merged.rds")
gamma_c_mu1<-lipids_merged$HDL.x
gamma_c_mu2<-lipids_merged$HDL.y

#Stock prices over 10 years

#ASX2018<-read.table("20180102.txt", header=FALSE, sep=",")
#ASX1998<-read.table("19980102.txt", header=FALSE, sep=",")
#ASX98_18<-merge(ASX1998,ASX2018,by="V1")
ASX98_18<-readRDS("ASX98_18.rds")
gamma_c_mu1<-ASX98_18$V6.x
gamma_c_mu2<-ASX98_18$V6.y

#Avocado prices
#avo<-read.table("avocado prices.csv", header=T, sep=",")
avo<-readRDS("avo.rds")
gamma_c_mu1<-avo[avo$Date=="4/01/2015","AveragePrice"]
gamma_c_mu2<-avo[avo$Date=="25/03/2018","AveragePrice"]

# c.Setting up as longitiduinal structured data
patient<-as.factor(seq(1:length(gamma_c_mu1)))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")
dataset<-dataset[order(dataset$patient),]
dist="GA"


############## 0.1 Simulation ##################
source("common_functions.R")

# a. Simulation parameters
set.seed(1000);options(scipen=999);
#dist="NO";a=1; b=2; c=0.75; mu1=1; mu2=2; n=1000
dist="GA";a=.25; b=1.75; c=NA; mu1=10; mu2=12; n=1000
#dist="GA";a=.2; b=.2; c=NA; mu1=10; mu2=12; n=1000
#dist="PO";a=NA; b=1; c=.1; mu1=5; mu2=5; n=1000 ## Highly skewed
#dist="PO";a=NA; b=.5; c=9; mu1=5; mu2=5; n=1000 ## Not highly skewed

dataset <- generateBivDist(n,a,b,c,mu1,mu2,dist)

plotDist(dataset,dist)

results<-fitBivModels(data=dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE)
if(dist=="NO"){clean_results<-results}else{clean_results<-cbind(results,round(exp(results[,c(1,2)]),4));
clean_results<-cbind(clean_results[,c(9,10)],clean_results[,1:8]);colnames(clean_results)<-c("mu_1","mu_2",colnames(clean_results[,c(3:10)]))}
clean_results

###########Fitting individual models##############

library(gamlss)
library(lme4)
library(gee)
library(gamlss.mx)

model_glm <- gamlss(formula=random_variable~as.factor(time==1), data=dataset, family=GA()) 
model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")

model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")

model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                         , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)
summary(model_re_np)
summary(model_glm)
summary(model_gee)
summary(model_lme4)
summary(model_re_nosig)



#########Estimating marginal means from predictions ################



model_glm <- gamlss(formula=random_variable~-1+as.factor(time==1), data=dataset, family=GA()) 
model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")
model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
model_re_nosig_ident <- gamlss(formula=1+random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA(mu.link="identity")) #UNBIASED
model_re_no <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO(mu.link = "identity",sigma.link="log")) #UNBIASED
model_re_no_log <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO(mu.link = "log",sigma.link="log")) #BIASED but better AIC

model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))

mean(predict(model_re_nosig)[dataset$time==0]+model_re_nosig$mu.coefSmo[[1]]$coef)
mean(predict(model_re_nosig)[dataset$time==1]+model_re_nosig$mu.coefSmo[[1]]$coef)



mean(predict(model_re_nosig,type="response")[dataset$time==0])
mean(predict(model_re_nosig,type="response")[dataset$time==1])

hist(predict(model_re_nosig)[dataset$time==0])
hist(model_re_nosig$mu.coefSmo[[1]]$coef)

hist(predict(model_re_nosig)[dataset$time==0]+model_re_nosig$mu.coefSmo[[1]]$coef)


hist(predict(model_re_nosig)[dataset$time==0])
hist(model_re_nosig$mu.coefSmo[[1]]$coef)

model_re_nosig$mu.coefSmo[[1]]$coef

#install.packages("emmeans")
#library(emmeans)
#EMM<-emmeans(model_lme4,specs = all.vars(terms(model_lme4)[[3]]) )
#EMM
#pairs(EMM)

install.packages("ggeffects")
library(ggeffects)
predict_response(model_lme4)

#install.packages("marginaleffects")
library(marginaleffects)
summary(model_lme4)
glm_mod<-model_lme4


predictions(
  glm_mod,
  by = "time",type = "response")
avg_predictions(
  glm_mod,
  by = "time",type = "response")
predictions(
  glm_mod,
  by = "time",type = "link")
avg_predictions(
  glm_mod,
  by = "time",type = "link")
predictions(
  glm_mod,
  by = "time",
  newdata = datagrid(grid_type = "balanced"), type="response")
predictions(
  glm_mod,
  by="time",
  newdata = datagrid(model=glm_mod, time=0:1, grid_type = "counterfactual"), type="response")


re_gamlss<-model_re_nosig$mu.coefSmo[[1]]$coef
re_lme4<-getME(model_lme4,"b")[,1]

pred_lme4<- predict(model_lme4,type="link", se.fit=TRUE)
pred_gamlss<- predict(model_re_nosig,type="link", se.fit=TRUE)

#######CONDITIONAL MODEL RESULTS

exp(mean((pred_lme4$fit)[dataset$time==0])) #######This is the result from the conditional model
exp(mean((pred_lme4$fit)[dataset$time==1]))
exp(mean((pred_gamlss$fit)[dataset$time==0]))
exp(mean((pred_gamlss$fit)[dataset$time==1]))

(mean((pred_lme4$fit)[dataset$time==0])) #######This is the result from the conditional model
(mean((pred_lme4$fit)[dataset$time==1]))
(mean((pred_gamlss$fit)[dataset$time==0]))
(mean((pred_gamlss$fit)[dataset$time==1]))

mean(exp(pred_lme4$fit)[dataset$time==0]) ##########But this is the true marginal result
mean(exp(pred_lme4$fit)[dataset$time==1])
mean(exp(pred_gamlss$fit)[dataset$time==0])
mean(exp(pred_gamlss$fit)[dataset$time==1])


pred <- predictions(
  glm_mod,
  newdata = datagrid(patient = NA,
                     time = 0:1),
  re.form = NA)

#glm_mod <- glm(vs ~ hp + am, data = mtcars, family = binomial)

avg_predictions(glm_mod)$estimate
#> [1] 0.06308965

## Step 1: predict on the link scale
p <- predictions(glm_mod, type = "link")$estimate

predictions(glm_mod, type = "response",by="time")
## Step 2: average
mean_p <- mean(p)
## Step 3: backtransform
glm_mod@resp$family$linkinv(mean_p)
#> [1] 0.06308965


#mod<-model_re_nosig
pred <- predictions(mod, type = "response")
#Equivalent to
mean(predict(model_re_nosig,type="response")[dataset$time==0])
mean(predict(model_re_nosig,type="response")[dataset$time==1])

mean(pred$random_variable[pred$time==0])
mean(pred$random_variable[pred$time==1])

avg_predictions(mod,by="time") #unbiased
avg_predictions(mod,by="time",type="link") #Biased
avg_predictions(mod, by = "time",
                vcov = "HC3",
                conf_level = .9)


predictions(
  mod,
  by="time",type="link")


predictions(
  mod,
  type = "link",
  by = "time",
  newdata = datagrid(time = 0:1, grid_type = "counterfactual"))


plot(model_lme4)

##########Looking at the random effects shape ##########
re_gamlss<-model_re_nosig$mu.coefSmo[[1]]$coef
re_lme4<-getME(model_lme4,"b")[,1]

re_gamlss_logno<-model_re_nosig_logno$mu.coefSmo[[1]]$coef

par(mfrow=c(2,2))
hist(re_gamlss,main="GAMLSS Random Effect")
hist(re_lme4,main="LME4 Random Effect")
hist(exp(re_gamlss),main="EXP GAMLSS Random Effect")
hist(exp(re_lme4),main="EXP LME4 Random Effect")

fitDist(re_gamlss,type="realline")$fits ###EGB2
fitDist(re_lme4,type="realline")$fits

par(mfrow=c(3,2))
#histDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(re_gamlss)),nbins=100,family="GA")
histDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(re_lme4)),nbins=20,family="GA",main="GA: Time 1")
#histDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(re_gamlss)),nbins=100,family="GA")
histDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(re_lme4)),nbins=20,family="GA",main="GA: Time 2")
histDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(re_lme4)),nbins=20,family="GG",main="GG: Time 1")
histDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(re_lme4)),nbins=20,family="GG",main="GG: Time 2")
histDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(re_lme4)),nbins=20,family="GB2",main="GB2: Time 1")
histDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(re_lme4)),nbins=20,family="GB2",main="GB2: Time 2")

fit1<-fitDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(model_re_nosig$mu.coefSmo[[1]]$coef)))
fit2<-fitDist(exp(log(dataset[dataset[,"time"]==0,"random_variable"])-(getME(model_lme4,"b")[,1])))
fit3<-fitDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(model_re_nosig$mu.coefSmo[[1]]$coef)))
fit4<-fitDist(exp(log(dataset[dataset[,"time"]==1,"random_variable"])-(getME(model_lme4,"b")[,1])))

fit1$fits
fit2$fits
fit3$fits
fit4$fits

model_glm <- glm(random_variable~-1+as.factor(time==1), data=dataset, family=Gamma(link = "log"), maxit=1000)
model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")
model_re_nosig <- gamlss(formula=100+random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA(mu.link="identity")) 


model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=BCPE()) 
summary(model_re_nosig)



model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                 , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)

model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))


model_glm<-glm(random_variable~-1+as.factor(time==1), data=dataset, family=gaussian, maxit=1000)

model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset)
summary(model_lme4)

as.data.frame(VarCorr(model_lme4))

(dotplot(ranef(model_lme4, postVar=TRUE)))

fixef(model_lme4)
coef(summary(model_lme4))
AIC(model_lme4)
AIC(model_glm)




mean(dataset[dataset[,"time"]==1,"random_variable"])/mean(dataset[dataset[,"time"]==0,"random_variable"])
mean(dataset[dataset[,"time"]==1,"random_variable"]/dataset[dataset[,"time"]==0,"random_variable"])
hist(dataset[dataset[,"time"]==1,"random_variable"]/dataset[dataset[,"time"]==0,"random_variable"],breaks=20)

library(copula)
exchTest(cbind(dataset[dataset[,"time"]==0,"random_variable"],dataset[dataset[,"time"]==1,"random_variable"]))

#############Investigating gamlss (5) fit errors

require(gamlss.mx)
model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                         , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)

model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~1, data=dataset, family=NBI())

model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NBI(), method=CG(1000))

summary(model_re)
summary(model_re_nosig)
summary(model_re_np)

library(gamlss.mx)

model_re1 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
         , g.control = gamlss.control(trace = TRUE), mixture="np")
model_re2 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                     , g.control = gamlss.control(trace = TRUE), mixture="gq")
model_re3 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                     , g.control = gamlss.control(trace = TRUE), mixture="np",K=20)
model_re4 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                     , g.control = gamlss.control(trace = TRUE), mixture="gq",K=20)
model_re12 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                      , g.control = gamlss.control(trace = TRUE), mixture="np",K=2)
model_re22 <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                      , g.control = gamlss.control(trace = TRUE), mixture="gq",K=2)

summary(model_re1)
summary(model_re2)
summary(model_re3)
summary(model_re4)
summary(model_re12)
summary(model_re22)

summary(model_re22)[4]

################Investigating copula fit


require(GJRM)
eq.mu.1 <- formula(random_variable~1)
eq.mu.2 <- formula(random_variable.1~1)
fl <- list(eq.mu.1, eq.mu.2)

gamma_c_mu1<-dataset[dataset$time==0,]
gamma_c_mu2<-dataset[dataset$time==1,]

if(dist=="NO"){margin_dist="N"}
if(dist=="GA"){margin_dist="GA"}

model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
summary(model_copula_n)












################# 1. Investigating the dependence structure and best copula fit
################## 1.1 bias case plotting #######
require(GJRM)
library(MASS)
library(psych)
library(copula)
library(VineCopula)
library(gamlss)
library(moments)
require(ggpubr)
require(ggplot2)
require(dglm)
library(latex2exp)

####UKNOWN FIT u,v

gamma_c_mu1<-dataset[dataset$time==0,"random_variable"]
gamma_c_mu2<-dataset[dataset$time==1,"random_variable"]

fit1<-fitdistr(gamma_c_mu1,"gamma")
fit2<-fitdistr(gamma_c_mu2,"gamma")
u<-pgamma(gamma_c_mu1,shape=fit1$estimate[1],rate=fit1$estimate[2])
v<-pgamma(gamma_c_mu1,shape=fit2$estimate[1],rate=fit2$estimate[2])

fittedClayton=rCopula(n, archmCopula(family="clayton",BiCopSelect(u,v,family=3)$par,dim=2))
fittedTDist=rCopula(n, tCopula(BiCopSelect(u,v,family=2)$par,dim=2,df=BiCopSelect(u,v,family=2)$par2))

# Plot density as points
z<-ggplot(data = as.data.frame(gamma_c_mu1)) +
  geom_histogram(data = as.data.frame(gamma_c_mu1), aes(x=gamma_c_mu1, y=..density..),bins=70) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu1, y=dgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])), color="blue", size = .75) +
  ylim(0,.6) +
  xlim(0,40) +
  labs(x=TeX("$Y_1$")) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
x<-ggplot(data = as.data.frame(gamma_c_mu2)) +
  geom_histogram(data = as.data.frame(gamma_c_mu2), aes(x=gamma_c_mu2, y=..density..),bins=70) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu2, y=dgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])), color="blue", size = .75) +
  ylim(0,.6) +
  labs(x=TeX("$Y_2$")) +
  xlim(0,40) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
c<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
  geom_point(size=0.5,color="black") + 
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$")) +
  xlim(0,40) +
  ylim(0,40) +
  geom_smooth(method="loess", level=.99) 
d<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=15,color="black") + 
  scale_fill_brewer() +
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$"),fill="density")+
  xlim(0,1) +
  ylim(0,1) 
e<-ggplot(data=as.data.frame(fittedClayton),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=15,color="black") +
  scale_fill_brewer() +
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$"))+
  xlim(0,1) +
  ylim(0,1) 
f<-ggplot(data=as.data.frame(fittedTDist),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d( contour_var="density",bins=15,color="black") + 
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$"))+
  xlim(0,1) +
  ylim(0,1) 
ggarrange(z,x,c,nrow=1)
#ggsave(file="example_bias_case_margin_plots.jpeg",last_plot(),width=14,height=4,dpi=300)
ggarrange(d,e,f,nrow=1)
#ggsave(file="example_bias_case_contour_plots.jpeg",last_plot(),width=14,height=4,dpi=300)

#plot(u,v,main="Uniform transform of both marginals",xlab="Time 1 Marginal Gamma (Uniform Transform)",ylab="Time 2 Marginal Gamma (Uniform Transform)")
#plot(fittedClayton,main="Simulation of Fitted Clayton Copula",xlab="Fitted Time 1 Marginal Gamma",ylab="Fitted Time 2 Marginal Gamma")
#plot(fittedTDist,main="Simulation of Fitted Normal Copula",xlab="Fitted Time 1 Marginal Gamma",ylab="Fitted Time 2 Marginal Gamma")

#par(mfrow=c(1,3))
#persp(kde2d(u,v,h=.4,n=65),main="Uniform transform of marginals",zlim=c(0,4))
#persp(kde2d(fittedClayton[,1],fittedClayton[,2],h=.4,n=65),main="Simulated fitted copula",zlim=c(0,4))
#persp(kde2d(fittedTDist[,1],fittedTDist[,2],h=.4,n=65),main="Simulated fitted normal",zlim=c(0,4))

################ 1.2 Applications case plotting #############

library(moments)
skewness(gamma_c_mu1)
skewness(gamma_c_mu2)

require(MASS)  
require(copula)
require(gamlss)
require(dglm)
require(ggplot2)
require(ggpubr)
library(latex2exp)
u<-0
v<-0

fit <- dglm(gamma_c_mu1~1, family=Gamma(link="log"))
mu=exp(fit$coefficients)
shape=exp(-1*fit$dispersion.fit$coefficients)
scale <- mu/shape

u<-pgamma(gamma_c_mu1,shape=shape,scale=scale)

fit <- dglm(gamma_c_mu2~1, family=Gamma(link="log"))
mu=exp(fit$coefficients)
shape=exp(-1*fit$dispersion.fit$coefficients)
scale <- mu/shape

v<-pgamma(gamma_c_mu2,shape=shape,scale=scale)

#ASX# bins=30; binsc=15; limx=50; limy=0.15
#lipids 
bins=20; binsc=15; limx=3; limy=2
# Plot density as points
a<-ggplot(data = as.data.frame(gamma_c_mu1)) +
  geom_histogram(data = as.data.frame(gamma_c_mu1), aes(x=gamma_c_mu1, y=..density..),bins=bins) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu1, y=dgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])), color="blue", size = .75) +
  ylim(0,limy) +
  xlim(0,limx) +
  labs(x=TeX("$Y_1$")) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
b<-ggplot(data = as.data.frame(gamma_c_mu2)) +
  geom_histogram(data = as.data.frame(gamma_c_mu2), aes(x=gamma_c_mu2, y=..density..),bins=bins) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu2, y=dgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])), color="blue", size = .75) +
  ylim(0,limy) +
  labs(x=TeX("$Y_2$")) +
  xlim(0,limx) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
c<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
  geom_point(size=0.5,color="black") + 
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$")) +
  geom_smooth(method="loess", level=.95) +
  xlim(0,limx) +
  ylim(0,limx) 
d<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=binsc,color="black") + 
  scale_fill_brewer() +
  labs(x = TeX("$Y_1$"), y=TeX("$Y_2$"),fill="density")
#e<-ggplot(data=as.data.frame(fittedClayton),aes(x=V1,y=V2)) + 
#  #geom_point(size=0.25,color="black") + 
#  geom_density_2d(contour_var="density",bins=20,color="black") +
#  scale_fill_brewer() +
#  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")
#f<-ggplot(data=as.data.frame(fittedTDist),aes(x=V1,y=V2)) + 
#  #geom_point(size=0.25,color="black") + 
#  geom_density_2d( contour_var="density",bins=20,color="black") + 
#  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted normal copula")
ggarrange(a,b,d,nrow=1)
#ggsave(file="applications_asx.jpeg",last_plot(),width=14,height=4,dpi=300)
#ggarrange(d,e,f,common.legend = TRUE,nrow=1,legend="right")
##ggsave(file="example_bias_case_contour_plots.jpeg",last_plot(),width=14,height=4,dpi=300)

################# Getting standard errors for bias case  #########

getSE(10,12,.25,1.75,1000)


getSE <- function(origmu1,origmu2,a,b,n) {

  results <- c(0,0,0,0,0)
  #B_2 parameterisation
  mu1=log(a*1/origmu1); mu2=log(a*1/origmu2);  ###Changes a lot based on y1/y2
  par=c(mu1,mu2,a,b,exp(mu1)/a,exp(mu2)/a)
  results[1:4]=sqrt(numericalDerivativeSE(par,parameterisation="B2")/n) #sqrt(numericalDerivativeSE(par))/sqrt(n)
  
  #B_t parameterisation
  mu1=log(a*1/origmu1); mu2=log((1/origmu2)/(1/origmu1));  ###Changes a lot based on y1/y2
  par=c(mu1,mu2,a,b,exp(mu1)/a,exp(mu1+mu2)/a)
  results[5]=sqrt(numericalDerivativeSE(par,parameterisation="Bt")/n)[2] #sqrt(numericalDerivativeSE(par))/sqrt(n)
  
  return(results)

}


rbind(summary_glm, summary_gee,summary_re_nosig,summary_re,summary_cop,summary_cop_n,actuals)

########## 5. Investigating GLMM fits #######

library(gamlss)
library(mgcv)
set.seed(1000)

####Benchmark Gamma GLM
gamlss_fit_0<-gamlss(formula=random_variable~as.factor(time==1)
       , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
       , method=CG(10000))
plot(gamlss_fit_0)

####Gamma + Random Effect
gamlss_fit_GA<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=GA()
                      #, method=CG(10000)
                      )
plot(gamlss_fit_GA)
####Gamma + Random Effect
#gamlss_fit_GA_s<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
#                      , sigma.formula=~as.factor(time==1)
#                      , data=dataset, family=GA(),method=CG(1000))
#plot(gamlss_fit_GA)

#####Look at the fitted random effect
#So - 1. the arndom effect doesn't look normal, 2. 
par(mfrow=c(2,2))
hist(log(gamma_c_mu1),xlim=c(-8,0),main="Margin 1 (log)",breaks =20)
hist(log(gamma_c_mu2),xlim=c(-8,0),main="Margin 2 (log)",breaks =20)
hist(gamlss_fit_GA$mu.coefSmo[[1]]$coef,xlim=c(-4,4),main="Random effect",breaks =20)
hist(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef,xlim=c(-8,0),main="Margin 2 minus Random effect",breaks =20)

####Looking at the 
plot.new()
par(mfrow=c(2,2))
histDist(gamma_c_mu1,main = "Original margin 1 (gamma)",family=GA(),ylim=c(0,14),nbins=100,xlim=c(0,.5))
histDist(gamma_c_mu2,main = "Original margin 2 (gamma)",family=GA(),ylim=c(0,14),nbins=50,xlim=c(0,.5))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),ylim=c(0,14),nbins=400,xlim=c(0,.5),main="Margin 1 minus random effect")
histDist(exp(log(gamma_c_mu2)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),ylim=c(0,14),nbins=600,xlim=c(0,.5),main="Margin 2 minus random effect")
###So adding RE makes it MORE skewed

#Not a substantial improvement
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),nbins=300,xlim=c(0,.5),main="GA",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GB2(),nbins=300,xlim=c(0,.5),main="GB2",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=BCTo(),nbins=300,xlim=c(0,.5),main="BCTo",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GP(),nbins=300,xlim=c(0,.5),main="GP",ylim=c(0,14))

histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),nbins=300,xlim=c(0,.5),main="GA",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GB2(),nbins=300,xlim=c(0,.5),main="LOGNO",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=BCTo(),nbins=300,xlim=c(0,.5),main="LOGNO2",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GP(),nbins=300,xlim=c(0,.5),main="GP",ylim=c(0,14))


fitDist(exp(log(gamma_c_mu2)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),type="realplus")$fits
#GB2      BCTo        GP   PARETO2  PARETO2o     BCPEo     BCCGo        GG    LOGNO2     LOGNO      WEI3      WEI2       WEI 
#-2555.342 -2554.560 -2553.936 -2553.936 -2553.936 -2553.029 -2550.251 -2549.281 -2528.517 -2528.517 -2457.936 -2457.936 -2457.936 
#GIG        GA       EXP        IG    IGAMMA 
#-2439.929 -2382.074 -2285.360 -2054.721 -1926.681 

fitDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),type="realplus")$fits
#BCTo      BCPEo    PARETO2         GP   PARETO2o        GB2      BCCGo         GG     LOGNO2      LOGNO       WEI2       WEI3 
#-2244.0461 -2244.0391 -2243.6023 -2243.6023 -2243.6023 -2243.0053 -2238.9137 -2237.5020 -2215.2062 -2215.2062 -2150.5406 -2150.5406 
#WEI        GIG         GA        EXP     IGAMMA         IG 
#-2150.5406 -2076.3217 -2067.8069 -1910.3904  -967.5800  -790.0647 


#Won't coverge
#gamlss_fit_GB2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
#                      , data=dataset, family=GB2()
#                      , method=CG(1000)
#                      )
#
plot(gamlss_fit_GB2)

#Plots look ridiculous but global deviance is better 
gamlss_fit_BCTo<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=BCTo()
                      , method=RS(100)
                      )

plot(gamlss_fit_BCTo)

##This is the only 2-parameter distribution with a better AIC fit
gamlss_fit_GP<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                        , data=dataset, family=GP()
                        #,method=RS(1000)
                        , method = CG(10000)
                        #, method=mixed(10,1000)
)

plot(gamlss_fit_GP) ###Fit looks ridiculous

gamlss_fit_PARETO2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=PARETO2()
                      #,method=RS(1000)
                      #, method = CG(10000)
                      #, method=mixed(10,1000)
)

plot(gamlss_fit_PARETO2) ###Fit looks ridiculous

gamlss_fit_LOGNO2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                           , data=dataset, family=LOGNO2()
                           #,method=RS(1000)
                           #, method = CG(10000)
                           #, method=mixed(10,1000)
)

plot(gamlss_fit_LOGNO2) ###Fit looks ridiculous

gamlss_fit_LOGNO<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                          , data=dataset, family=LOGNO()
                          #,method=RS(1000)
                          #, method = CG(10000)
                          #, method=mixed(10,1000)
)

plot(gamlss_fit_LOGNO) ###Fit looks ridiculous

gamlss_fit_WEI<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                         , data=dataset, family=WEI2()
                         #,method=RS(1000)
                         #, method = CG(10000)
                         #, method=mixed(10,1000)
)

plot(gamlss_fit_WEI) ###Fit looks ridiculous
summary(gamlss_fit_WEI)


gamlss_fit_GA$aic
gamlss_fit_BCTo$aic
gamlss_fit_GP$aic
gamlss_fit_PARETO2$aic
gamlss_fit_LOGNO2$aic
gamlss_fit_LOGNO$aic
gamlss_fit_WEI$aic


