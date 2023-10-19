#Importing required packages
require(gamlss); require(gee); require(lme4);
require(mgcv); require(geepack)

#Setting parameters for extreme case example
set.seed(1)
a=0.25; b=1.75; mu1=1; mu2=2; n=100

#Simulating bivariate gamma of Nadarajah and Gupta
w<-rbeta(n,a,b)
gamma_c_mu1<-w * rgamma(n,shape=a+b,scale=1/mu1)
gamma_c_mu2<-w * rgamma(n,shape=a+b,scale=1/mu2)

#Setting up data in correct format for random effect model
patient<-as.factor(seq(1:n))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0)
                             ,cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")

#Running GLM, GEE and multiple GLMM packages
model_glm <- glm(random_variable~as.factor(time==1),data=dataset
                 ,family=Gamma(link = "log"),maxit=10000)
model_gee<-geese(random_variable~as.factor(time==1), id=patient
                 , data=dataset, family=Gamma(link="log")
                 , mean.link = "log", corstr = "exchangeable"
                 , control=geese.control(trace=TRUE, maxit=10000))
model_lme4<-glmer(random_variable~as.factor(time==1) + (1 | patient)
                  , data=dataset, family=Gamma(link="log"))
model_gamm<-gamm(random_variable~as.factor(time==1)
                 , random = list(patient=~1)
                 ,data=dataset, family=Gamma(link="log")
                 ,niterPQL=1000)
model_re_nosig <- gamlss(random_variable~as.factor(time==1)
                         + random(as.factor(patient))
                         , data=dataset, family=GA(), method=RS())
model_re <- gamlss(formula=random_variable~as.factor(time==1)
                   + random(as.factor(patient))
                   , sigma.formula=~as.factor(time==1)
                   , data=dataset, family=GA() , method=RS())

#Running GJRM - note GJRM package is loaded after gamlss
require(GJRM)
eq.mu.1 <- gamma_c_mu1~1
eq.mu.2 <- gamma_c_mu2~1
fl <- list(eq.mu.1, eq.mu.2)
model_copula<-gjrm(fl, margins = c("GA" , "GA"), copula = "C0"
                   , data=data.frame(gamma_c_mu1,gamma_c_mu2), model="B")
model_copula_n<-gjrm(fl, margins = c("GA" , "GA"), copula = "N"
                     , data=data.frame(gamma_c_mu1,gamma_c_mu2), model="B")
