theta=4
a_1=rgamma(10000,shape=.25*10,scale=10)
a_2=rgamma(10000,shape=.25*12,scale=12)
re<-rnorm(10000,mean=0,sd=theta)

par(mfrow=c(3,3))
#hist(a_1);hist(a_2);hist(re)
hist(log(a_1));hist(log(a_2));hist(re);
hist((a_1));hist((a_2));hist(exp(re));

#fits<-fitDist(exp(re))
#fits$fits
histDist(exp(re),family="LOGNO",nbins = 100)

#SO the marginal distribution of the model is actually....?
a_1_cond_re<-exp(log(a_1)-(re))
a_2_cond_re<-exp(log(a_2)-(re))

fits<-fitDist(a_1_cond_re)
fits$fits
par(mfrow=c(3,3))
histDist(a_1_cond_re,family="GG",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCCG",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCCGo",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCPE",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCPEo",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCTo",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="BCT",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="GB2",nbins=200,ylim = c(0,.3),xlim=c(0,50))
histDist(a_1_cond_re,family="GA",nbins=200,ylim = c(0,.3),xlim=c(0,50))

fits2<-fitDist(a_2_cond_re)
fits2$fits
histDist(a_2_cond_re,family="BCTo",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="BCT",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="GB2",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="BCPEo",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="BCPE",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="GG",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="BCCG",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="BCCGo",nbins=100,ylim = c(0,.5))
histDist(a_2_cond_re,family="GA",nbins=100,ylim = c(0,.5))



hist(exp(norm_data+log(a_1)))



gamma_fit_1 <- gamlss(exp(norm_data+log(a_1))~1,family=GA())
gamma_fit_2 <- gamlss(exp(norm_data+log(a_2))~1,family=GA())

margin_1_unif<-pnorm(gamma_fit_1$residuals)
margin_2_unif<-pnorm(gamma_fit_2$residuals)

hist(margin_1_unif)
hist(margin_2_unif)

plot(margin_1_unif,margin_2_unif)




##############Random effect
m1=2;s1=2
m2=4;s2=4
theta=2

a<-rnorm(100000,mean=m1,sd=sqrt(s1^2-theta^2))
b<-rnorm(100000,mean=m2,sd=sqrt(s2^2-theta^2))
re<-rnorm(100000,mean=0,sd=theta)

c(mean(a),sd(a),m1,sqrt(s1^2-theta^2))
c(mean(b),sd(b),m2,sqrt(s2^2-theta^2))

c<-a+re
d<-b+re

c(mean(c),sd(c),m1,sqrt(s1^2))
c(mean(d),sd(d),m2,sqrt(s2^2))

#Cov
c(theta^2, cov(c,d))

#Cor
#c((theta^2)/((sqrt(s1^2+theta^2))*sqrt(s2^2+theta^2)),cor(c,d))
c((theta^2)/(s1*s2),cor(c,d))




#######Normal with clayton copula
library(copula)

tau <- 0.75
theta <- iTau(claytonCopula(), tau = tau)
d <- 2
cc <- claytonCopula(theta, dim = d)
n <- 1000
set.seed(271)

clayton_copula<-rCopula(1000,claytonCopula(param=iTau(claytonCopula(), tau = .5),dim=2))
r_1<-qnorm(clayton_copula[,1],mean=1,sd=1)
r_2<-qnorm(clayton_copula[,2],mean=2,sd=1)

dataset<-rbind(cbind(1:length(r_1),rep(1,length(r_1)),r_1),cbind(1:length(r_2),rep(2,length(r_2)),r_2))
colnames(dataset)<-c("patient","time","random_variable")
dataset<-dataset[order(dataset[,"patient"]),]

library(gamlss)
library(lme4)
gamlss(random_variable~-1+as.factor(time==2),data=as.data.frame(dataset))
model_gamlss4 <- gamlss(random_variable~-1+as.factor(time==2)+random(as.factor(patient)),data=as.data.frame(dataset))
summary(model_gamlss4)
model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==2) + (1|patient), data=as.data.frame(dataset))
summary(model_lme4)       


#######Normal with log link... is OK?
library(copula)

clayton_copula<-rCopula(1000,normalCopula(.75))
r_1<-qnorm(clayton_copula[,1],mean=1,sd=1)
r_2<-qnorm(clayton_copula[,2],mean=2,sd=1)

dataset<-rbind(cbind(1:length(r_1),rep(1,length(r_1)),r_1),cbind(1:length(r_2),rep(2,length(r_2)),r_2))
colnames(dataset)<-c("patient","time","random_variable")
dataset<-dataset[order(dataset[,"patient"]),]

library(gamlss)
library(lme4)
gamlss(random_variable+100~-1+as.factor(time==2),data=as.data.frame(dataset),family=NO(mu.link="log"))
model_gamlss4 <- gamlss(random_variable+100~-1+as.factor(time==2)+random(as.factor(patient)),data=as.data.frame(dataset),family=NO(mu.link="log"))
summary(model_gamlss4)
#model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==2) + (1|patient), data=as.data.frame(dataset))
#summary(model_lme4)       









