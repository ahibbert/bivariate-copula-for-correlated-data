
generateBivDist_withCov <- function(n,a,b,c,mu1,mu2,dist,x1,x2,x3) {
  
  if(dist=="GA") {
    #Simulating bivariate random variable according to functional input
    w<-rbeta(n,a,b)
    margin_1<-w*rgamma(n,shape=a+b,scale=mu1)
    margin_2<-w*rgamma(n,shape=a+b,scale=mu2)
    
  }
  if(dist=="NO") {
    require(MASS)
    normData<-mvrnorm(n,mu=c(mu1,mu2),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    margin_1<-normData[,1]
    margin_2<-normData[,2]
    
    sex <- sample(0:1, length(margin_1), replace=TRUE)
    age <- runif(length(margin_1), min=0, max=100)
    trt <- sample(0:1, length(margin_1), replace=TRUE)
    
    time_1=margin_1 + x1*sex + x2*((age)) 
    time_2=margin_2 + x1*sex + x2*((age)) + x3*trt
    
  }
  
  if(dist=="LO") {
    
    require(MASS)
    a=1;b=1
    normData<-mvrnorm(n,mu=c(0,0),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    margin_1<-as.numeric(pnorm(normData[,1])<=mu1)
    margin_2<-as.numeric(pnorm(normData[,2])<=mu2)
    
  }
  
  if(dist=="PO") {
    
    #Compound multiple poisson of Stein & Juritz, 1987
    mixing_dist<-rgamma(n,shape=c,scale=b)
    
    margin_1=vector(length = n) 
    margin_2=vector(length = n) 
    for (i in 1:n) {
      margin_1[i]=rpois(1,mu1*mixing_dist[i])
    }
    for (i in 1:n) {
      margin_2[i]=rpois(1,mu2*mixing_dist[i])
    }
  }
  
  #Transforming data to format required for random effect models
  patient<-as.factor(seq(1:n))
  dataset<-as.data.frame(rbind(cbind(patient,time_1,0,sex,age,trt)
                               ,cbind(patient,time_2,1,sex,age,trt)))
  colnames(dataset)<-c("patient","random_variable","time","sex","age","trt")
  
  dataset<-dataset[order(dataset$patient),]
  
  return(dataset)
}

base_data=generateBivDist_withCov(n=1000,a=1,b=1,c=1,mu1=1,mu2=2,dist="NO",x1=0.5,x2=.01,x3=1)
library(gamlss)
model_gamlss=gamlss(formula=random_variable~as.factor(time==1)+as.factor(sex)+age+as.factor(time==1)*trt,data=base_data)
summary(model_gamlss)
plot(model_gamlss)
term.plot(model_gamlss)


#Next lets try all the fits for normal including GJRM, then automate this






data_1=cbind(base_data[,c("sex","age","trt","time_1")],1)
data_2=cbind(base_data[,c("sex","age","trt","time_2")],2)
colnames(data_1)=colnames(data_2)=c("sex","age","trt","response","time")
          
rbind(data_1,data_2) -> dataset

library(gamlss)
fit=gamlss(response~time+sex+poly(age,2)+as.factor(trt*time),data=dataset)
summary(fit)
plot(fit)


glm(base_data$time_2~base_data$time_1 + base_data$sex + base_data$age + base_data$trt)


invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=gaussian, maxit=1000)))
invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO())))
#invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                 , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))

model_lme4 <- lmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)

model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=gaussian)