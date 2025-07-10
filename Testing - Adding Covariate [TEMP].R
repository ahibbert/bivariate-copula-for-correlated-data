## This is a test script to run the models with covariates

source("common_functions.R")
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="NO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=.3,mu2=.7,dist="LO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="GA",x1=1,x2=1)
dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="PO",x1=1,x2=1)

results=fitBivModels_Bt_withCov(dataset,dist="PO",include="ALL",a=1,b=1,c=.5,mu1=1,mu2=1,calc_actuals=FALSE)

library("glmtoolbox")
invisible(capture.output(model_gee<-overglm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family="nb1(log)", maxiter=25, corstr = "exchangeable")))

#invisible(capture.output(model_glm <- glm.nb(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, maxit=1000)))
#model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, init.beta=model_glm$coefficients,
#                family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")

invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))

#invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
#invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
df=summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),]
df

model_gee

model_gee$logLik

###TEMP PLOTTING

df=results$coefficients
plot.new(); par(mfrow=c(ncol(df),2),las=3)
true_vals = c(1,1,1,1)
factor_names=c("time 1 intercept","time 2 intercept","binary","linear")
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", ylab = "Est",xlab="",main= factor_names[i],ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  #abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}

df=results$ses
true_vals = df[2,]
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", ylab = "SE", xlab="",main= factor_names[i],ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  #abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}

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
