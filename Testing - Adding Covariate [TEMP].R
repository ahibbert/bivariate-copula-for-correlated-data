## This is a test script to run the models with covariates

source("common_functions.R")
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="NO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=.3,mu2=.7,dist="LO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="GA",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="PO",x1=1,x2=1)

calcTrueCovariateValues = function(n,a,b,c,mu1,mu2,dist,x1,x2) {
  
  #n=1000;a=1;b=1;c=.5;mu1=1;mu2=2;dist="NO";x1=1;x2=1;s1=s2=1
  #n=1000;a=1;b=1;c=.5;mu1=.3;mu2=.7;dist="LO";x1=1;x2=1;s1=s2=1
  
  dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2); 
  
  linkFunction=function(input_rv,dist="NO") {
    if(dist=="NO") {
      output_rv=input_rv  
    } else if (dist=="PO"|dist=="NB"|dist=="GA") {
      output_rv=log(input_rv)
    } else if (dist=="LO") {
      output_rv=logit(input_rv)
    }
    return(output_rv)
  }
  
  linkInvFunction=function(input_rv,dist="NO") {
    if(dist=="NO") {
      output_rv=input_rv  
    } else if (dist=="PO"|dist=="NB"|dist=="GA") {
      output_rv=exp(input_rv)
    } else if (dist=="LO") {
      output_rv=logit_inv(input_rv)
    }
    return(output_rv)
  }
  
  library(gamlss)
  pdf=if(dist=="NO") {dNO
  } else if (dist=="PO"){dNBI 
  } else if (dist=="GA"){dGA
  } else if (dist=="LO"){dLO}
  
  mu_1_eta=linkFunction(mu1,dist=dist)
  mu_2_eta=linkFunction(mu1,dist=dist)
  
  getLogLik=function(par_input,dataset,pdf,linkFunction,linkInvFunction,dist) {
    
    mu1=linkInvFunction(par_input[1],dist=dist);
    mu2=linkInvFunction(par_input[2],dist=dist);
    x1=par_input[3];
    x2=par_input[4];
    s1=exp(par_input[5]);
    s2=exp(par_input[6]);
    
    mean_estimates=cbind(dataset$time,linkInvFunction(dataset$sex*x1+dataset$age*x2-(dataset$time-1)*linkFunction(mu1,dist)+(dataset$time)*linkFunction(mu2,dist),dist))  

    if(dist=="LO") {
      time1=pdf(dataset$random_variable[dataset$time==0],mean_estimates[mean_estimates[,1]==0,2])
      time2=pdf(dataset$random_variable[dataset$time==1],mean_estimates[mean_estimates[,1]==1,2])
    } else {
      time1=pdf(dataset$random_variable[dataset$time==0],mean_estimates[mean_estimates[,1]==0,2],sigma=(s1))
      time2=pdf(dataset$random_variable[dataset$time==1],mean_estimates[mean_estimates[,1]==1,2],sigma=(s2))
    }
    
    return(-sum(log(time1))-sum(log(time2)))
  
  }
  #getLogLik(par=c(mu1_eta,mu2_eta,x1,x2,s1,s2),dataset,pdf,linkFunction,linkInvFunction,dist)
    
  #Write an optimisation that chooses optim_par such that getLogLik is maximised
  optim_par=optim(par=c(runif(1),runif(1),runif(1),runif(1),runif(1),runif(1))
                  , fn=getLogLik, dataset=dataset, pdf=pdf, linkFunction=linkFunction,linkInvFunction=linkInvFunction, dist=dist
                  , control=list(maxit=10000, reltol=1e-8, trace=0))
  
  return(optim_par)

}

simCovariateMLEs = function(sims,n,a,b,c,mu1,mu2,dist,x1,x2,trace) {
  optim_cov_outputs=matrix(NA,ncol=6,nrow=sims)
  counter=0
  while (counter < sims) {
    optim_est=calcTrueCovariateValues(n,a,b,c,mu1,mu2,dist,x1,x2)
    if(optim_est$convergence==0) {
      counter=counter+1
      optim_cov_outputs[counter,]=optim_est$par
    }
    if(trace==TRUE) {
      print(counter)  
    }
  }
  colnames(optim_cov_outputs)=c("mu1","mu2","x1","x2","s1","s2")
  return(sqrt(diag(cov(optim_cov_outputs[,1:4]))*n)/sqrt(n))
}

true=simCovariateMLEs(sims=1000,n=1000,a=.5,b=.5,c=.75,mu1=1,mu2=2,dist="NO",x1=1,x2=1,trace=TRUE)
true=simCovariateMLEs(sims=1000,n=1000,a=.5,b=.5,c=.75,mu1=.3,mu2=.7,dist="LO",x1=1,x2=1,trace=TRUE)
true=simCovariateMLEs(sims=1000,n=1000,a=.5,b=.5,c=.75,mu1=1,mu2=2,dist="GA",x1=1,x2=1,trace=TRUE)
true=simCovariateMLEs(sims=1000,n=1000,a=.5,b=.5,c=.75,mu1=1,mu2=2,dist="PO",x1=1,x2=1,trace=TRUE)



true=calcTrueCovariateValues(dataset,dist="PO",a=1,b=1,c=.5,mu1=1,mu2=2,x1=1,x2=1)

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
