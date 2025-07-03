## This is a test script to run the models with covariates

generateBivDist_withCov <- function(n,a,b,c,mu1,mu2,dist,x1,x2) {
  
  sex <- sample(0:1, n, replace=TRUE)
  age <- runif(n, min=-1, max=1)
  
  if(dist=="GA") {
    #Simulating bivariate random variable according to functional input
    #n=1000;a=1;b=1;c=0.5;mu1=2;mu2=3;dist="GA";x1=1;x2=1
    w<-rbeta(n,a,b)
    
    mu1_long=mu1+sex*x1+age*x2
    mu2_long=mu2+sex*x1+age*x2
    
    time_1<-w*rgamma(n,shape=a+b,scale=mu1_long)
    time_2<-w*rgamma(n,shape=a+b,scale=mu2_long)
    
  }
  if(dist=="NO") {
    require(MASS)
    normData<-mvrnorm(n,mu=c(mu1,mu2),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    margin_1<-normData[,1]
    margin_2<-normData[,2]
    
    #trt <- sample(0:1, length(margin_1), replace=TRUE)
    
    time_1=margin_1 + x1*sex + x2*age
    time_2=margin_2 + x1*sex + x2*age
    
  }
  
  if(dist=="LO") {
    
    require(MASS)
    a=1;b=1
    #c=0.5;n=1000 
    normData<-mvrnorm(n,mu=c(0,0),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    
    margin_1<-normData[,1] + x1*sex + x2*age
    margin_2<-normData[,2] + x1*sex + x2*age
    
    time_1<-as.numeric(pnorm(margin_1,mean=mean(margin_1),sd=sd(margin_1))<=mu1)
    time_2<-as.numeric(pnorm(margin_2,mean=mean(margin_2),sd=sd(margin_2))<=mu2)
    
  }
  
  if(dist=="PO") {
    
    #Compound multiple poisson of Stein & Juritz, 1987
    mixing_dist<-rgamma(n,shape=c,scale=b)
    
    mu1_long=mu1+sex*x1+age*x2
    mu2_long=mu2+sex*x1+age*x2
    
    time_1=vector(length = n) 
    time_2=vector(length = n) 
    for (i in 1:n) {
      time_1[i]=rpois(1,mu1_long[i]*mixing_dist[i])
    }
    for (i in 1:n) {
      time_2[i]=rpois(1,mu2_long[i]*mixing_dist[i])
    }
  }
  
  #Transforming data to format required for random effect models
  patient<-as.factor(seq(1:n))
  dataset<-as.data.frame(rbind(cbind(patient,time_1,0,sex,age)
                               ,cbind(patient,time_2,1,sex,age)))
  colnames(dataset)<-c("patient","random_variable","time","sex","age")
  
  dataset<-dataset[order(dataset$patient),]
  
  return(dataset)
}

fitBivModels_Bt_withCov <-function(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=FALSE) {
  
  #dist="NO"; include="ALL"; calc_actuals=FALSE
  
  n=nrow(dataset[dataset$time==0,])
  
  #Data Setup
  gamma_c_mu1<-dataset[dataset$time==0,]
  gamma_c_mu2<-dataset[dataset$time==1,]
  
  #Calculating actuals for parameters where available
  if(calc_actuals==FALSE) {actuals<-c(NA,NA,NA,NA,NA,NA,NA,NA)} else {
    
    library(e1071)
    
    if(dist=="GA"){
      actuals<-c( log(a*mu1)
                  , log(a*mu2)
                  , sqrt(1/a)/sqrt(n)
                  , (2*sqrt(1/a)-2*(log((1/(a*mu1+a*mu2)) * (sqrt(b)/(a*(a+b+1)) )+1)))/sqrt(n)
                  , NA
                  ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
                  ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
                  ,(skewness(gamma_c_mu1$random_variable)+skewness(gamma_c_mu2$random_variable))*10000/2
      )
    }
    if(dist=="LO"){
      actuals<-c( log(mu1/(1-mu1))
                  , log(mu2/(1-mu2))
                  , NA
                  , NA
                  , NA
                  ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
                  ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
                  ,(skewness(gamma_c_mu1$random_variable)+skewness(gamma_c_mu2$random_variable))*10000/2
      )
    }
    if(dist=="NO"){
      actuals<-c( 
        mu1
        , mu2
        , (a*sqrt(1-c^2))/sqrt(n)
        , sqrt(a^2+b^2-2*a*b*c)/sqrt(n)
        , NA
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
        ,(skewness(gamma_c_mu1$random_variable)+skewness(gamma_c_mu2$random_variable))*10000/2
      )
    }
    if(dist=="PO"){
      
      e_x1 = mu1*c*b
      e_x2 = mu2*c*b
      v_x1 = (((mu1^2)*(c*b^2)+(mu1*c*b))/((mu1*c*b)^2))
      v_x2 = (((mu2^2)*(c*b^2)+(mu2*c*b))/((mu2*c*b)^2))
      
      actuals<-c( 
        e_x1
        , e_x2
        , sqrt(v_x1)     /sqrt(n)
        , sqrt(
          (v_x2 + v_x1)
          - log((mu1*mu2*c)/(e_x1*e_x2))
        ) /sqrt(n)
        , NA
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
        ,(skewness(gamma_c_mu1$random_variable)+skewness(gamma_c_mu2$random_variable))*10000/2
      )
    }
    
  }
  
  #Non-GJRM model regression first as once GJRM is loaded it breaks GAMLSS 
  if(include=="ALL" || include=="non-GJRM" ) {
    
    require(gamlss)
    require(gee)
    require(lme4)
    require(MASS)
    require(gamlss.mx)
    library(mgcv)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="GA") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1)+as.factor(sex)+age, data=dataset, family=Gamma(link = "log"), maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+as.factor(sex)+age+random(as.factor(patient)), data=dataset, family=GA()) ))
      #model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=GA(), method=CG(1000))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                       , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
      
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset, family=Gamma(link="log"))))
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=Gamma(link="log"))
      
    }
    if(dist=="NO") {
      ############UPDATED TO WITH COVARIATES
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1)+as.factor(sex)+age, data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+as.factor(sex)+age+random(as.factor(patient)), data=dataset, family=NO())))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      model_lme4 <- lmer(formula=random_variable~as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)
      model_gamm = gamm(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=gaussian)
      
    }
    
    if(dist=="LO") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1)+as.factor(sex)+age, data=dataset, family=binomial, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+as.factor(sex)+age+random(as.factor(patient)), data=dataset, family=BI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer(formula=random_variable~as.factor(time==1) +as.factor(sex)+age+ (1|patient), data=dataset,family=binomial)
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=binomial)
    }
    
    if(dist=="PO"||dist=="NB") {
      invisible(capture.output(model_glm <- glm.nb(random_variable~as.factor(time==1)+as.factor(sex)+age, data=dataset, maxit=1000)))
      #invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
      library(geeM)
      model_gee<-geem(random_variable~as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, init.beta=model_glm$coefficients,
                      family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")
      
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+as.factor(sex)+age+random(as.factor(patient)), data=dataset, family=NBI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer.nb(formula=random_variable~as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=nb(link="log"))
      
    }
    
    results_table=list()
    
    results_table[[1]]=summary(model_glm)$coeff[,1:2]
    if (dist=="PO")
    {
      results_table[[2]]=cbind(summary(model_gee)[[1]],summary(model_gee)[[3]])
    } else {
      results_table[[2]]=summary(model_gee)$coeff[,c(1,4)]
    }
    results_table[[3]]=cbind(summary(model_re_nosig)[1:4],summary(model_re_nosig)[6:9])
    results_table[[4]]=cbind(summary(model_re_np)[1:4],summary(model_re_np)[8:11])
    results_table[[5]]=summary(model_lme4)$coefficients[,c(1,2)]
    results_table[[6]]=cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix)))
    
    logLiks=c(
      logLik(model_glm)
      , NA
      , logLik(model_re_nosig)
      , logLik(model_re_np)
      , logLik(model_lme4)
      , logLik(model_gamm$lme)
    )
    
    dfs=c( ####Come back to DF
      n*2-df.residual(model_glm)
      , n*2-df.residual(model_glm)+1
      , model_re_nosig$df.fit
      , model_re_np$df.fit
      , NA
      , NA
    )
    
  }
  
  #GJRM model regressions
  if(include=="ALL" || include=="GJRM" ) {
    
    require(GJRM)
    
    #Setting up GJRM equations
    eq.mu.1 <- formula(random_variable~as.factor(sex)+age)
    eq.mu.2 <- formula(random_variable.1~as.factor(sex)+age)
    fl <- list(eq.mu.1, eq.mu.2)
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    if(dist=="PO"){margin_dist="NBI"}
    if(dist=="LO"){margin_dist="logit"}
    
    copula_models=list()
    copula_models_results=list()
    i=1
    for (copula in c("C0","N","J0","G0","F","AMH","FGM","PL","HO","T")) {
      copula_models[[i]] <- gjrm(fl, margins = c(margin_dist,margin_dist), copula = copula, data=data.frame(gamma_c_mu1,gamma_c_mu2), model="B")
      copula_models_results[[i]]=list()
      cm=copula_models[[i]]
      copula_models_results[[i]][[1]]=summary(cm)$tableP1[,1:2]
      copula_models_results[[i]][[2]]=summary(cm)$tableP2[,1:2]
      copula_models_results[[i]][[3]]=logLik(cm)
      copula_models_results[[i]][[4]]=summary(cm)$t.edf #Number of parameters
      copula_models_results[[i]][[5]]=solve(cm$He) #Inverse of hessian - variance matrix
      i=i+1
    }
  }
  
  ########### 4. Combining results #########
  
  
  coefficients_table= rbind(
    results_table[[1]][,1]
    , results_table[[2]][,1]
    , results_table[[3]][,1]
    , results_table[[4]][,1]
    , results_table[[5]][,1]
    , results_table[[6]][,1]
  )
  ses_table= rbind(
    results_table[[1]][,2]
    , results_table[[2]][,2]
    , results_table[[3]][,2]
    , results_table[[4]][,2]
    , results_table[[5]][,2]
    , results_table[[6]][,2]
  )
  loglik_table= rbind(
    c(logLiks[1],dfs[1])
    , c(logLiks[2],dfs[2])
    , c(logLiks[3],dfs[3])
    , c(logLiks[4],dfs[4])
    , c(logLiks[5],dfs[5])
    , c(logLiks[6],dfs[6])
  )
  
  #standardising copula results for coeffs
  std_copula_models_results=std_copula_se_results=matrix(NA,ncol=4,nrow=10)
  std_copula_models_logliks=matrix(NA,ncol=2,nrow=10)
  
  for (i in 1:length(copula_models_results)) {
    temp_cop_results=c(copula_models_results[[i]][[1]][,1],copula_models_results[[i]][[2]][,1])
    num_cov=length(temp_cop_results)
    std_copula_models_results[i,c(1)]=temp_cop_results[c(1)] #Intercept
    std_copula_models_results[i,c(2)]=temp_cop_results[c(4)]-temp_cop_results[c(1)] #Time
    std_copula_models_results[i,c(3)]=(temp_cop_results[c(2)]+temp_cop_results[c(5)])/2 #Sex / binary
    std_copula_models_results[i,c(4)]=(temp_cop_results[c(3)]+temp_cop_results[c(6)])/2 #Age / linear
    
    std_copula_se_results[i,c(1)]=sqrt(copula_models_results[[i]][[5]][1,1])
    std_copula_se_results[i,c(2)]=sqrt(copula_models_results[[i]][[5]][4,4] + copula_models_results[[i]][[5]][1,1] + 2*copula_models_results[[i]][[5]][4,1])
    std_copula_se_results[i,c(3)]=sqrt(copula_models_results[[i]][[5]][2,1] + copula_models_results[[i]][[5]][5,1] + 2*copula_models_results[[i]][[5]][2,5])
    std_copula_se_results[i,c(4)]=sqrt(copula_models_results[[i]][[5]][3,1] + copula_models_results[[i]][[5]][6,1] + 2*copula_models_results[[i]][[5]][3,6])
    std_copula_models_logliks[i,1]=copula_models_results[[i]][[3]]
    std_copula_models_logliks[i,2]=copula_models_results[[i]][[4]]
    
    std_copula_se_results[is.na(std_copula_se_results)]<-0 #Replace NAs with 0s
  }
  coefficients_table=rbind(coefficients_table,std_copula_models_results)
  ses_table=rbind(ses_table,std_copula_se_results)
  loglik_table=rbind(loglik_table,std_copula_models_logliks)
  
  rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=c("glm" ,"gee","re_nosig","re_np","lme4","gamm"   ,"cop","cop_n","cop_j","cop_g","cop_f","cop_amh","cop_fgm","cop_pl","cop_h","cop_t")
  
  output_list=list(
    coefficients=coefficients_table
    , ses=ses_table
    , logliks=loglik_table
    , actuals=actuals
  )
  
  return(output_list)
  ###################### 5. End of function #########
  
}

#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=.5,mu2=.5,dist="LO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="NO",x1=1,x2=1)
#dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="GA",x1=1,x2=1)
dataset=generateBivDist_withCov(n=1000,a=1,b=1,c=.5,mu1=1,mu2=2,dist="PO",x1=1,x2=1)

#library(gamlss); fits=fitDist(dataset$random_variable,type="count")

results=fitBivModels_Bt_withCov(dataset,dist="PO",include="ALL",a=1,b=1,c=.5,mu1=1,mu2=2,calc_actuals=FALSE)

###TEMP PLOTTING

df=results$coefficients
plot.new(); par(mfrow=c(2,ncol(df)))
true_vals = c(1,1,1,1)
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", xlab = "Values", ylab = "Name",main= colnames(df)[i],ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  #abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}

df=results$ses
true_vals = df[2,]
for(i in 1:ncol(df)) {
  # Create the plot, with y-axis suppressed
  plot(df[,i]~seq_len(nrow(df)), xaxt = "n", xlab = "Values", ylab = "Name",main= colnames(df)[i],ylim=c(min(df[,i], true_vals[i]) - 0.5, max(df[,i], true_vals[i]) + 0.5))
  #abline(h=true_vals[i],col="red")
  # Add the custom y-axis with names
  axis(1, at = seq_len(nrow(df)), labels = rownames(df))
}

#############Proving we can combine multiple gammas
shape1=2
shape2=2
scale1=3
scale2=5

gamma1=rgamma(100000,shape=shape1,scale=scale1)
gamma2=rgamma(100000,shape=shape2,scale=scale2)

meanguess=(shape1*scale1 + shape2*scale2)/2
varguess=(scale2*shape2*(scale2^2)+scale1*shape1*(scale1^2))/(scale1+scale2)
scaleguess=varguess/meanguess
shapeguess=(meanguess^2)/varguess

matrix(data=c(
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
  ,(mean(c(gamma1,gamma2))^2)/var(c(gamma1,gamma2))/shapeguess
  ,var(c(gamma1,gamma2))/mean(c(gamma1,gamma2))/scaleguess
  ,mean(c(gamma1,gamma2))/meanguess #shape*scale=mean of 9 - so we add the scales and divide by 2
  ,var(c(gamma1,gamma2))/varguess #This is the scale parameter (4.3????)
  
),ncol=4,byrow=TRUE,dimnames=list(c(paste("gamma 1:"),paste("gamma 2:"),paste("gamma (1,2) calc"),"Par method","Diff"),c("shape","scale","mean","var")))

###VARIANCE GUESS IS CORRECT OK SO WE CAN ESTIMATE SCALE


##########OK What about negbin 





