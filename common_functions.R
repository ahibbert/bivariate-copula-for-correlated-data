logit <- function(x) {
  return(log(x/(1-x)))
}
#' @export
logit_inv <- function(x) {
  return(
    if(all(is.nan(exp(x)/(1+exp(x))))) {
      return(1)
    } else {
      return(exp(x)/(1+exp(x)))
    }
  )
}

generateBivDist <- function(n,a,b,c,mu1,mu2,dist) {
  
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
  dataset<-as.data.frame(rbind(cbind(patient,margin_1,0)
                               ,cbind(patient,margin_2,1)))
  colnames(dataset)<-c("patient","random_variable","time")
  
  dataset<-dataset[order(dataset$patient),]
  
  return(dataset)
}

fitBivModels <-function(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE) {
  
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
                  , NA
                  , NA
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
        , (b*sqrt(1-c^2))/sqrt(n)
        , sqrt(a^2+b^2-2*a*b*c)/sqrt(n)
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
        , sqrt(v_x2)     /sqrt(n)
        , sqrt(
            (v_x2 + v_x1)
            - log((mu1*mu2*c)/(e_x1*e_x2))
        ) /sqrt(n)
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
        ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
        ,(skewness(gamma_c_mu1$random_variable)+skewness(gamma_c_mu2$random_variable))*10000/2
      )
    }
    
  }
  
  if(include=="ALL" || include=="non-GJRM" ) {

    library(gee)
    library(gamlss)
    library(lme4)
    library(gamlss.mx)
    library(mgcv)
    library(MASS)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="GA") {
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1), data=dataset, family=Gamma(link = "log"), maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA()) ))
      #model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=GA(), method=CG(1000))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                              , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
      
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))))
    
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1), random=list(patient=~1), data=dataset, family=Gamma(link="log"))
      
    }
    if(dist=="LO") {
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1), data=dataset, family=binomial, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=BI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                              , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset,family=binomial)
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1), random=list(patient=~1), data=dataset, family=binomial)
    }
    if(dist=="NO") {
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1), data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                              , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset)
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1), random=list(patient=~1), data=dataset, family=gaussian)
    }
    
    if(dist=="PO"||dist=="NB") {
      
      invisible(capture.output(model_glm <- glm.nb(random_variable~-1+as.factor(time==1), data=dataset, maxit=1000)))
      #invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
      library(geeM)
      model_gee<-geem(random_variable~-1+as.factor(time==1), id=patient, data=dataset, init.beta=model_glm$coefficients,
                       family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")

      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NBI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                              , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer.nb(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset)
      
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1), random=list(patient=~1), data=dataset, family=nb(link="log"))
    }
    
    ###Capturing coefficient values and errors from each model
    summary_glm<-c( summary(model_glm)$coeff[1]
                    ,summary(model_glm)$coeff[2]
                    ,summary(model_glm)$coeff[3]
                    ,summary(model_glm)$coeff[4]
                    , logLik(model_glm)
                    , AIC(model_glm)
                    , BIC(model_glm)
                    , 3
    )
    if(dist=="PO"||dist=="NB"){model_sum_gee=summary(model_gee)
                                summary_gee<-c( 
                                 model_sum_gee$beta[1]
                               , model_sum_gee$beta[2]
                               , model_sum_gee$se.model[1]#### 
                               , model_sum_gee$se.model[2]####
                               , NA
                               , NA
                               , NA
                               , 4)} else{
    summary_gee<-c( summary(model_gee)$coeff[1]
                    , summary(model_gee)$coeff[2]
                    , summary(model_gee)$coeff[7]#### 
                    , summary(model_gee)$coeff[8]####
                    , NA
                    , NA
                    , NA
                    , 4
    )}
    
    invisible(capture.output(
      summary_re_nosig<-c( summary(model_re_nosig)[1]
                           ,summary(model_re_nosig)[2]
                           ,if(dist=="LO"){summary(model_re_nosig)[3]}else{summary(model_re_nosig)[4]}
                           ,if(dist=="LO"){summary(model_re_nosig)[4]}else{summary(model_re_nosig)[5]}
                           ,logLik(model_re_nosig)
                           ,AIC(model_re_nosig)
                           ,BIC(model_re_nosig)
                           , model_re_nosig$df.fit
      )
    ))
    
    #invisible(capture.output(
      #summary_re<-c( summary(model_re)[1]
      #               ,summary(model_re)[2]
      #               ,if(dist=="PO"){summary(model_re)[3]}else{summary(model_re)[5]}
      #               ,if(dist=="PO"){summary(model_re)[4]}else{summary(model_re)[6]}
      #               , AIC(model_re)
      #               , BIC(model_re)
      #               , model_re$df.fit
      #)
    #))
    
    invisible(capture.output(
      summary_re_np<-c( summary(model_re_np)[1]
                     ,summary(model_re_np)[2]
                     ,if(dist=="LO"){summary(model_re_np)[5]}else{summary(model_re_np)[6]}
                     ,if(dist=="LO"){summary(model_re_np)[6]}else{summary(model_re_np)[7]}
                     ,logLik(model_re_np)
                     , AIC(model_re_np)
                     , BIC(model_re_np)
                     , model_re_np$df.fit
      )
    ))
    
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      , AIC(model_lme4)
                      , BIC(model_lme4)
                      , 4)
    
    ###Calculating effective degrees of freedom from Donohue
    X<-getME(model_lme4,name="X")
    Z<-getME(model_lme4,name="Z")
    U<-cbind(X,Z)
    W<-model_lme4@resp$sqrtrwt #weights(model_lme4,type = "working")
    UWU=(t(as.matrix(U))%*%(diag(as.vector(W)))%*%as.matrix(U))
    dim(UWU)
    D<-getME(model_lme4,name="Lambda")
    
    if(sum(D)==0) {lme_EDF=summary_lme4[length(summary_lme4)]} else {
      D_inv<-solve(D)
      dinv_plus_00<-c(0,0,diag(D_inv))
      lme_EDF=sum(diag(UWU%*%solve(UWU+diag(dinv_plus_00))))
      
    }
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      ,-2*logLik(model_lme4)+2*lme_EDF
                      , BIC(model_lme4)
                      ,lme_EDF)
    
    summary_gamm<-c( summary(model_gamm$lme)$coefficients[[1]][1]
                     , summary(model_gamm$lme)$coefficients[[1]][2]
                     , sqrt(diag(model_gamm$lme$varFix))[1]
                     , sqrt(diag(model_gamm$lme$varFix))[2]
                     ,logLik(model_gamm$lme)
                     ,AIC(model_gamm$lme)
                     ,BIC(model_gamm$lme)
                     ,lme_EDF
    )
    
  }
  
  if(include=="ALL" || include=="GJRM" ) {
  
    require(GJRM)
    
    #Setting up GJRM equations
    eq.mu.1 <- formula(random_variable~1)
    eq.mu.2 <- formula(random_variable.1~1)
    fl <- list(eq.mu.1, eq.mu.2)
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    if(dist=="PO"){margin_dist="NBI"}
    if(dist=="LO"){margin_dist="logit"}

    model_copula<-    gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_j<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_g<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "G0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_f<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "F",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_amh<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "AMH",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_fgm<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "FGM",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_pl<- gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "PL",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_h<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "HO",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_t<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "T",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    
    summary_cop<-c( model_copula$coefficients[1]
                    , model_copula$coefficients[2]
                    , summary(model_copula)$tableP1[2] #SE for time 0
                    , summary(model_copula)$tableP2[2] #SE for time 1
                    ,logLik(model_copula)
                    , 2*5-2*logLik(model_copula)
                    ,BIC(model_copula)
                    , 5
    )
    
    summary_cop_n<-c( model_copula_n$coefficients[1]
                      , model_copula_n$coefficients[2] 
                      , summary(model_copula_n)$tableP1[2] #SE for time 0
                      , summary(model_copula_n)$tableP2[2] #SE for time 1
                      , logLik(model_copula_n)
                      , 2*5-2*logLik(model_copula_n)
                      ,BIC(model_copula_n)
                      , 5
                      
    )
    summary_cop_j<-c( model_copula_j$coefficients[1]
                      , model_copula_j$coefficients[2]
                      , summary(model_copula_j)$tableP1[2] #SE for time 0
                      , summary(model_copula_j)$tableP2[2] #SE for time 1
                      , logLik(model_copula_j)
                      , 2*5-2*logLik(model_copula_j)
                      ,BIC(model_copula_j)
                      , 5
                      
    )
    summary_cop_g<-c( model_copula_g$coefficients[1]
                      , model_copula_g$coefficients[2] 
                      , summary(model_copula_g)$tableP1[2] #SE for time 0
                      , summary(model_copula_g)$tableP2[2] #SE for time 1
                      , logLik(model_copula_g)
                      , 2*5-2*logLik(model_copula_g)
                      ,BIC(model_copula_g)
                      , 5
                      
    )
    summary_cop_f<-c( model_copula_f$coefficients[1]
                      , model_copula_f$coefficients[2]
                      , summary(model_copula_f)$tableP1[2] #SE for time 0
                      , summary(model_copula_f)$tableP2[2] #SE for time 1
                      , logLik(model_copula_f)
                      , 2*5-2*logLik(model_copula_f)
                      ,BIC(model_copula_f)
                      , 5
                      
    )
    summary_cop_amh<-c( model_copula_amh$coefficients[1]
                        , model_copula_amh$coefficients[2]
                        , summary(model_copula_amh)$tableP1[2] #SE for time 0
                        , summary(model_copula_amh)$tableP2[2] #SE for time 1
                        , logLik(model_copula_amh)
                        , 2*5-2*logLik(model_copula_amh)
                        ,BIC(model_copula_amh)
                        , 5
                        
    )
    summary_cop_fgm<-c( model_copula_fgm$coefficients[1]
                        , model_copula_fgm$coefficients[2]
                        , summary(model_copula_fgm)$tableP1[2] #SE for time 0
                        , summary(model_copula_fgm)$tableP2[2] #SE for time 1
                        , logLik(model_copula_fgm)
                        , 2*5-2*logLik(model_copula_fgm)
                        ,BIC(model_copula_fgm)
                        , 5
                        
    )
    summary_cop_pl<-c( model_copula_pl$coefficients[1]
                       , model_copula_pl$coefficients[2]
                       , summary(model_copula_pl)$tableP1[2] #SE for time 0
                       , summary(model_copula_pl)$tableP2[2] #SE for time 1
                       , logLik(model_copula_pl)
                       , 2*5-2*logLik(model_copula_pl)
                       ,BIC(model_copula_pl)
                       , 5
                       
    )
    summary_cop_h<-c( model_copula_h$coefficients[1]
                      , model_copula_h$coefficients[2]
                      , summary(model_copula_h)$tableP1[2] #SE for time 0
                      , summary(model_copula_h)$tableP2[2] #SE for time 1
                      , logLik(model_copula_h)
                      , 2*5-2*logLik(model_copula_h)
                      ,BIC(model_copula_h)
                      , 5
                      
    )
    summary_cop_t<-c( model_copula_t$coefficients[1]
                      , model_copula_t$coefficients[2]
                      , summary(model_copula_t)$tableP1[2] #SE for time 0
                      , summary(model_copula_t)$tableP2[2] #SE for time 1
                      , logLik(model_copula_t)
                      , AIC(model_copula_t)
                      ,BIC(model_copula_t)
                      , 6
                      
    )
  }

  ########### 4. Combining results #########
  
  if(include=="ALL") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,summary_gamm,summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)  
  }
  if(include=="GJRM") {
    results_exbias<- rbind(summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)
  }
  if(include=="non-GJRM") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,summary_gamm,actuals)
  }
  
  results_exbias[,1:4]<- round(results_exbias[,1:4],4)
  results_exbias[,5:8]<- round(results_exbias[,5:8],0)
  
  colnames(results_exbias)<-c("b_1","b_2","se_b1","se_b2","LogLik","AIC","BIC","EDF")
  
  return(results_exbias)
  
}

fitBivModels_Bt <-function(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE) {
  
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
  
  if(include=="ALL" || include=="non-GJRM" ) {
    
    require(gamlss)
    require(gee)
    require(lme4)
    require(MASS)
    require(gamlss.mx)
    library(mgcv)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="GA") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=Gamma(link = "log"), maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA()) ))
      #model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=GA(), method=CG(1000))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                       , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
      
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))))
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=Gamma(link="log"))
      
    }
    if(dist=="NO") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- lmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=gaussian)
    }
    
    if(dist=="LO") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=binomial, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=BI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset,family=binomial)
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=binomial)
    }
    
    if(dist=="PO"||dist=="NB") {
      invisible(capture.output(model_glm <- glm.nb(random_variable~as.factor(time==1), data=dataset, maxit=1000)))
      #invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
      library(geeM)
      model_gee<-geem(random_variable~as.factor(time==1), id=patient, data=dataset, init.beta=model_glm$coefficients,
                      family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")
      
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NBI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer.nb(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)
      
      model_gamm = gamm(formula=random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=nb(link="log"))
      
    }
    
    ###Capturing coefficient values and errors from each model
    summary_glm<-c( summary(model_glm)$coeff[1]
                    ,summary(model_glm)$coeff[2]
                    ,summary(model_glm)$coeff[3]
                    ,summary(model_glm)$coeff[4]
                    , logLik(model_glm)
                    , AIC(model_glm)
                    , BIC(model_glm)
                    , 3
    )
    
    if(dist=="PO"||dist=="NB"){model_sum_gee=summary(model_gee)
    summary_gee<-c( 
      model_sum_gee$beta[1]
      , model_sum_gee$beta[2]
      , model_sum_gee$se.model[1]#### 
      , model_sum_gee$se.model[2]####
      , NA
      , NA
      , NA
      , 4)} else{
        summary_gee<-c( summary(model_gee)$coeff[1]
                        , summary(model_gee)$coeff[2]
                        , summary(model_gee)$coeff[7]#### 
                        , summary(model_gee)$coeff[8]####
                        , NA
                        , NA
                        , NA
                        , 4
        )}
    
    invisible(capture.output(
      summary_re_nosig<-c( summary(model_re_nosig)[1]
                           ,summary(model_re_nosig)[2]
                           ,if(dist=="LO"){summary(model_re_nosig)[3]}else{summary(model_re_nosig)[4]}
                           ,if(dist=="LO"){summary(model_re_nosig)[4]}else{summary(model_re_nosig)[5]}
                           ,logLik(model_re_nosig)
                           ,AIC(model_re_nosig)
                           ,BIC(model_re_nosig)
                           , model_re_nosig$df.fit
      )
    ))
    
    #invisible(capture.output(
    #summary_re<-c( summary(model_re)[1]
    #               ,summary(model_re)[2]
    #               ,if(dist=="PO"){summary(model_re)[3]}else{summary(model_re)[5]}
    #               ,if(dist=="PO"){summary(model_re)[4]}else{summary(model_re)[6]}
    #               , AIC(model_re)
    #               , BIC(model_re)
    #               , model_re$df.fit
    #)
    #))
    
    invisible(capture.output(
      summary_re_np<-c( summary(model_re_np)[1]
                        ,summary(model_re_np)[2]
                        ,if(dist=="LO"){summary(model_re_np)[5]}else{summary(model_re_np)[6]}
                        ,if(dist=="LO"){summary(model_re_np)[6]}else{summary(model_re_np)[7]}
                        ,logLik(model_re_np)
                        , AIC(model_re_np)
                        , BIC(model_re_np)
                        , model_re_np$df.fit
      )
    ))
    
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      , AIC(model_lme4)
                      , BIC(model_lme4)
                      , 4)
    
    
    ###Calculating effective degrees of freedom from Donohue
    X<-getME(model_lme4,name="X")
    Z<-getME(model_lme4,name="Z")
    U<-cbind(X,Z)
    W<-model_lme4@resp$sqrtrwt #weights(model_lme4,type = "working")
    UWU=(t(as.matrix(U))%*%(diag(as.vector(W)))%*%as.matrix(U))
    dim(UWU)
    D<-getME(model_lme4,name="Lambda")
    
    if(sum(D)==0) {lme_EDF=summary_lme4[length(summary_lme4)]} else {
      D_inv<-solve(D)
      dinv_plus_00<-c(0,0,diag(D_inv))
      lme_EDF=sum(diag(UWU%*%solve(UWU+diag(dinv_plus_00))))
      
    }
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      ,-2*logLik(model_lme4)+2*lme_EDF
                      , BIC(model_lme4)
                      ,lme_EDF)
    
  
    summary_gamm<-c( summary(model_gamm$lme)$coefficients[[1]][1]
                     , summary(model_gamm$lme)$coefficients[[1]][2]
                     , sqrt(diag(model_gamm$lme$varFix))[1]
                     , sqrt(diag(model_gamm$lme$varFix))[2]
                     ,logLik(model_gamm$lme)
                     ,AIC(model_gamm$lme)
                     ,BIC(model_gamm$lme)
                     ,lme_EDF
    )
    
  }
  
  if(include=="ALL" || include=="GJRM" ) {
    
    require(GJRM)
    
    #Setting up GJRM equations
    eq.mu.1 <- formula(random_variable~1)
    eq.mu.2 <- formula(random_variable.1~1)
    fl <- list(eq.mu.1, eq.mu.2)
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    if(dist=="PO"){margin_dist="NBI"}
    if(dist=="LO"){margin_dist="logit"}
    
    model_copula<-    gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_j<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_g<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "G0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_f<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "F",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_amh<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "AMH",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_fgm<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "FGM",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_pl<- gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "PL",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_h<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "HO",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_t<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "T",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    
    slv_cop_He=solve(model_copula$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop<-c( model_copula$coefficients[1]
                    , model_copula$coefficients[2]
                    , summary(model_copula)$tableP1[2] #SE for time 0
                    , se_2
                    ,logLik(model_copula)
                    , 2*5-2*logLik(model_copula)
                    ,BIC(model_copula)
                    , 5
    )
    slv_cop_He=solve(model_copula_n$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_n<-c( model_copula_n$coefficients[1]
                      , model_copula_n$coefficients[2] 
                      , summary(model_copula_n)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_n)
                      , 2*5-2*logLik(model_copula_n)
                      ,BIC(model_copula_n)
                      , 5
                      
    )
    slv_cop_He=solve(model_copula_j$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_j<-c( model_copula_j$coefficients[1]
                      , model_copula_j$coefficients[2]
                      , summary(model_copula_j)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_j)
                      , 2*5-2*logLik(model_copula_j)
                      ,BIC(model_copula_j)
                      , 5
                      
    )
    slv_cop_He=solve(model_copula_g$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_g<-c( model_copula_g$coefficients[1]
                      , model_copula_g$coefficients[2] 
                      , summary(model_copula_g)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_g)
                      , 2*5-2*logLik(model_copula_g)
                      ,BIC(model_copula_g)
                      , 5
                      
    )
    slv_cop_He=solve(model_copula_f$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_f<-c( model_copula_f$coefficients[1]
                      , model_copula_f$coefficients[2]
                      , summary(model_copula_f)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_f)
                      , 2*5-2*logLik(model_copula_f)
                      ,BIC(model_copula_f)
                      , 5
                      
    )
    slv_cop_He=solve(model_copula_amh$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_amh<-c( model_copula_amh$coefficients[1]
                        , model_copula_amh$coefficients[2]
                        , summary(model_copula_amh)$tableP1[2] #SE for time 0
                        , se_2
                        , logLik(model_copula_amh)
                        , 2*5-2*logLik(model_copula_amh)
                        ,BIC(model_copula_amh)
                        , 5
                        
    )
    slv_cop_He=solve(model_copula_fgm$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_fgm<-c( model_copula_fgm$coefficients[1]
                        , model_copula_fgm$coefficients[2]
                        , summary(model_copula_fgm)$tableP1[2] #SE for time 0
                        , se_2
                        , logLik(model_copula_fgm)
                        , 2*5-2*logLik(model_copula_fgm)
                        ,BIC(model_copula_fgm)
                        , 5
                        
    )
    slv_cop_He=solve(model_copula_pl$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_pl<-c( model_copula_pl$coefficients[1]
                       , model_copula_pl$coefficients[2]
                       , summary(model_copula_pl)$tableP1[2] #SE for time 0
                       , se_2
                       , logLik(model_copula_pl)
                       , 2*5-2*logLik(model_copula_pl)
                       ,BIC(model_copula_pl)
                       , 5
                       
    )
    slv_cop_He=solve(model_copula_h$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_h<-c( model_copula_h$coefficients[1]
                      , model_copula_h$coefficients[2]
                      , summary(model_copula_h)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_h)
                      , 2*5-2*logLik(model_copula_h)
                      ,BIC(model_copula_h)
                      , 5
                      
    )
    slv_cop_He=solve(model_copula_t$He);se_2=sqrt(slv_cop_He[1,1]+slv_cop_He[2,2]-2*slv_cop_He[1,2])
    summary_cop_t<-c( model_copula_t$coefficients[1]
                      , model_copula_t$coefficients[2]
                      , summary(model_copula_t)$tableP1[2] #SE for time 0
                      , se_2
                      , logLik(model_copula_t)
                      , AIC(model_copula_t)
                      ,BIC(model_copula_t)
                      , 6
                      
    )
  }
  
  ########### 4. Combining results #########
  
  if(include=="ALL") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,summary_gamm,summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)  
  }
  if(include=="GJRM") {
    results_exbias<- rbind(summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)
  }
  if(include=="non-GJRM") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,summary_gamm,actuals)
  }
  
  results_exbias[,1:4]<- round(results_exbias[,1:4],4)
  results_exbias[,5:8]<- round(results_exbias[,5:8],0)
  
  colnames(results_exbias)<-c("b_1","b_2","se_b1","se_b2","LogLik","AIC","BIC","EDF")
  
  return(results_exbias)
  
}

generateBivDist_withCov <- function(n,a,b,c,mu1,mu2,dist,x1,x2) {
  
  if((n%%10)!=0) {
    stop("n must be a multiple of 100")
  }
  sex <- c(rep(0,n/2),rep(1,n/2))
  age <- rep(1:10,n/10)
  
  if(dist=="GA") {
    #Simulating bivariate random variable according to functional input
    #n=1000;a=1;b=1;c=0.5;mu1=2;mu2=3;dist="GA";x1=1;x2=1
    w<-rbeta(n,a,b)
    
    mu1_long=exp(log(mu1)+sex*x1+age*x2)
    mu2_long=exp(log(mu2)+sex*x1+age*x2)
    
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
    
    margin_1<-pnorm(normData[,1],mean=0,sd=a) 
    margin_2<-pnorm(normData[,2],mean=0,sd=b)
  
    #time_1<-as.numeric(pnorm(margin_1,mean=mean(margin_1),sd=sd(margin_1))<=mu1)
    #time_2<-as.numeric(pnorm(margin_2,mean=mean(margin_2),sd=sd(margin_2))<=mu2)
    
    time_1<-as.numeric(logit_inv(logit(margin_1) - x1*sex - x2*age)<=mu1)
    time_2<-as.numeric(logit_inv(logit(margin_2) - x1*sex - x2*age)<=mu2)
    
  }
  if(dist=="PO") {
    
    #Compound multiple poisson of Stein & Juritz, 1987
    mixing_dist<-rgamma(n,shape=c,scale=b)
    
    mu1_long=exp(log(mu1)+sex*x1+age*x2)
    mu2_long=exp(log(mu2)+sex*x1+age*x2)
    
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

fitBivModels_Bt_withCov <-function(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE,cv=FALSE) {
  
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
    library(glmtoolbox)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="GA") {
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=GA)))
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=Gamma(link = "log"), maxit=1000)))
      invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~-1+as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                                                       , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset, family=Gamma(link="log"))))
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=Gamma(link="log"))
    } else if(dist=="NO") {
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NO)))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~-1+as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      invisible(capture.output(model_lme4 <- lmer(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)))
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=gaussian)
    } else if(dist=="PO") {
      invisible(capture.output(model_glm <- glm.nb(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, maxit=1000)))
      #invisible(capture.output(model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
      model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, init.beta=model_glm$coefficients,
                        family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable")
      #invisible(capture.output(model_gee<-overglm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family="nb1(log)", maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~-1+as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      model_lme4 <- glmer.nb(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age + (1|patient), data=dataset)
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=nb(link="log"),control=list(niterEM=1000))
    } else if (dist=="LO") {
      invisible(capture.output(model_glm <- glm(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, data=dataset, family=binomial, maxit=1000)))
      invisible(capture.output(model_gee<-glmgee(random_variable~-1+as.factor(time==1)+as.factor(sex)+age, id=patient, data=dataset, family=binomial, maxiter=25, corstr = "exchangeable")))
      #model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=BI)
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age+re(random=~1|patient), data=dataset, family=NBI)))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, sigma.formula=~-1+as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= BI()
                                                       , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) +as.factor(sex)+age+ (1|patient), data=dataset,family=binomial)))
      
      model_gamm = gamm(formula=random_variable~-1+as.factor(time==1)+as.factor(sex)+age, random=list(patient=~1), data=dataset, family=binomial)
    }
    
    
    results_table=list()
    
    results_table[[1]]=summary(model_glm)$coeff[,1:2]
    results_table[[2]]=summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2]
    results_table[[3]]=cbind(summary(model_re_nosig)[1:4],summary(model_re_nosig)[6:9])
    results_table[[4]]=cbind(summary(model_re_np)[1:4],summary(model_re_np)[8:11])
    results_table[[5]]=summary(model_lme4)$coefficients[,c(1,2)]
    results_table[[6]]=cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix)))
    
    logLiks=c(
      logLik(model_glm)
      , model_gee$logLik
      , logLik(model_re_nosig)
      , logLik(model_re_np)
      , logLik(model_lme4)
      , logLik(model_gamm$lme)
    )
    
    dfs=c( ####Come back to DF
      n*2-df.residual(model_glm)
      , n*2-model_gee$df.residual
      , model_re_nosig$df.fit
      , model_re_np$df.fit
      , NA
      , NA
    )
    
    #Creating summary tables
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
    
    sigmas=rbind(
      rep(summary(model_glm)$dispersion,2)
      , rep(model_gee$phi,2)
      , rep((model_re_nosig$sigma.coefficients),2)
      , model_re_np$sigma.coefficients
      , rep(summary(model_lme4)$sigma,2)
      , rep(model_gamm$lme$sigma,2)
    )
    
    #Extract correlations - for random effect models this is the random effect sd
    correlations=rbind(
    0
    ,(model_gee$corr[2,1])
    ,as.numeric(VarCorr(getSmo(model_re_nosig))[[1]]) 
    ,0
    ,summary(model_lme4)$varcor$patient[1,1]
    ,var(ranef(model_gamm$lme)[[1]])
    )
    rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=
      rownames(sigmas)=rownames(correlations)=c("glm" ,"gee","re_nosig","re_np","lme4","gamm")
    
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
      copula_models[[i]] <- gjrm(fl, margins = c(margin_dist,margin_dist), copula = copula, data=data.frame(dataset[dataset$time==0,],dataset[dataset$time==1,]), model="B")
      copula_models_results[[i]]=list()
      cm=copula_models[[i]]
      copula_models_results[[i]][[1]]=summary(cm)$tableP1[,1:2]
      copula_models_results[[i]][[2]]=summary(cm)$tableP2[,1:2]
      copula_models_results[[i]][[3]]=logLik(cm)
      copula_models_results[[i]][[4]]=summary(cm)$t.edf #Number of parameters
      copula_models_results[[i]][[5]]=solve(cm$He) #Inverse of hessian - variance matrix
      copula_models_results[[i]][[6]]=c(cm$sigma1,cm$sigma2) #Sigmas)
      copula_models_results[[i]][[7]]=cm$theta #Thetas)
      i=i+1
    }
    
    #standardising copula results for coeffs
    std_copula_models_results=std_copula_se_results=matrix(NA,ncol=4,nrow=10)
    std_copula_models_logliks=matrix(NA,ncol=2,nrow=10)
    std_copula_models_sigmas=matrix(NA,ncol=2,nrow=10)
    std_copula_models_correlations=matrix(NA,ncol=1,nrow=10)
    
    for (i in 1:length(copula_models_results)) {
      temp_cop_results=c(copula_models_results[[i]][[1]][,1],copula_models_results[[i]][[2]][,1])
      num_cov=length(temp_cop_results)
      std_copula_models_results[i,c(1)]=temp_cop_results[c(1)] #Intercept
      std_copula_models_results[i,c(2)]=temp_cop_results[c(4)]#-temp_cop_results[c(1)] #Time
      std_copula_models_results[i,c(3)]=(temp_cop_results[c(2)]+temp_cop_results[c(5)])/2 #Sex / binary
      std_copula_models_results[i,c(4)]=(temp_cop_results[c(3)]+temp_cop_results[c(6)])/2 #Age / linear
      
      std_copula_se_results[i,c(1)]=sqrt(copula_models_results[[i]][[5]][1,1])
      std_copula_se_results[i,c(2)]=sqrt(copula_models_results[[i]][[5]][4,4]) #+ copula_models_results[[i]][[5]][1,1] + 2*copula_models_results[[i]][[5]][4,1])
      std_copula_se_results[i,c(3)]=sqrt((copula_models_results[[i]][[5]][2,2] + copula_models_results[[i]][[5]][5,5] + 2*copula_models_results[[i]][[5]][2,5]))/2
      std_copula_se_results[i,c(4)]=sqrt((copula_models_results[[i]][[5]][3,3] + copula_models_results[[i]][[5]][6,6] + 2*copula_models_results[[i]][[5]][3,6]))/2
      std_copula_models_logliks[i,1]=copula_models_results[[i]][[3]]
      std_copula_models_logliks[i,2]=copula_models_results[[i]][[4]]
      std_copula_models_sigmas[i,c(1,2)]=copula_models_results[[i]][[6]]
      std_copula_models_correlations[i,1]=copula_models_results[[i]][[7]]
      
      std_copula_se_results[is.na(std_copula_se_results)]<-0 #Replace NAs with 0s
    }
  }
  
  ########### 4. Combining results #########
  
  if(include=="ALL") {
    coefficients_table=rbind(coefficients_table,std_copula_models_results)
    ses_table=rbind(ses_table,std_copula_se_results)
    loglik_table=rbind(loglik_table,std_copula_models_logliks)
    sigmas=rbind(sigmas,std_copula_models_sigmas)
    correlations=rbind(correlations,std_copula_models_correlations)
    rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=rownames(sigmas)=rownames(correlations)=c("glm" ,"gee","re_nosig","re_np","lme4","gamm"   ,"cop","cop_n","cop_j","cop_g","cop_f","cop_amh","cop_fgm","cop_pl","cop_h","cop_t")
  } else if (include == "GJRM") {
    coefficients_table=std_copula_models_results
    ses_table=std_copula_se_results
    loglik_table=std_copula_models_logliks
    sigmas=std_copula_models_sigmas
    correlations=std_copula_models_correlations
    rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=rownames(sigmas)=rownames(correlations)=c("cop","cop_n","cop_j","cop_g","cop_f","cop_amh","cop_fgm","cop_pl","cop_h","cop_t")
  } else if (include == "non-GJRM") {
    coefficients_table=coefficients_table
    ses_table=ses_table
    loglik_table=loglik_table
    rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=rownames(sigmas)=rownames(correlations)=c("glm" ,"gee","re_nosig","re_np","lme4","gamm")
  }
    
  output_list=list(
    coefficients=coefficients_table
    , ses=ses_table
    , logliks=loglik_table
    , actuals=actuals
    , sigmas=sigmas
    , correlations=correlations
    , y=dataset$random_variable
    , dist=dist
  )
  
  return(output_list)
  ###################### 5. End of function #########
  
}

# simRE = function (n,dist,mu1,mu2, x1, x2, sigma1,sigma2,b_var,sims=1000) {
#   
#   #n=n;dist=dist;mu1=mu1_in;mu2=mu2_in;x1=x1_in; x2=x2_in;sigma1=sigma1_in;sigma2=sigma2_in;b_var=b_var_in;sims=100
#   
#   library(gamlss)
#   
#   #sex <- c(rep(0,n/2),rep(1,n/2))
#   #age <- rep(1:10,n/10)
#   
#   out_outer=matrix(ncol=5,nrow=sims)
#   colnames(out_outer)=c("mean_1","mean_2","var_1","var_2","corr")
#   
#   for (i in 1:sims) {
#     if(dist=="GA") {
#       re_val=rnorm(n,mean=0,sd=sqrt(sqrt(b_var)))
#       t1=rGA(n,mu=exp(log(mu1)),sigma=sqrt(sigma1))
#       t2=rGA(n,mu=exp(log(mu2)),sigma=sqrt(sigma2))
#     }
#     
#     #mean - mu1
#     #sigma_1^2 * mu1^2 - variance
#     
#     t1_plus_re=exp(log(t1)+re_val)
#     t2_plus_re=exp(log(t2)+re_val)
#     
#     out_outer[i,]=c(mean(t1_plus_re)
#                     ,mean(t2_plus_re)
#                     ,var(t1_plus_re)
#                     ,var(t2_plus_re)
#                     ,cor(t1_plus_re,t2_plus_re))
#   }
#   
#   medians=c(median(out_outer[,1])
#             ,median(out_outer[,2])
#             ,(median(out_outer[,3]))
#             ,(median(out_outer[,4]))
#             ,median(out_outer[,5]))
#   names(medians)=colnames(out_outer)
#   print(medians)
#   print(colMeans(out_outer))
#   
#   #return(colMeans(out_outer))
#   return(medians)
# }

sim_model <- function(model,dist,n,coefficients,sigmas,correlations) {
  
  sex <- c(rep(0,n/2),rep(1,n/2))
  age <- rep(1:10,n/10)
  
  #Calculate linear predictor
  lp_1=coefficients[model,1]+coefficients[model,3]*sex + coefficients[model,4]*age
  lp_2=coefficients[model,2]+coefficients[model,3]*sex + coefficients[model,4]*age
  
  if (model == "glm") {
    #model="glm"
    #Get linear predictor
    
    if(dist=="NO") {
      t1=rNO(n, mu=(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=rNO(n, mu=(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="GA") {
      t1=rGA(n, mu=exp(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=rGA(n, mu=exp(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="LO") {
      t1=rBI(n, mu=logit_inv(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=rBI(n, mu=logit_inv(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="PO") {
      t1=rNBI(n, mu=exp(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=rNBI(n, mu=exp(lp_2), sigma=sqrt(sigmas[model,2]))
    }
    
  } else if (model == "gee") {
    library(VineCopula)
    #model="gee"
    #Generate basic correlation structure based on correlation parameter
    c0=BiCopSim(n,family=1,par=correlations[model,1])
    
    if(dist=="NO") {
      t1=qNO(c0[,1], mu=(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=qNO(c0[,2], mu=(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="GA") {
      t1=qGA(c0[,1], mu=exp(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=qGA(c0[,2], mu=exp(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="LO") {
      t1=qBI(c0[,1], mu=logit_inv(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=qBI(c0[,2], mu=logit_inv(lp_2), sigma=sqrt(sigmas[model,2]))
    } else if (dist=="PO") {
      t1=qNBI(c0[,1], mu=exp(lp_1), sigma=sqrt(sigmas[model,1]))
      t2=qNBI(c0[,2], mu=exp(lp_2), sigma=sqrt(sigmas[model,2]))
    }
    
  } else if (model == "re_nosig") {
    #model="re_nosig"
    
    #Generate random effects
    b_var=correlations[model,1]
    re_val=rnorm(n,mean=0,sd=((sqrt(b_var))))
    
    if(dist=="NO") {
      t1=rNO(n, mu=(lp_1+re_val), sigma=exp(sigmas[model,1]))
      t2=rNO(n, mu=(lp_2+re_val), sigma=exp(sigmas[model,2]))
    } else if (dist=="GA") {
      t1=rGA(n, mu=exp(lp_1+re_val), sigma=exp(sigmas[model,1]))
      t2=rGA(n, mu=exp(lp_2+re_val), sigma=exp(sigmas[model,2]))
    } else if (dist=="LO") {
      t1=rBI(n, mu=logit_inv(lp_1+re_val), sigma=exp(sigmas[model,1]))#*exp(re_val)####COME BACK FOR THIS
      t2=rBI(n, mu=logit_inv(lp_2+re_val), sigma=exp(sigmas[model,2]))#*exp(re_val)
    } else if (dist=="PO") {
      t1=rNBI(n, mu=exp(lp_1+re_val), sigma=exp(sigmas[model,1]))
      t2=rNBI(n, mu=exp(lp_2+re_val), sigma=exp(sigmas[model,2]))
    }
  } else if (model == "lme4" | model == "gamm") {
    #model="lme4"; model="gamm"
    
    #Generate random effects
    b_var=correlations[model,1]
    re_val=rnorm(n,mean=0,sd=((sqrt(b_var))))
    
    if(dist=="NO") {
      t1=rNO(n, mu=(lp_1+re_val), sigma=(sigmas[model,1]))
      t2=rNO(n, mu=(lp_2+re_val), sigma=(sigmas[model,2]))
    } else if (dist=="GA") {
      t1=rGA(n, mu=exp(lp_1+re_val), sigma=(sigmas[model,1]))
      t2=rGA(n, mu=exp(lp_2+re_val), sigma=(sigmas[model,2]))
    } else if (dist=="LO") {
      t1=rBI(n, mu=logit_inv(lp_1+re_val), sigma=(sigmas[model,1]))#*exp(re_val)####COME BACK FOR THIS
      t2=rBI(n, mu=logit_inv(lp_2+re_val), sigma=(sigmas[model,2]))#*exp(re_val)
    } else if (dist=="PO") {
      t1=rNBI(n, mu=exp(lp_1+re_val), sigma=(sigmas[model,1]))
      t2=rNBI(n, mu=exp(lp_2+re_val), sigma=(sigmas[model,2]))
    }
  } else if (model == "re_np") {
    t1=rep(NA,n)
    t2=rep(NA,n)
  } else if (startsWith(model, "cop")) {
    #model="cop"
    
    library(VineCopula)
    cop_model_names=c("cop","cop_n","cop_j","cop_g","cop_f","cop_amh","cop_fgm","cop_pl","cop_h","cop_t")
    vine_cop_family=c(3,    1,      6,      4,      5,      NA,       NA,       NA,      NA,     NA)
    cop_vine=vine_cop_family[grep(paste("\\b",model,"\\b",sep=""),cop_model_names)]

    if(is.na(cop_vine)) {
      t1=rep(NA,n)
      t2=rep(NA,n)
    } else {
        simCop=BiCopSim(N=n, family=cop_vine,par=correlations[model,])
        
        if(dist=="NO") {
          t1=qNO(simCop[,1], mu=(lp_1), sigma=(sigmas[model,1]))
          t2=qNO(simCop[,2], mu=(lp_2), sigma=(sigmas[model,2]))
        } else if (dist=="GA") {
          t1=qGA(simCop[,1], mu=exp(lp_1), sigma=(sigmas[model,1]))
          t2=qGA(simCop[,2], mu=exp(lp_2), sigma=(sigmas[model,2]))
        } else if (dist=="LO") {
          t1=qBI(simCop[,1], mu=logit_inv(lp_1), sigma=(sigmas[model,1]))
          t2=qBI(simCop[,2], mu=logit_inv(lp_2), sigma=(sigmas[model,2]))
        } else if (dist=="PO") {
          t1=qNBI(simCop[,1], mu=exp(lp_1), sigma=(sigmas[model,1]))
          t2=qNBI(simCop[,2], mu=exp(lp_2), sigma=(sigmas[model,2]))
        }
    }
  }
  return(as.vector(rbind(t1, t2)))
}

evaluateModels <- function(fits,model_list=rownames(fits$correlations),vg_sims=100) {
  
  #Extract model info
  logliks=fits$logliks
  actuals=fits$actuals
  dist=fits$dist
  ses=fits$ses
  y=fits$y
  n=length(y)/2
  
  coefficients=fits$coefficients
  correlations=fits$correlations
  sigmas=fits$sigmas
  model_list_complete=model_list
  
  # Simulate vg_sims times from each model formulation
  
  #vg_sims=100
  #model_list=rownames(fits$correlations);vg_sims=100
  sim_model_out=list()
  for(model in model_list) {
    sim_model_out[[model]]=matrix(NA,ncol=vg_sims,nrow=length(y))
    for(i in 1:vg_sims) {
      sim_model_out[[model]][,i]=sim_model(model,dist,n,coefficients,sigmas,correlations)
    }
    if(all(is.na(sim_model_out[[model]]))) {
      model_list_complete=model_list_complete[model_list_complete!=model]
    }
  }
  
  # Calculate Variogram scores
  print("CALCULATING VARIOGRAM SCORES (2. CALCULATION)")
  library(scoringRules)
  
  w_vs=matrix(1,ncol=length(y),nrow=length(y))
  w_vs_0=matrix(0,ncol=length(y),nrow=length(y))
  # Set values where row = col + n or col = row + n to 10
  # For every pair of adjacent observations
  for (i in seq(1, n*2, by=2)) {
    w_vs[i, i+1] <- ((n^2-2*(n-1))/(2*(n-1)))^2
    w_vs[i+1, i] <- ((n^2-2*(n-1))/(2*(n-1)))^2
    w_vs_0[i, i+1] <- 1
    w_vs_0[i+1, i] <- 1
  }
  
  vs2_wt=vs2=es=vs2_wt_coronly=vs1=rep(NA,length(model_list_complete))
  names(vs2_wt)=names(vs2)=names(es)=names(vs1)=names(vs2_wt_coronly)=model_list_complete
  
  for (model in model_list_complete) {
    vs2_wt[model]=        vs_sample(y=y,dat=sim_model_out[[model]],p=2,w_vs=w_vs)
    vs2_wt_coronly[model]=vs_sample(y=y,dat=sim_model_out[[model]],p=2,w_vs=w_vs_0)
    vs2[model]=           vs_sample(y=y,dat=sim_model_out[[model]],p=2)
    vs1[model]=           vs_sample(y=y,dat=sim_model_out[[model]],p=1)
    es[model]=            es_sample(y=y,dat=sim_model_out[[model]])
  }
  
  vs2
  
  return(list(
    coefficients=coefficients
    , ses=ses
    , sigmas=sigmas
    , correlations=correlations
    , logliks=logliks
    , vs2_wt=vs2_wt
    , vs2=vs2
    , es=es
    , vs1=vs1
    ,vs2_wt_coronly=vs2_wt_coronly
  ))
}
# 
# evaluateModels <- function(fits,model_list=rownames(fits$correlations),vg_sims=100) {
#   
#   #model_list=rownames(fits$correlations);vg_sims=100
#   print("EXTRACTING & TRANSFORMING MODEL INFORMATION")
#   
#   logliks=fits$logliks
#   actuals=fits$actuals
#   dist=fits$dist
#   y=fits$y
#   n=length(y)/2
#   
#   
#   ##These are updated in this function 
#   coefficients_in=fits$coefficients
#   ses_in=fits$ses
#   sigmas_in=fits$sigmas
#   correlations_in=fits$correlations
#   
#   sigmas=sigmas_in*0
#   correlations=correlations_in*0
#   
#   ###### FOR RANDOM EFFECT MODELS WE HAVE TO ESTIMATE MARGINAL PARAMETERS HERE
#   #"re_nosig","re_np","lme4","gamm"
#   #This is because the coefficients are not the marginal parameters
#   
#   
#   #####Normal
#   
#   if(dist=="NO") {
#     coefficients=coefficients_in
#     ses=ses_in
#     for(model in model_list) {
#       if(model == "glm" | model == "gee") {
#         sigmas[model,]=sqrt(sigmas_in[model,])
#         correlations[model,]=correlations_in[model,]
#       } else if (model == "re_nosig") {
#         sigmas[model,]=sqrt(exp(sigmas_in[model,])^2+correlations_in[model,])
#         correlations[model,]=correlations_in[model,]/(sigmas[model,1]*sigmas[model,2])
#       } else if (model=="re_np") {
#         sigmas[model,]=exp(sigmas_in[model,])
#         correlations[model,]=correlations_in[model,]/(sigmas[model,1]*sigmas[model,2])
#       } else if (model == "gamm" | model == "lme4") {
#         sigmas[model,]=sqrt((sigmas_in[model,])^2+correlations_in[model,])
#         correlations[model,]=correlations_in[model,]/(sigmas[model,1]*sigmas[model,2])
#       }
#       else if (startsWith(model,"cop")) {
#         sigmas[model,]=(sigmas_in[model,])
#         if(model=="cop_n") {
#           correlations[model,]=correlations_in[model,]  
#         } else {
#           correlations[model,]=simCopulaCorrelation(model=model,par=c(coefficients[model,],sigmas[model,],correlations_in[model,]),dist=dist)
#         }
#         
#       }
#       
#       mu1_all=coefficients[,1]
#       mu2_all=coefficients[,2]
#       x1_all=coefficients[,3]
#       x2_all=coefficients[,4]
#       
#       a_all=sigmas[,1]
#       b_all=sigmas[,2]
#       c_all=correlations[,1]
#       
#     }
#   } else if (dist=="GA") {
# 
#     #ake in fit information and setup matrixes
#     a_vector=correlations*0
#     b_vector=correlations*0
#     
#     mu_link=log
#     sigma_link=log
#     mu_linkinv=exp
#     sigma_linkinv=exp
#     
#     ses=ses_in
#     
#     coefficients=coefficients_in*0
#     coefficients[,c(1,2)]=mu_linkinv(coefficients_in[,c(1,2)])
#     coefficients[,c(3,4)]=coefficients_in[,c(3,4)]
#     
#     for(model in model_list) {
#       if(model == "glm" | model == "gee") {
#         sigmas[model,]=sqrt(sigmas_in[model,])
#         correlations[model,]=correlations_in[model,]
#       } else if (model == "re_nosig") {
#         
#         b_var_in=correlations_in[model,]
#         sigma1_in=exp(sigmas_in[model,1])
#         sigma2_in=exp(sigmas_in[model,2])
#         
#         mu1_in=exp(coefficients_in[model,1])
#         mu2_in=exp(coefficients_in[model,2])
#         x1_in=coefficients[model,3]
#         x2_in=coefficients[model,4]
#         
#         re_sim_out=simRE(n,dist,mu1=mu1_in,mu2=mu2_in,x1=x1_in, x2=x2_in,sigma1=sigma1_in,sigma2=sigma2_in,b_var=b_var_in,sims=100)
# 
#         correlations[model,]=re_sim_out["corr"]
#         sigmas[model,]=sqrt(sqrt(re_sim_out[c("var_1","var_2")])/c(re_sim_out["mean_1"],re_sim_out["mean_2"]))
#         
#         coefficients[model,c(1,2)]=c(re_sim_out["mean_1"],re_sim_out["mean_2"])
#         
#       } else if (model=="re_np") {
#         
#         b_var_in=correlations_in[model,]
#         sigma1_in=exp(sigmas_in[model,1])
#         sigma2_in=exp(sigmas_in[model,2])
#         
#         mu1_in=exp(coefficients_in[model,1])
#         mu2_in=exp(coefficients_in[model,2])
#         x1_in=coefficients[model,3]
#         x2_in=coefficients[model,4]
#         
#         re_sim_out=simRE(n,dist,mu1=mu1_in,mu2=mu2_in,x1=x1_in, x2=x2_in,sigma1=sigma1_in,sigma2=sigma2_in,b_var=b_var_in,sims=100)
#         
#         correlations[model,]=re_sim_out["corr"]
#         sigmas[model,]=sqrt(sqrt(re_sim_out[c("var_1","var_2")])/c(re_sim_out["mean_1"],re_sim_out["mean_2"]))
#         
#         coefficients[model,c(1,2)]=c(re_sim_out["mean_1"],re_sim_out["mean_2"])
#         
#         
#       } else if (model == "gamm" | model == "lme4") {
#         
#         b_var_in=correlations_in[model,]
#         sigma1_in=(sigmas_in[model,1])
#         sigma2_in=(sigmas_in[model,2])
#         
#         mu1_in=exp(coefficients_in[model,1])
#         mu2_in=exp(coefficients_in[model,2])
#         x1_in=coefficients[model,3]
#         x2_in=coefficients[model,4]
#         
#         re_sim_out=simRE(n,dist,mu1=mu1_in,mu2=mu2_in,x1=x1_in, x2=x2_in,sigma1=sigma1_in,sigma2=sigma2_in,b_var=b_var_in,sims=100)
#         
#         correlations[model,]=re_sim_out["corr"]
#         sigmas[model,]=sqrt(sqrt((re_sim_out[c("var_1","var_2")])/(c(re_sim_out["mean_1"],re_sim_out["mean_2"])^2)))
#         sigmas[model,]
#         
#         coefficients[model,c(1,2)]=c(re_sim_out["mean_1"],re_sim_out["mean_2"])
#         
#         #a=1/sigmasq
#         #musq*ssq=var
#         
#       }
#       else if (startsWith(model,"cop")) {
#         sigmas[model,]=(sigmas_in[model,])
#         if(model=="cop_n") {
#           correlations[model,]=correlations_in[model,]  
#         } else {
#           correlations[model,]=simCopulaCorrelation(model=model,par=c(coefficients[model,],sigmas[model,],correlations_in[model,]),dist=dist)
#         }
#         
#       }
#       
# 
#     
#       #mu1=exp(coefficients[model,1])
#       #mu2=exp(coefficients[model,2])
#       #x1=coefficients[model,3]
#       #x2=coefficients[model,4]
#       
#       #Var(Y_t)=sigma_t^2 * mu_t^2 | sigma=1/sqrt(alpha)
#       #Var(Y_t)=(1/sqrt(alpha))^2 * mu_t^2=1/alpha *  | alpha = 1/Var(Y_t)
#       #a=
#       #b=
#       #c=NA
#       
#     }
#     mu1_all=coefficients[,1]
#     mu2_all=coefficients[,2]
#     x1_all=coefficients[,3]
#     x2_all=coefficients[,4]
#     
#     a_all=1/(sigmas[,2]^2)
#     b_all=rep(1,length(a_all))
#     c_all=correlations[,1]
#   }
#   
#   c_all[c_all<0]=0
#   
#   sig_cor_out=cbind(sigmas,correlations)
#   colnames(sig_cor_out)=c("sigma1","sigma2","cor")
#   sig_cor_out
#   
#   # Function to compute b given a and c based on the equation:
#   b = (-sqrt(-4*a_all^2*c_all^2 + a_all^2 - 4*a_all*c_all^2) - 2*a_all*c_all^2 + a_all - 2*c_all^2) / (2*c_all^2)
# 
#   compute_b <- function(b, a, c) {
#     # Function whose root you want to find
#     value = ((sqrt(a * b) / (a + b + 1)) - c)
#     return(value)
#   }
#   
#   for (i in c(2,3,5,6,7:11)) { # Skip 1 (glm) and 4 (re_np)
#     # Objective function: squared value to find root
#     objective <- function(b) {
#       compute_b(b, a_all[i], c_all[i])^2
#     }
#     opt <- optim(par = b_all[i], fn = objective, method = "BFGS")
#     b_all[i] <- opt$par
#   } 
#   for (i in c(1,4,12:16)) {
#     b_all[i]=0
#   }
#   
#   names(b_all)=names(a_all)=names(c_all)=model_list
#   
#   if(any(is.na(correlations))|any(is.na(sigmas))) {
#     print("Some correlations or sigmas could not be calculated");
#     rownames(sig_cor_out)[apply(sig_cor_out, 1, function(x) any(is.na(x)))]
#     #print(sig_cor_out)
#   }
#   
#   
#   model_list_complete=rownames(sig_cor_out)[apply(sig_cor_out, 1, function(x) !any(is.na(x)))]
#   
#   ##########EXTRACT PARAMETERS FROM INPUT ##########
#   print("CALCULATING VARIOGRAM SCORES (1. SIMULATION)")
#   vg_x=vg_x_mean=list(); vg=list(); datasets=list()
#   for(model in model_list_complete) {
#     print(model)
#     for (run in 1:vg_sims) {
#       a=a_all[model];b=b_all[model];c=c_all[model];
#       mu1=mu1_all[model];mu2=mu2_all[model];
#       x1=x1_all[model];x2=x2_all[model];
#       datasets[[model]][[run]]=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=fits$dist,x1=x1,x2=x2)$random_variable
#     }
#   }
#   
#   print("CALCULATING VARIOGRAM SCORES (2. CALCULATION)")
#   library(scoringRules)
#   
#   w_vs=matrix(1,ncol=length(y),nrow=length(y))
#   w_vs_0=matrix(0,ncol=length(y),nrow=length(y))
#   # Set values where row = col + n or col = row + n to 10
#   # For every pair of adjacent observations
#   for (i in seq(1, n*2, by=2)) {
#     w_vs[i, i+1] <- ((n^2-2*(n-1))/(2*(n-1)))^2
#     w_vs[i+1, i] <- ((n^2-2*(n-1))/(2*(n-1)))^2
#     w_vs_0[i, i+1] <- 1
#     w_vs_0[i+1, i] <- 1
#   }
#   
#   vs2_wt=vs2=es=vs2_wt_coronly=vs1=rep(NA,length(model_list_complete))
#   names(vs2_wt)=names(vs2)=names(es)=names(vs1)=names(vs2_wt_coronly)=model_list_complete
#   
#   for (model in model_list_complete) {
#     vs2_wt[model]=vs_sample(y=y,dat=as.matrix(do.call(cbind, datasets[[model]])),p=2,w_vs=w_vs)
#     vs2_wt_coronly[model]=vs_sample(y=y,dat=as.matrix(do.call(cbind, datasets[[model]])),p=2,w_vs=w_vs_0)
#     vs2[model]=vs_sample(y=y,dat=as.matrix(do.call(cbind, datasets[[model]])),p=2)
#     vs1[model]=vs_sample(y=y,dat=as.matrix(do.call(cbind, datasets[[model]])),p=1)
#     es[model]=es_sample(y=y,dat=as.matrix(do.call(cbind, datasets[[model]])))
#   }
#   
#   return(list(
#     coefficients=coefficients
#     , ses=ses
#     , sigmas=sigmas
#     , correlations=correlations
#     , logliks=logliks
#     , vs2_wt=vs2_wt
#     , vs2=vs2
#     , es=es
#     , vs1=vs1
#     ,vs2_wt_coronly=vs2_wt_coronly
#   ))
#   
# }
# 
# 

simCopulaCorrelation <-function(model,par,n=1000,sims=100,dist) {
  library(VineCopula)
  #par=c(coefficients[model,],sigmas[model,],correlations_in[model,])
  print(model)
  cop_model_names=c("cop","cop_n","cop_j","cop_g","cop_f","cop_amh","cop_fgm","cop_pl","cop_h","cop_t")
  vine_cop_family=c(3,    1,      6,      4,      5,      NA,       NA,       NA,      NA,     NA)
  cop_vine=vine_cop_family[grep(paste("\\b",model,"\\b",sep=""),cop_model_names)]
  
  print(cop_vine)
  
  if(is.na(cop_vine)) {
    return(NA)
  } else {
    corr_mat=matrix(NA,nrow=sims,ncol=1)
    for (i in 1:sims) {
      simCop=BiCopSim(N=n, family=cop_vine,par=par[length(par)])
      qFUN=if(dist=="NO") {qNO} else if (dist=="GA") {qGA} else if (dist=="PO" | dist=="NB") {qNBI}
      t1 = qFUN(simCop[,1],mu=par["as.factor(time == 1)FALSE"],sigma=par[length(par)-2])
      t2 = qFUN(simCop[,2],mu=par["as.factor(time == 1)TRUE"],sigma=par[length(par)-1])
      
      #Would this work without parameters...? For the normal maybe, but otherwise I don't think it will be scale invariant
      #t1 = qFUN(simCop[,1])
      #t2 = qFUN(simCop[,2])
      
      corr_mat[i,1]=cor(t1,t2)
    }
    return(c(t1,t2))
  }
}

plotDist <- function (dataset,dist) {
  
  require(gamlss)  
  require(latex2exp)
  require(ggplot2)
  require(ggpubr)
  
  num_margins=length(unique(dataset[,"time"]))
  
  margin_data=list()
  margin_unif=list()
  margin_fit=list()
  
  for (i in 1:num_margins) {
    margin_data[[i]]<-(dataset[dataset[,"time"]==i-1,"random_variable"])
    margin_fit[[i]]<-gamlss(margin_data[[i]]~1,family=dist)
    margin_unif[[i]]<-pnorm(margin_fit[[i]]$residuals)
  }
  
  ##plot.new()
  #par(mfrow=c(1,num_margins))
  
  #for (i in 1:num_margins) {histDist(margin_data[[i]],family=dist,xlab=TeX(paste("$Y_",i,"$")),main=paste("Histogram of margin",i,"and fitted",dist))}
  #invisible(readline(prompt="Press [enter] to continue"))
  
  plots=list()
  
  z=1
  for (i in 1:(num_margins)) {
      for (j in 1:(num_margins)) {
        if(i==j) {
          input_data=data.frame(margin_data[[i]])
          colnames(input_data)<-"X1"
          p <- ggplot(input_data, aes(x=X1)) + 
            geom_histogram() + 
            labs(x = TeX(paste("$Y_",i,"$")))
        }
        if(i!=j) {
          input_data=data.frame(cbind(margin_unif[[i]],margin_unif[[j]]))
          p=ggplot(data=input_data,aes(x=X1,y=X2)) +
            #geom_point(size=0.25,color="black") + 
            geom_density_2d(contour_var="density",bins=20,color="black") + 
            scale_fill_brewer() +
            labs(x = TeX(paste("$Y_",i,"$")), y=TeX(paste("$Y_",j,"$")),fill="density")
        }
        
        plots[[z]]=p
        z=z+1
      }
  }
  ggarrange(plotlist=plots,ncol=num_margins,nrow=num_margins)
  
}

generateMvtDist<-function(dist,mu_vector,sigma_vector,rho_vector) {
  
  if(dist=="NO") {
    require(MASS)
    cor_matrix<-diag(rep(1,length(sigma_vector)))
    if (length(rho_vector)==1) {
      cor_matrix[lower.tri(cor_matrix,diag=FALSE)]<-rho_vector[1]
      cor_matrix[upper.tri(cor_matrix,diag=FALSE)]<-rho_vector[1]
      
      cov_matrix<-diag(sigma_vector) %*% cor_matrix %*% diag(sigma_vector)
    } else {
      cor_matrix=cor_matrix+rho_vector
      cov_matrix<-diag(sigma_vector) %*% cor_matrix %*% diag(sigma_vector)
    }
    
    data<-mvrnorm(n,mu=mu_vector,Sigma = cov_matrix)
  }
  
  data_output<-cbind(data[,1],0)
  for (i in 2:ncol(data)) {
    data_output<-rbind(data_output,cbind(data[,i],i-1)) 
  }
  colnames(data_output) <- cbind("random_variable","time")
  return(data_output)
}

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
  return_list=list(colMeans(optim_cov_outputs),sqrt(diag(cov(optim_cov_outputs))*n)/sqrt(n))
  names(return_list)=c("coefficients","ses")
  return(return_list)
}