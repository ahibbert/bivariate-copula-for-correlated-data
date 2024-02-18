generateBivDist <- function(n,a,b,c,mu1,mu2,dist) {
  
  set.seed(100)
  
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
  
  if(dist=="PO") {
    
    #Compound multiple poisson of Stein & Juritz, 1987
    c<-rgamma(n,c)
  
    margin_1=vector(length = n) 
    margin_2=vector(length = n) 
    for (i in 1:n) {
      margin_1[i]=rpois(1,mu1*c[i])
    }
    for (i in 1:n) {
      margin_2[i]=rpois(1,mu2*c[i])
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

fitBivModels <-function(data,dist,include="ALL") {
  
  #Data Setup
  gamma_c_mu1<-dataset[dataset$time==0,]
  gamma_c_mu2<-dataset[dataset$time==1,]
  
  if(dist=="GA"){
    actuals<-c( log(a/mu1)
                , log(a/mu2)
                , 0#model_copula$tau
                , 0
                , NA
                ,NA
                ,NA
    )
  }
  if(dist=="NO"){
    actuals<-c( 
      mu1
      , mu2
      , (a*sqrt(1-c^2))/sqrt(n)
      , (b*sqrt(1-c^2))/sqrt(n)
      , sqrt(a^2+b^2-2*a*b*c)/sqrt(n)
      ,NA
      ,NA 
    )
  }
  
  if(include=="ALL" || include=="non-GJRM only" ) {
  
    require(gamlss)
    require(gee)
    require(lme4)
    library(MASS)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    
    if(dist=="GA") {
      model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=Gamma(link = "log"), maxit=1000)
      model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")
      model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA()) 
      model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=GA(), method=CG(1000))
      model_lme4 <- glmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
    }
    if(dist=="NO") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO())))
      invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      model_lme4 <- lmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)
    }
    
    ###Capturing coefficient values and errors from each model
    summary_glm<-c( summary(model_glm)$coeff[1]
                    ,summary(model_glm)$coeff[2]
                    ,summary(model_glm)$coeff[3]
                    ,summary(model_glm)$coeff[4]
                    , AIC(model_glm)
                    , BIC(model_glm)
                    , 2
    )
    summary_gee<-c( summary(model_gee)$coeff[1]
                    , summary(model_gee)$coeff[2]
                    , summary(model_gee)$coeff[3] 
                    , summary(model_gee)$coeff[4] 
                    , NA
                    , NA
                    , 3
    )
    
    
    invisible(capture.output(
      summary_re_nosig<-c( summary(model_re_nosig)[1]
                           ,summary(model_re_nosig)[2]
                           ,summary(model_re_nosig)[4]
                           ,summary(model_re_nosig)[5]
                           , AIC(model_re_nosig)
                           ,BIC(model_re_nosig)
                           , model_re_nosig$df.fit
      )
    ))
    
    invisible(capture.output(
      summary_re<-c( summary(model_re)[1]
                     ,summary(model_re)[2]
                     ,summary(model_re)[5]
                     ,summary(model_re)[6]
                     , AIC(model_re)
                     , BIC(model_re)
                     , model_re$df.fit
      )
    ))
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,AIC(model_lme4)
                      , BIC(model_lme4)
                      ,4)
    
    ###Calculating effective degrees of freedom from Donohue
    X<-getME(model_lme4,name="X")
    Z<-getME(model_lme4,name="Z")
    U<-cbind(X,Z)
    W<-model_lme4@resp$sqrtrwt #weights(model_lme4,type = "working")
    UWU=(t(as.matrix(U))%*%(diag(as.vector(W)))%*%as.matrix(U))
    dim(UWU)
    D<-getME(model_lme4,name="Lambda")
    D_inv<-solve(D)
    dinv_plus_00<-c(0,0,diag(D_inv))
    lme_EDF=sum(diag(UWU%*%solve(UWU+diag(dinv_plus_00))))
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,-2*logLik(model_lme4)+2*lme_EDF
                      , BIC(model_lme4)
                      ,lme_EDF)
    
  }
  
  if(include=="ALL" || include=="GJRM only" ) {
  
    require(GJRM)
    
    #Setting up GJRM equations
    eq.mu.1 <- formula(random_variable~1)
    eq.mu.2 <- formula(random_variable.1~1)
    fl <- list(eq.mu.1, eq.mu.2)
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    
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
                    , 2*5-2*logLik(model_copula)
                    ,BIC(model_copula)
                    , 5
    )
    
    summary_cop_n<-c( model_copula_n$coefficients[1]
                      , model_copula_n$coefficients[2] 
                      , summary(model_copula_n)$tableP1[2] #SE for time 0
                      , summary(model_copula_n)$tableP2[2] #SE for time 1
                      , 2*5-2*logLik(model_copula_n)
                      ,BIC(model_copula_n)
                      , 5
                      
    )
    summary_cop_j<-c( model_copula_j$coefficients[1]
                      , model_copula_j$coefficients[2]
                      , summary(model_copula_j)$tableP1[2] #SE for time 0
                      , summary(model_copula_j)$tableP2[2] #SE for time 1
                      , 2*5-2*logLik(model_copula_j)
                      ,BIC(model_copula_j)
                      , 5
                      
    )
    summary_cop_g<-c( model_copula_g$coefficients[1]
                      , model_copula_g$coefficients[2] 
                      , summary(model_copula_g)$tableP1[2] #SE for time 0
                      , summary(model_copula_g)$tableP2[2] #SE for time 1
                      , 2*5-2*logLik(model_copula_g)
                      ,BIC(model_copula_g)
                      , 5
                      
    )
    summary_cop_f<-c( model_copula_f$coefficients[1]
                      , model_copula_f$coefficients[2]
                      , summary(model_copula_f)$tableP1[2] #SE for time 0
                      , summary(model_copula_f)$tableP2[2] #SE for time 1
                      , 2*5-2*logLik(model_copula_f)
                      ,BIC(model_copula_f)
                      , 5
                      
    )
    summary_cop_amh<-c( model_copula_amh$coefficients[1]
                        , model_copula_amh$coefficients[2]
                        , summary(model_copula_amh)$tableP1[2] #SE for time 0
                        , summary(model_copula_amh)$tableP2[2] #SE for time 1
                        , 2*5-2*logLik(model_copula_amh)
                        ,BIC(model_copula_amh)
                        , 5
                        
    )
    summary_cop_fgm<-c( model_copula_fgm$coefficients[1]
                        , model_copula_fgm$coefficients[2]
                        , summary(model_copula_fgm)$tableP1[2] #SE for time 0
                        , summary(model_copula_fgm)$tableP2[2] #SE for time 1
                        , 2*5-2*logLik(model_copula_fgm)
                        ,BIC(model_copula_fgm)
                        , 5
                        
    )
    summary_cop_pl<-c( model_copula_pl$coefficients[1]
                       , model_copula_pl$coefficients[2]
                       , summary(model_copula_pl)$tableP1[2] #SE for time 0
                       , summary(model_copula_pl)$tableP2[2] #SE for time 1
                       , 2*5-2*logLik(model_copula_pl)
                       ,BIC(model_copula_pl)
                       , 5
                       
    )
    summary_cop_h<-c( model_copula_h$coefficients[1]
                      , model_copula_h$coefficients[2]
                      , summary(model_copula_h)$tableP1[2] #SE for time 0
                      , summary(model_copula_h)$tableP2[2] #SE for time 1
                      , 2*5-2*logLik(model_copula_h)
                      ,BIC(model_copula_h)
                      , 5
                      
    )
    summary_cop_t<-c( model_copula_t$coefficients[1]
                      , model_copula_t$coefficients[2]
                      , summary(model_copula_t)$tableP1[2] #SE for time 0
                      , summary(model_copula_t)$tableP2[2] #SE for time 1
                      , AIC(model_copula_t)
                      ,BIC(model_copula_t)
                      , 5
                      
    )
  }

  ########### 4. Combining results #########
  
  if(include=="ALL") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re,summary_lme4,summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)  
  }
  if(include=="GJRM only") {
    results_exbias<- rbind(summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)
  }
  if(include=="non-GJRM only") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re,summary_lme4,actuals)
  }
  
  return(results_exbias)
  
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
