simulateCorrelatedVarNOGJRM <- function(n,a,b,mu1,mu2)  {
  
  set.seed(100)
  
  #Loading required packages
  require(gamlss)
  require(MASS)
  require(gee)
  require(VGAM)
  
  #Simulating bivariate random variable according to functional input
  w<-rbeta(n,a,b)
  gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1)
  gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2)
  
  #Transforming data to format required for random effect models
  patient<-as.factor(seq(1:n))
  dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0)
                               ,cbind(patient,gamma_c_mu2,1)))
  colnames(dataset)<-c("patient","random_variable","time")
  
  #Running generalised linear model for the realisation
  model_glm <- glm(random_variable~as.factor(time==1)
                   , data=dataset
                   , family=Gamma(link = "log")
                   , maxit=1000)
  
  #Running GLMM with no sigma time variable
  
  model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                           , data=dataset
                           , family=GA()
                           ) 
  

  #Running GLMM with sigma time variable
  model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                     , sigma.formula=~as.factor(time==1)
                     , data=dataset
                     , family=GA()
                     , start.from = model_re_nosig #Optional
                     , method=CG(1000))
  
  
  #Running GEE
  model_gee<-gee(random_variable~as.factor(time==1)
                 , id=patient
                 , data=dataset
                 , family=Gamma(link = "log")
                 , maxiter=25
                 , corstr = "exchangeable")
  
  #Extracting coefficient estimates from each of the models
  summary_glm<-c( summary(model_glm)$coeff[1]
                  ,summary(model_glm)$coeff[2]
                  ,summary(model_glm)$coeff[3]
                  ,summary(model_glm)$coeff[4]
  )
  summary_gee<-c( summary(model_gee)$coeff[1]
                  ,summary(model_gee)$coeff[2]
                  ,summary(model_gee)$coeff[3] 
                  ,summary(model_gee)$coeff[4] 
  )
  
  invisible(capture.output(
    summary_re_nosig<-c( summary(model_re_nosig)[1]
                         ,summary(model_re_nosig)[2]
                         ,summary(model_re_nosig)[4]
                         ,summary(model_re_nosig)[5]
    )
  ))
  invisible(capture.output(
    summary_re<-c( summary(model_re)[1]
                   ,summary(model_re)[2]
                   ,summary(model_re)[5]
                   ,summary(model_re)[6]
    )
  ))
  
  summary_cop<-c( 0,0,0,0 #Blank as this is combined with GJRM function at a later point
  )
  summary_cop_n<-c( 0,0,0,0 #Blank as this is combined with GJRM function at a later point
  )
  
  #Calculating true simulated estimates for the distribution based on the parameters
  actuals<-c( log(a*(1/mu1))
              , -log(a*(1/mu1))+log(a*(1/mu2))
              , 0
              , 0
  )
  
  #Combining simulation estimates into a single table
  output<-rbind(summary_glm
                , summary_gee
                , summary_re_nosig
                , summary_re
                , summary_cop
                , summary_cop_n
                , actuals)
  
  colnames(output)<-c("Time 1 Intercept","Time 2 Intercept","Time 1 SE","Time 2 SE")
  
  return(output)
}

simulateCorrelatedVarGJRM <- function(n,a,b,mu1,mu2)  {
  
  set.seed(100)
  
  #Loading required packages
  require(GJRM)
  require(MASS)
  
  #Simulating bivariate random variable according to functional input
  w<-rbeta(n,a,b)
  gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1)
  gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2)
  
  ####TEMP
  #u<-pgamma(gamma_c_mu1,shape=a,scale=1/mu1)
  #v<-pgamma(gamma_c_mu2,shape=a,scale=1/mu2)
  #summary(BiCopSelect(u,v,familyset = c(0,1,2,3,4,5,6,13,16)))
  ###TEMP
  
  eq.mu.1 <- gamma_c_mu1~1
  eq.mu.2 <- gamma_c_mu2~1
  fl <- list(eq.mu.1, eq.mu.2)
  model_copula <- gjrm(fl
                       , margins = c("GA" , "GA") 
                       , copula = "C0"
                       , data=data.frame(gamma_c_mu1,gamma_c_mu2)
                       , model="B")
  model_copula_n <- gjrm(fl
                         , margins = c("GA" , "GA") 
                         , copula = "N"
                         , data=data.frame(gamma_c_mu1,gamma_c_mu2)
                         , model="B")
  
  
  #Extracting coefficient estimates from each of the models
  #First four models are set as blank and run in a separate function
  summary_glm<-c( 0,0,0,0
  )
  summary_gee<-c( 0,0,0,0
  )
  summary_re_nosig<-c( 0,0,0,0
  )
  summary_re<-c( 0,0,0,0
  )  
  summary_cop<-c( model_copula$coefficients[1]
                  , model_copula$coefficients[2] - model_copula$coefficients[1]
                  , summary(model_copula)$tableP1[2] #SE for time 0
                  , summary(model_copula)$tableP2[2] #SE for time 1
  )
  summary_cop_n<-c( model_copula_n$coefficients[1]
                    , model_copula_n$coefficients[2] - model_copula_n$coefficients[1]
                    , summary(model_copula_n)$tableP1[2] #SE for time 0
                    , summary(model_copula_n)$tableP2[2] #SE for time 1
  )
  actuals<-c( 0
              , 0
              , cor(gamma_c_mu1,gamma_c_mu2,method="kendall") 
              , 0
  )
  
  output<-rbind(summary_glm
                , summary_gee
                , summary_re_nosig
                , summary_re
                , summary_cop
                , summary_cop_n
                , actuals)
  
  colnames(output)<-c("Time 1 Intercept","Time 2 Intercept","Time 1 SE","Time 2 SE")
  
  return(output)
}

##############Run simulations for non-GJRM models
results<-list()
#a=.1+.1*1:20; b=.1+.1*1:20; mu1=1; mu2=2; n=100
a=.1+.1*1:20; b=.1+.1*1:20; mu1=1; mu2=2; n=1000

#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1;
start=Sys.time()
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    for (k in 1:length(mu1)) {
      for (l in 1:length(mu2)) {
        set.seed(z)
        results[[z]] <- rbind(tryCatch({
          simulateCorrelatedVarNOGJRM(n,a[i],b[j],mu1[k],mu2[l])}
          , finally={simulateCorrelatedVarNOGJRM(n,a[i],b[j],mu1[k],mu2[l])})
          , c(a[i],b[j],mu1[k],mu2[l]))
        print(c(z,length(a)*length(b)*length(mu1)*length(mu2)
                ,z/(length(a)*length(b)*length(mu1)*length(mu2)) 
                , (Sys.time()-start)
                , (Sys.time()-start) / (z/(length(a)*length(b)*length(mu1)*length(mu2))))
        )
        z = z + 1
      }
    }
  }
}

results_re <- results
save(results_re,file="results_NOGJRM_n1000_geefix.rds")

##############Run simulations for GJRM models
results<-list()
#a=.1+.1*1:20; b=.1+.1*1:20; mu1=1; mu2=2; n=100
a=.1+.1*1:20; b=.1+.1*1:20; mu1=1; mu2=2; n=1000

#Code to iterate through various shapes of the bivariate distribution and fit the GJRM models
i=1; j=1; k=1; l=1; z=1;
start=Sys.time()
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    for (k in 1:length(mu1)) {
      for (l in 1:length(mu2)) {
        set.seed(z)
        results[[z]] <- rbind(tryCatch({
          simulateCorrelatedVarGJRM(n,a[i],b[j],mu1[k],mu2[l])}
          , finally={simulateCorrelatedVarGJRM(n,a[i],b[j],mu1[k],mu2[l])})
          , c(a[i],b[j],mu1[k],mu2[l]))
        print(c(z,length(a)*length(b)*length(mu1)*length(mu2)
                , z/(length(a)*length(b)*length(mu1)*length(mu2)) 
                , (Sys.time()-start)
                , (Sys.time()-start) / (z/(length(a)*length(b)*length(mu1)*length(mu2))))
        )
        z = z + 1
      }
    }
  }
}
results_gjrm<-results
save(results_gjrm,file="results_GJRM_C0_N_n1000.rds")


#load(file="results_GJRM.rds")
#load(file="results_GJRM_AMH.rds")

#######Combine GJRM and non-GJRM results into one table for comparison
load(file="results_GJRM_C0_N_n1000.rds")
load(file="results_NOGJRM_n1000_geefix.rds")

results=results_gjrm
for(i in 1:length(results_gjrm)) {
  results[[i]][1:7,]=results_gjrm[[i]][1:7,]+results_re[[i]][1:7,]
}

save(results,file="results_combined_N_C0_n1000_geefix.rds")
