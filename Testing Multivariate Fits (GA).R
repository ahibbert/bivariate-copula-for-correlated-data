library(VineCopula); library(MASS); library(gamlss); source("common_functions.R")

vs2_all=matrix(0,nrow=num_outer_sims,ncol=6)
simFit=function(n,rho=0,sims=100,data_in=NULL,mean_in=0,sigma_in=1,x1_in=0,x2_in=0,coef_in,dist_name=NA){
    coef_out=matrix(0,nrow=sims,ncol=3)
    dispersion_out=rep(0,sims)
    x1=x1_in
    x2=x2_in

    for(i in 1:sims) {
        simvine=simRVine(n=n,rho=rho,data_in=data_in)
        out=simvine[[1]]  
        loglik_vine=simvine[[2]]
        if(dist_name=="GA") {
          out_adj=qGA(out,mu=exp(mean_in+matrix(rep(x1*coef_in[1],d),ncol=d)+matrix(rep(x2*coef_in[2],d),ncol=d)),sigma=(sigma_in))
        } else if (dist_name =="NB") {
          out_adj=matrix(qNBI(out
                ,mu=exp(mean_in+matrix(rep(x1*coef_in[1],d),ncol=d)+matrix(rep(x2*coef_in[2],d),ncol=d))
                ,sigma=matrix(rep((sigma_in),n*d),ncol=d))
                ,ncol=d)
        } else if (dist_name=="LO") {
          out_adj=matrix(qBI(out
                ,mu=logit_inv((mean_in)+matrix(rep(x1*coef_in[1],d),ncol=d)+matrix(rep(x2*coef_in[2],d),ncol=d)))
                ,ncol=d)
        } else {
           stop("Distribution not recognized")
        }
        

        random_variable=c(out_adj[,1],out_adj[,2],out_adj[,3],out_adj[,4],out_adj[,5])
        patient=as.factor(rep(1:n,d))
        x1_long=rep(x1,d)
        x2_long=as.factor(rep(x2,d))
        data_long=data.frame(random_variable,patient,x1_long,x2_long)

        data_long=data_long[order(data_long$patient),]

        glm_fit= if(dist_name=="GA") {glm(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=Gamma(link="log"))}
        else if (dist_name =="NB") {glm.nb(random_variable~1+x1_long+as.factor(x2_long), data=data_long,maxit=1000)}
        else if (dist_name =="LO") {glm(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=binomial(link="logit"))}
        else {stop("Distribution not recognized")}

        coef_out[i,]=glm_fit$coefficients
        dispersion_out[i]=summary(glm_fit)$dispersion
        #capture GLM outputs into a matrix
    }
    result=cbind(colMeans(coef_out),apply(coef_out,2,sd))
    colnames(result)=c("Mean","SD")
    rownames(result)=c("Intercept","X1","X2")
    return(list(result,loglik_vine,par=simvine[[3]],glm_dispersion=mean(dispersion_out)))
}

set.seed(1000)
library(glmtoolbox)
library(gamlss)
library(mgcv)
library(lme4)
library(mvtnorm)

#Estimate true values
print("Estimating true parameter values...")
x1=runif(n); x2=as.numeric(runif(n)>0.5)
true_sim=simFit(n=1000,rho=rho,sims=true_sims,mean_in=sim_mean,sigma_in=(sim_sigma),x1_in=x1,x2_in=x2,coef_in=c(1,1),dist_name=dist_name)[[1]]

######################################### GAUSSIAN

# helper to safely fit a model and record failures without stopping the simulation loop
failure_log <- data.frame(sim=integer(), model=character(), error=character(), stringsAsFactors=FALSE)
safe_fit <- function(model_name, sim_id, fit_fn) {
  tryCatch(
    list(ok=TRUE, model=fit_fn()),
    error=function(e) {
      failure_log <<- rbind(failure_log, data.frame(sim=sim_id, model=model_name, error=conditionMessage(e), stringsAsFactors=FALSE))
      list(ok=FALSE, model=NULL)
    }
  )
}

empty_coef_se <- function() {
  out <- matrix(NA_real_, nrow=3, ncol=2)
  rownames(out) <- c("(Intercept)","x1_long","as.factor(x2_long)1")
  colnames(out) <- c("Estimate","Std. Error")
  out
}

sims=list()

coef_count=NULL; ses_count=NULL; coef_sum=NULL; ses_sum=NULL; loglik_sum=NULL; loglik_count=NULL

# Initialize timing tracking
model_times <- list(
  GLM = numeric(num_outer_sims),
  GEE = numeric(num_outer_sims),
  RENOSIG = numeric(num_outer_sims),
  LME4 = numeric(num_outer_sims),
  GAMM = numeric(num_outer_sims),
  VineCopula = numeric(num_outer_sims)
)
conv_rates <- list(
  GLM = numeric(num_outer_sims),
  GEE = numeric(num_outer_sims),
  RENOSIG = numeric(num_outer_sims),
  LME4 = numeric(num_outer_sims),
  GAMM = numeric(num_outer_sims)
)

vs2_all=matrix(0,nrow=num_outer_sims,ncol=6)
vs2_wt_all=matrix(0,nrow=num_outer_sims,ncol=6)


for (j in 1:num_outer_sims) {
  print(paste("Outer Simulation ",j," of ",num_outer_sims,sep=""))
  out=pnorm(mvtnorm::rmvnorm(n=n,sigma=matrix(c(1,rho,rho^2,rho^3,rho^4,
                      rho,1,rho,rho^2,rho^3,
                      rho^2,rho,1,rho,rho^2,
                      rho^3,rho^2,rho,1,rho,
                      rho^4,rho^3,rho^2,rho,1),nrow=5)))
  x1=runif(n)-0.5
  x2=as.numeric(runif(n)>0.5)

  if(dist_name=="GA") {
    out_adj=qGA(out,mu=exp(sim_mean+x1+x2),sigma = rep(sim_sigma,n))
  } else if (dist_name=="NB") {
    out_adj=matrix(qNBI(out,mu=exp(sim_mean+matrix(rep(x1,d),ncol=d)+matrix(rep(x2,d),ncol=d)),sigma=matrix(rep(sim_sigma,n*d),ncol=d)),ncol=d)
  } else if (dist_name=="LO") {
    out_adj=matrix(qBI(out
           ,mu=logit_inv((sim_mean)+matrix(rep(x1,d),ncol=d)+matrix(rep(x2,d),ncol=d)))
           ,ncol=d)
  } else {
    stop("Distribution not recognized")
  }
  library(e1071)
  skew_out=round(skewness(out_adj),2)

  random_variable=c(out_adj[,1],out_adj[,2],out_adj[,3],out_adj[,4],out_adj[,5])
  patient=as.factor(rep(1:n,d))
  x1_long=rep(x1,d)
  x2_long=as.factor(rep(x2,d))
  data_long=data.frame(random_variable,patient,x1_long,x2_long)

  data_long=data_long[order(data_long$patient),]

  print("Fitting models...")

  # Helper to time and fit models safely
  safe_fit_timed <- function(model_name, sim_id, fit_fn) {
    time_taken <- system.time({
      result <- safe_fit(model_name, sim_id, fit_fn)
    })[3]
    return(list(result = result, time = time_taken))
  }

  if(dist_name=="GA") {
    fit_result_glm <- safe_fit_timed("GLM", j, function() suppressMessages(glm(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=Gamma(link="log"))))
    fit_glm <- fit_result_glm$result
    model_times$GLM[j] <- fit_result_glm$time
    
    invisible(capture.output(fit_result_gee <- safe_fit_timed("GEE", j, function() suppressMessages(glmgee(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=Gamma(link="log"), id=patient, corstr="Unstructured")))))
    fit_gee <- fit_result_gee$result
    model_times$GEE[j] <- fit_result_gee$time
    
    invisible(capture.output(fit_result_renosig <- safe_fit_timed("RE-NOSIG", j, function() suppressMessages(gamlss(formula=random_variable~1+x1_long+as.factor(x2_long)+random(as.factor(patient)), data=as.data.frame(out_adj), family=GA(),method=RS(1000))))))
    fit_re_nosig <- fit_result_renosig$result
    model_times$RENOSIG[j] <- fit_result_renosig$time
    
    invisible(capture.output(fit_result_lme4 <- safe_fit_timed("LME4", j, function() suppressMessages(glmer(formula=random_variable~1+x1_long+as.factor(x2_long)+ (1|patient), data=as.data.frame(out_adj),family = Gamma(link="log"),control=glmerControl(optCtrl = list(maxfun=200000)))))))
    fit_lme4 <- fit_result_lme4$result
    model_times$LME4[j] <- fit_result_lme4$time
    
    invisible(capture.output(fit_result_gamm <- safe_fit_timed("GAMM", j, function() suppressMessages(gamm(formula=random_variable~1+x1_long+as.factor(x2_long), random=list(patient=~1), data=as.data.frame(out_adj), family=Gamma(link="log"))))))
    fit_gamm <- fit_result_gamm$result
    model_times$GAMM[j] <- fit_result_gamm$time
    
    invisible(capture.output(fit_glm_gamlss <- safe_fit("GAMLSS-GLM", j, function() suppressMessages(gamlss(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=GA)))))
  } else if (dist_name=="NB") {
    fit_result_glm <- safe_fit_timed("GLM", j, function() glm.nb(random_variable~1+x1_long+as.factor(x2_long), data=data_long,maxit=1000))
    fit_glm <- fit_result_glm$result
    model_times$GLM[j] <- fit_result_glm$time
    model_glm <- fit_glm$model
    
    invisible(capture.output(fit_result_gee <- safe_fit_timed("GEE", j, function() glmgee(random_variable~1+x1_long+as.factor(x2_long), data=data_long, init.beta=model_glm$coefficients, family=neg.bin(theta=summary(model_glm)$theta),corstr = "exchangeable"))))
    fit_gee <- fit_result_gee$result
    model_times$GEE[j] <- fit_result_gee$time
    
    invisible(capture.output(fit_result_renosig <- safe_fit_timed("RE-NOSIG", j, function() gamlss(formula=random_variable~1+x1_long+as.factor(x2_long)+random(as.factor(patient)), data=as.data.frame(out_adj), family=NBI(),method=RS(1000)))))
    fit_re_nosig <- fit_result_renosig$result
    model_times$RENOSIG[j] <- fit_result_renosig$time
    
    invisible(capture.output(fit_result_lme4 <- safe_fit_timed("LME4", j, function() glmer.nb(formula=random_variable~1+x1_long+as.factor(x2_long)+ (1|patient), data=as.data.frame(out_adj),control=glmerControl(optCtrl = list(maxfun=200000))))))
    fit_lme4 <- fit_result_lme4$result
    model_times$LME4[j] <- fit_result_lme4$time
    
    invisible(capture.output(fit_result_gamm <- safe_fit_timed("GAMM", j, function() gamm(formula=random_variable~1+x1_long+as.factor(x2_long), random=list(patient=~1), data=as.data.frame(out_adj), family=nb(link="log"),control=list(niterEM=1000)))))
    fit_gamm <- fit_result_gamm$result
    model_times$GAMM[j] <- fit_result_gamm$time
    
    invisible(capture.output(fit_glm_gamlss <- safe_fit("GAMLSS-GLM", j, function() gamlss(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=NBI(),method=RS(1000)))))
  } else if (dist_name=="LO") {
    fit_result_glm <- safe_fit_timed("GLM", j, function() glm(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=binomial(link="logit")))
    fit_glm <- fit_result_glm$result
    model_times$GLM[j] <- fit_result_glm$time
    
    invisible(capture.output(fit_result_gee <- safe_fit_timed("GEE", j, function() glmgee(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=binomial(link="logit"), id=patient, corstr="Unstructured"))))
    fit_gee <- fit_result_gee$result
    model_times$GEE[j] <- fit_result_gee$time
    
    invisible(capture.output(fit_result_renosig <- safe_fit_timed("RE-NOSIG", j, function() gamlss(formula=random_variable~1+x1_long+as.factor(x2_long)+random(as.factor(patient)), data=as.data.frame(out_adj), family=BI(),method=RS(1000)))))
    fit_re_nosig <- fit_result_renosig$result
    model_times$RENOSIG[j] <- fit_result_renosig$time
    
    invisible(capture.output(fit_result_lme4 <- safe_fit_timed("LME4", j, function() glmer(formula=random_variable~1+x1_long+as.factor(x2_long)+ (1|patient), data=as.data.frame(out_adj),family = binomial(link="logit"),control=glmerControl(optCtrl = list(maxfun=200000))))))
    fit_lme4 <- fit_result_lme4$result
    model_times$LME4[j] <- fit_result_lme4$time
    
    invisible(capture.output(fit_result_gamm <- safe_fit_timed("GAMM", j, function() gamm(formula=random_variable~1+x1_long+as.factor(x2_long), random=list(patient=~1), data=as.data.frame(out_adj), family=binomial(link="logit")))))
    fit_gamm <- fit_result_gamm$result
    model_times$GAMM[j] <- fit_result_gamm$time
    
    invisible(capture.output(fit_glm_gamlss <- safe_fit("GAMLSS-GLM", j, function() gamlss(random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=BI()))))
  } else {
    stop("Distribution not recognized")
  }

  model_glm <- fit_glm$model
  model_glm_gamlss <- fit_glm_gamlss$model
  model_gee <- fit_gee$model
  model_re_nosig <- fit_re_nosig$model
  model_lme4 <- fit_lme4$model
  model_gamm <- fit_gamm$model
  #model_gamm=NULL

  residuals_matrix <- if (!is.null(model_glm)) {
    cbind(
      model_glm$residuals[as.character(1:n)],
      model_glm$residuals[as.character( (n+1):(2*n))],
      model_glm$residuals[as.character(((2*n)+1):(3*n))],
      model_glm$residuals[as.character(((3*n)+1):(4*n))],
      model_glm$residuals[as.character(((4*n)+1):(5*n))]
    )
  } else {
    matrix(NA_real_, nrow=n, ncol=d)
  }

  if (!is.null(model_glm_gamlss)) {
    coef_in=model_glm_gamlss$mu.coefficients[2:3]
    if (dist_name=="LO") {
      time_vine <- system.time(simvinefit<-simFit(n=1000,rho=0,sims=100,data_in=pnorm(residuals_matrix),mean_in=(model_glm_gamlss$mu.coefficients[1])
      ,sigma_in=NA,x1_in=x1,x2_in=x2,coef_in=coef_in,dist_name= dist_name))[3]
      model_times$VineCopula[j] <- time_vine
    } else if (dist_name=="NB") {
      
      #Using method from Mitskopoulos1 et al (2022) to convert NBI residuals to uniform

      #hist(residuals_matrix, breaks=50, main=paste("Histogram of NBI residuals - Simulation ",j,sep=""), xlab="Residuals")
      pnorm_residuals_in=pnorm(residuals_matrix)
      time_vine <- system.time(simvinefit<-simFit(n=1000,rho=0,sims=100,data_in=pnorm_residuals_in,mean_in=model_glm_gamlss$mu.coefficients[1]
      ,sigma_in=exp(model_glm_gamlss$sigma.coefficients),x1_in=x1,x2_in=x2,coef_in=coef_in,dist_name = dist_name))[3]
      model_times$VineCopula[j] <- time_vine
    
    } else {
      time_vine <- system.time(simvinefit<-simFit(n=1000,rho=0,sims=100,data_in=pnorm(residuals_matrix),mean_in=model_glm_gamlss$mu.coefficients[1]
      ,sigma_in=exp(model_glm_gamlss$sigma.coefficients),x1_in=x1,x2_in=x2,coef_in=coef_in,dist_name = dist_name))[3]
      model_times$VineCopula[j] <- time_vine
    }

    
  } else {
    simvinefit <- empty_coef_se()
    model_times$VineCopula[j] <- 0
  }

    results_table=list()
    invisible(capture.output(results_table[[1]] <- if (!is.null(model_glm)) summary(model_glm)$coeff[,1:2] else empty_coef_se() ))
    invisible(capture.output(results_table[[2]] <- if (!is.null(model_gee)) summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2] else empty_coef_se() ))
    invisible(capture.output(results_table[[3]] <- if (!is.null(model_re_nosig)) cbind(summary(model_re_nosig)[1:3],summary(model_re_nosig)[5:7]) else empty_coef_se() ))
    invisible(capture.output(results_table[[5]] <- if (!is.null(model_lme4)) summary(model_lme4)$coefficients[,c(1,2)] else empty_coef_se() ))
    invisible(capture.output(results_table[[6]] <- if (!is.null(model_gamm)) cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix))) else empty_coef_se() ))

    rownames(results_table[[1]])=rownames(results_table[[2]])=rownames(results_table[[3]])=rownames(results_table[[5]])=rownames(results_table[[6]])=c("(Intercept)","x1_long","as.factor(x2_long)1")
    colnames(results_table[[1]])=colnames(results_table[[2]])=colnames(results_table[[3]])=colnames(results_table[[5]])=colnames(results_table[[6]])=c("Estimate","Std. Error")
    names(results_table)=c("GLM","GEE","GAMLSS","RE-NP","LME4","GAMM")

    print("Calculating effective degrees of freedom for LME4...")

    lme_EDF <- tryCatch({
      ###Calculating effective degrees of freedom from Donohue
      X<-getME(model_lme4,name="X")[,1:2]
      Z<-getME(model_lme4,name="Z")
      U<-cbind(X,Z)
      W<-model_lme4@resp$sqrtrwt #weights(model_lme4,type = "working")
      UWU=(t(as.matrix(U))%*%(diag(as.vector(W)))%*%as.matrix(U))
      dim(UWU)
      D<-getME(model_lme4,name="Lambda")

      if(sum(D)==0) {summary_lme4[length(summary_lme4)]} else {
      D_inv<-solve(D)
      dinv_plus_00<-c(0,0,diag(D_inv))
      sum(diag(UWU%*%solve(UWU+diag(dinv_plus_00))))
      }
    }, error=function(e) {
      warning(paste("Failed to calculate LME4 EDF for simulation", j, ":", conditionMessage(e)))
      NA_real_
    })

    print("Compiling results...")

        dfs=c( (n*d)-df.residual(model_glm)
        , (n*d)-model_gee$df.residual
        , model_re_nosig$df.fit
        , lme_EDF
        , lme_EDF
        , ((n*d)-df.residual(model_glm)) + d*(d-1)
      )
      logLiks=c(
        logLik(model_glm)
        , model_gee$logLik
        , logLik(model_re_nosig)
        , logLik(model_lme4)
        , if(!is.null(model_gamm)) {logLik(model_gamm$lme)} else {NA}
        , logLik(model_glm)+simvinefit[[2]]
      )

      #Creating summary tables
      coefficients_table= rbind(
        results_table[[1]][,1]
        , results_table[[2]][,1]
        , results_table[[3]][,1]
        , results_table[[5]][,1]
        , results_table[[6]][,1]
        , simvinefit[[1]][,1]
      )
      ses_table= rbind(
        results_table[[1]][,2]
        , results_table[[2]][,2]
        , results_table[[3]][,2]
        , results_table[[5]][,2]
        , results_table[[6]][,2]
        , simvinefit[[1]][,2]
      )
      loglik_table= cbind(
        logLiks,
        dfs,-2*logLiks+2*dfs,-2*logLiks+log(n*d)*dfs
      )

      conv_check=c(
        model_glm$converged,
        model_gee$converged,
        model_re_nosig$converged,
        if (!is.null(model_lme4)) !any( grepl("failed to converge", model_lme4@optinfo$conv$lme4$messages) ) else NA,
        if (!is.null(model_gamm)) !any( grepl("converge", warnings(model_gamm))) else NA,
        NA_real_        

      )

      sigmas=rbind(
        rep(summary(model_glm)$dispersion,d)
        , rep(model_gee$phi,d)
        , if(dist_name=="LO") { rep(1,d) } else {rep((model_re_nosig$sigma.coefficients),d)}
        , rep(summary(model_lme4)$sigma,d)
        , rep(model_gamm$lme$sigma,d)
        , rep(simvinefit$glm_dispersion,d)
      )

      #Extract correlations - for random effect models this is the random effect sd
      correlations=list(
        0
        ,(model_gee$corr)
        ,getSmo(model_re_nosig)$sigb
        ,summary(model_lme4)$varcor$patient[1,1]
        ,var(ranef(model_gamm$lme)[[1]])
        , simvinefit$par
      )
      
      # Store convergence rates
      conv_rates$GLM[j] <- as.numeric(if (!is.null(model_glm)) model_glm$converged else NA)
      conv_rates$GEE[j] <- as.numeric(if (!is.null(model_gee)) model_gee$converged else NA)
      conv_rates$RENOSIG[j] <- as.numeric(if (!is.null(model_re_nosig)) model_re_nosig$converged else NA)
      conv_rates$LME4[j] <- as.numeric(if (!is.null(model_lme4)) !any( grepl("failed to converge", model_lme4@optinfo$conv$lme4$messages) ) else NA)
      conv_rates$GAMM[j] <- as.numeric(if (!is.null(model_gamm)) !any( grepl("converge", warnings(model_gamm))) else NA)

  sims[[j]]=coefficients_table
  rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=rownames(sigmas)=c("GLM","GEE","GAMLSS","LME4","GAMM","VineCopula")
  colnames(loglik_table)=c("LogLik","DF","AIC","BIC")
  names(correlations)=c("GLM","GEE","GAMLSS","LME4","GAMM","VineCopula")

  ############## VARIOGRAM SCORE CALCULATION
   library(scoringRules)

  # Simulate from each fitted model
  sim_model_out=list()
  vs2=vs2_wt=rep(NA,length(rownames(sigmas)))
  names(vs2)=names(vs2_wt)=rownames(sigmas)

  print("Calculating variogram scores...")

  vg_sims=100

  w_vs=matrix(1,ncol=n*d,nrow=n*d)
  w_vs_0=matrix(0,ncol=n*d,nrow=n*d)
  # Set values where row = col + n or col = row + n to 10
  # For every pair of adjacent observations
  for (i in 1:nrow(w_vs)) {
    if (i + 1 <= nrow(w_vs)) {
      w_vs[i, i + 1] <- ((n^2 - 2 * (n - 1)) / (2 * (n - 1)))^2
      w_vs[i + 1, i] <- ((n^2 - 2 * (n - 1)) / (2 * (n - 1)))^2
      w_vs_0[i, i + 1] <- 1
      w_vs_0[i + 1, i] <- 1
    }
  }

  for(model_name in rownames(sigmas)) {
    sim_model_out[[model_name]]=matrix(nrow=n*d,ncol=vg_sims)
    for(i in 1:vg_sims) {
      set.seed(100+i+j*vg_sims)
      sim_model_out[[model_name]][,i] <- as.vector(sim_model_mvt(model=model_name, dist=dist_name, n=n
      , coefficients=coefficients_table, sigma=sigmas, correlations=correlations, x1=x1, x2=x2))
     }
    vs2[model_name] = vs_sample(y=as.vector(out_adj),dat=sim_model_out[[model_name]],p=2)
    vs2_wt[model_name] = vs_sample(y=as.vector(out_adj),dat=sim_model_out[[model_name]],p=2,w_vs=w_vs)
  }
  vs2_all[j,]=vs2
  vs2_wt_all[j,]=vs2_wt


  if (is.null(coef_sum)) {
    coef_sum <- ifelse(is.na(coefficients_table), 0, coefficients_table)
    ses_sum <- ifelse(is.na(ses_table), 0, ses_table)
    coef_count <- !is.na(coefficients_table)
    ses_count <- !is.na(ses_table)
    #now do the same thing for columns of the loglik table
    loglik_sum <- ifelse(is.na(loglik_table), 0, loglik_table)
    loglik_count <- !is.na(loglik_table)
    conv_check_sum <- conv_check
  } else {
    coef_sum <- coef_sum + ifelse(is.na(coefficients_table), 0, coefficients_table)
    ses_sum <- ses_sum + ifelse(is.na(ses_table), 0, ses_table)
    coef_count <- coef_count + !is.na(coefficients_table)
    ses_count <- ses_count + !is.na(ses_table)
    loglik_sum <- loglik_sum + ifelse(is.na(loglik_table), 0, loglik_table)
    loglik_count <- loglik_count + !is.na(loglik_table)
    conv_check_sum <- conv_check_sum + conv_check
  }


}

coefficients_table_extended <- coef_sum / coef_count
ses_table_extended <- ses_sum / ses_count
coefficients_table_extended[coef_count == 0] <- NA_real_
ses_table_extended[ses_count == 0] <- NA_real_
loglik_count_table_extended <- loglik_sum / loglik_count
loglik_count_table_extended[loglik_count == 0] <- NA_real_
conv_check_final <- conv_check_sum / num_outer_sims

# Data cleaning for variogram scores
for (i in 1:nrow(vs2_all)) {
  row_mean <- mean(vs2_all[i, ], na.rm = TRUE)
  vs2_all[i, vs2_all[i, ] > 100 * row_mean] <- NA
}

for (i in 1:nrow(vs2_wt_all)) {
  row_mean <- mean(vs2_wt_all[i, ], na.rm = TRUE)
  vs2_wt_all[i, vs2_wt_all[i, ] > 100 * row_mean] <- NA
}

# Print timing summary
cat("\n========== MODEL TIMING SUMMARY ==========\n")
cat("Average execution time (in seconds) and convergence rates across", num_outer_sims, "simulations:\n\n")

timing_summary <- data.frame(
  Model = c("GLM", "GEE", "RE-NOSIG", "LME4", "GAMM", "VineCopula"),
  Avg_Time_Sec = c(
    mean(model_times$GLM, na.rm=TRUE),
    mean(model_times$GEE, na.rm=TRUE),
    mean(model_times$RENOSIG, na.rm=TRUE),
    mean(model_times$LME4, na.rm=TRUE),
    mean(model_times$GAMM, na.rm=TRUE),
    mean(model_times$VineCopula, na.rm=TRUE)
  ),
  Convergence_Rate = c(
    mean(conv_rates$GLM, na.rm=TRUE),
    mean(conv_rates$GEE, na.rm=TRUE),
    mean(conv_rates$RENOSIG, na.rm=TRUE),
    mean(conv_rates$LME4, na.rm=TRUE),
    mean(conv_rates$GAMM, na.rm=TRUE),
    NA_real_
  )
)


if (nrow(failure_log) > 0) {
  print("Model fit failures captured during simulations:")
  print(failure_log)
}

######PLOTTING

library(ggplot2)
library(gridExtra)

# Assuming coefficients_table, ses_table, and final_high_cor already exist in your environment

coefficients_table_extended <- coefficients_table
ses_table_extended <- ses_table

# Prepare data for plotting
model_names <- rownames(coefficients_table_extended)

colnames(coefficients_table_extended) <- c("Intercept", "X1", "X2")

# Create data frames for each coefficient
plot_data_list <- lapply(1:ncol(coefficients_table_extended), function(i) {
  data.frame(
    Model = factor(model_names, levels = model_names),
    Coefficient = coefficients_table_extended[, i],
    SE = ses_table_extended[, i],
    Lower = coefficients_table_extended[, i] - 1.96 * ses_table_extended[, i],
    Upper = coefficients_table_extended[, i] + 1.96 * ses_table_extended[, i]
  )
})

# Create plots for each coefficient
plots <- lapply(1:ncol(coefficients_table_extended), function(i) {
  # Calculate y-axis limits with small margin
  y_min <- max( min(plot_data_list[[i]]$Lower,true_sim[i,1]-1.96*true_sim[i,2]),true_sim[i,1]-5*true_sim[i,2])
  y_max <- min(max(plot_data_list[[i]]$Upper,true_sim[i,1]+1.96*true_sim[i,2]),true_sim[i,1]+5*true_sim[i,2])
  y_range <- y_max - y_min
  y_margin <- y_range * 0.05  # 5% margin
  
  ggplot(plot_data_list[[i]], aes(x = Model, y = Coefficient)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
    labs(
      title = colnames(coefficients_table_extended)[i],
      x = NULL,
      y = "Coefficient Estimate"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    geom_hline(yintercept = true_sim[i,1], color = "red", linetype = "dashed") +
    geom_hline(yintercept = true_sim[i,1] - 1.96 * true_sim[i,2], color = "blue", linetype = "dashed") +
    geom_hline(yintercept = true_sim[i,1] + 1.96 * true_sim[i,2], color = "blue", linetype = "dashed") +
    scale_color_manual(values = c("True Coefficient" = "red", "True 95% CI" = "blue"), name = "") +
    scale_linetype_manual(values = c("True Coefficient" = "dashed", "True 95% CI" = "dashed"), name = "") +
    ylim(y_min - y_margin, y_max + y_margin)
})


if(dist_name=="GA") {dist_name_out="Gamma"} else if (dist_name=="NB") {dist_name_out="Negative Binomial"} else if (dist_name=="LO") {dist_name_out="Binomial"} else {dist_name_out="Unknown"}
# Display all three plots together
p=grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3
  , top = paste("Multivariate ",dist_name_out, " (T=5) Coefficient Estimates with 95% CI | rho=",rho," | skew=", skew_out, "| sims=",num_outer_sims,sep=""))

ggsave(p, filename = paste("Charts/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,".png",sep=""), width = 9, height = 4)


# Prepare data for ggplot
vs2_data <- data.frame(
  Score = c(vs2_all),
  Model = rep(rownames(sigmas), each = nrow(vs2_all)),
  Type = "Unweighted"
)

vs2_wt_data <- data.frame(
  Score = c(vs2_wt_all),
  Model = rep(rownames(sigmas), each = nrow(vs2_wt_all)),
  Type = "Weighted"
)

# Create ggplot boxplots with separate facets
p_vs2 <- ggplot(rbind(vs2_data, vs2_wt_data), aes(x = factor(Model, levels = rownames(sigmas)), y = Score, fill = Model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Type, scales = "free_y") +
  labs(
    title = paste("Variogram Scores across Models | rho=", rho, " | sims=", num_outer_sims, sep=""),
    x = "Model",
    y = "Variogram Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

print(p_vs2)

ggsave(p_vs2, filename = paste("Charts/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_VariogramScores.png",sep=""), width = 10, height = 6)



rownames(loglik_count_table_extended)=c("GLM","GEE","GAMLSS","LME4","GAMM","VineCopula")
colnames(loglik_count_table_extended)=c("LogLik","DF","AIC","BIC")
write.csv(loglik_count_table_extended, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_LogLik.csv",sep=""))
write.csv(coefficients_table_extended, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_Coef.csv",sep=""))
write.csv(ses_table_extended, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_SE.csv",sep=""))
write.csv(true_sim, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_TrueSim.csv",sep=""))
write.csv(timing_summary, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_Timing.csv",sep=""), row.names=FALSE)

write.csv(vs2_all, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_VariogramScores.csv",sep=""), row.names=FALSE)

print("Simulation complete.")
