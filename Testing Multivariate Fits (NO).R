library(VineCopula); source("common_functions.R")

simFit=function(n,rho=0,sims=100,data_in=NULL,mean_in=0,sigma_in=1,x1_in=0,x2_in=0,coef_in){
    coef_out=matrix(0,nrow=sims,ncol=3)
    dispersion_out=rep(0,sims)
    x1=x1_in
    x2=x2_in

    for(i in 1:sims) {
        simvine=simRVine(n=n,rho=rho,data_in=data_in)
        out=simvine[[1]]
        loglik_vine=simvine[[2]]
        out_adj=qnorm(out,mean=mean_in,sd=sigma_in)+matrix(rep(x1*coef_in[1],d),ncol=d)+matrix(rep(x2*coef_in[2],d),ncol=d)

        random_variable=c(out_adj[,1],out_adj[,2],out_adj[,3],out_adj[,4],out_adj[,5])
        patient=as.factor(rep(1:n,d))
        x1_long=rep(x1,d)
        x2_long=as.factor(rep(x2,d))
        data_long=data.frame(random_variable,patient,x1_long,x2_long)

        data_long=data_long[order(data_long$patient),]

        glm_fit=glm(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian())
        coef_out[i,]= glm_fit$coefficients
        dispersion_out[i]=summary(glm_fit)$dispersion
        #capture GLM outputs into a matrix
    }
    result=cbind(colMeans(coef_out),apply(coef_out,2,sd))
    colnames(result)=c("Mean","SD")
    rownames(result)=c("Intercept","X1","X2")
    return(list(result,loglik_vine,par=simvine[[3]],glm_dispersion=mean(dispersion_out)))
}

#Estimate true values
print("Estimating true parameter values...")
x1=runif(n); x2=as.numeric(runif(n)>0.5)
true_sim=simFit(n=1000,rho=rho,sims=true_sims,mean_in=sim_mean,sigma_in=sim_sigma,x1_in=x1,x2_in=x2,coef_in=c(1,1))[[1]]

sims=list()
coef_sum = ses_sum = coef_count = ses_count =loglik_count=loglik_sum= NULL

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

  print(paste("Simulation run:", j, "of", num_outer_sims))

  set.seed(100)
  library(glmtoolbox)
  library(gamlss)
  library(mgcv)
  library(lme4)
  out=pnorm(mvtnorm::rmvnorm(n=n,sigma=matrix(c(1,rho,rho^2,rho^3,rho^4,
                      rho,1,rho,rho^2,rho^3,
                      rho^2,rho,1,rho,rho^2,
                      rho^3,rho^2,rho,1,rho,
                      rho^4,rho^3,rho^2,rho,1),nrow=5)))
  x1=runif(n)
  x2=as.numeric(runif(n)>0.5)
  out_adj=qnorm(out,mean=1)+x1+x2

  random_variable=c(out_adj[,1],out_adj[,2],out_adj[,3],out_adj[,4],out_adj[,5])
  patient=as.factor(rep(1:n,d))
  x1_long=rep(x1,d)
  x2_long=as.factor(rep(x2,d))
  data_long=data.frame(random_variable,patient,x1_long,x2_long)

  data_long=data_long[order(data_long$patient),]

  print("Fitting models...")

  # GLM
  time_glm <- system.time(invisible(capture.output(model_glm<-           glm(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian()))))[3]
  model_times$GLM[j] <- time_glm
  
  # GEE
  time_gee <- system.time(invisible(capture.output(model_gee<-        glmgee(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian(), id=patient, corstr="Unstructured"))))[3]
  model_times$GEE[j] <- time_gee
  
  # RE-NOSIG
  time_renosig <- system.time(invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~1+x1_long+as.factor(x2_long)+random(as.factor(patient)), data=as.data.frame(out_adj), family=NO()))))[3]
  model_times$RENOSIG[j] <- time_renosig
  
  # LME4
  time_lme4 <- system.time(invisible(capture.output(model_lme4 <-       lmer(formula=random_variable~1+x1_long+as.factor(x2_long)+ (1|patient), data=as.data.frame(out_adj)))))[3]
  model_times$LME4[j] <- time_lme4
  
  # GAMM
  time_gamm <- system.time(invisible(capture.output(model_gamm <-        gamm(formula=random_variable~1+x1_long+as.factor(x2_long), random=list(patient=~1), data=as.data.frame(out_adj), family=gaussian))))[3]
  model_times$GAMM[j] <- time_gamm

  residuals_matrix <- cbind(
    model_glm$residuals[as.character(1:1000)],
    model_glm$residuals[as.character(1001:2000)],
    model_glm$residuals[as.character(2001:3000)],
    model_glm$residuals[as.character(3001:4000)],
    model_glm$residuals[as.character(4001:5000)]
  )

  coef_in=model_glm$coefficients
  time_vine <- system.time(simvinefit<-simFit(n=1000,rho=0,sims=100,data_in=pnorm(residuals_matrix),mean_in=coef_in[1],sigma_in=sigma(model_glm),x1_in=x1,x2_in=x2,coef_in=coef_in[2:3]))[3]
  model_times$VineCopula[j] <- time_vine
 
  results_table=list()

  print("Extracting results...")
  invisible(capture.output(results_table[[1]]<-summary(model_glm)$coeff[,1:2]))
  invisible(capture.output(results_table[[2]]<-summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2]))
  invisible(capture.output(results_table[[3]]<-cbind(summary(model_re_nosig)[1:3],summary(model_re_nosig)[5:7])))
  #esults_table[[4]]=cbind(summary(model_re_np)[1:4],summary(model_re_np)[8:11])
  invisible(capture.output(results_table[[5]]<-summary(model_lme4)$coefficients[,c(1,2)]))
  invisible(capture.output(results_table[[6]]<-cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix)))))

  rownames(results_table[[1]])=rownames(results_table[[2]])=rownames(results_table[[3]])=rownames(results_table[[5]])=rownames(results_table[[6]])=c("(Intercept)","x1_long","as.factor(x2_long)1")
  colnames(results_table[[1]])=colnames(results_table[[2]])=colnames(results_table[[3]])=colnames(results_table[[5]])=colnames(results_table[[6]])=c("Estimate","Std. Error")
  names(results_table)=c("GLM","GEE","GAMLSS","RE-NP","LME4","GAMM")

    ###Calculating effective degrees of freedom from Donohue
    X<-getME(model_lme4,name="X")[,1:2]
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
        , logLik(model_gamm$lme)
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

       conv_check=c(
        model_glm$converged,
        model_gee$converged,
        model_re_nosig$converged,
        if (!is.null(model_lme4)) !any( grepl("failed to converge", model_lme4@optinfo$conv$lme4$messages) ) else NA,
        if (!is.null(model_gamm)) !any( grepl("converge", warnings(model_gamm))) else NA,
        NA_real_        

      )

      # Store convergence rates
      conv_rates$GLM[j] <- as.numeric(model_glm$converged)
      conv_rates$GEE[j] <- as.numeric(model_gee$converged)
      conv_rates$RENOSIG[j] <- as.numeric(model_re_nosig$converged)
      conv_rates$LME4[j] <- as.numeric(!any( grepl("failed to converge", model_lme4@optinfo$conv$lme4$messages) ))
      conv_rates$GAMM[j] <- as.numeric(!any( grepl("converge", warnings(model_gamm))))



  sims[[j]]=coefficients_table

  rownames(coefficients_table)=rownames(ses_table)=rownames(loglik_table)=rownames(sigmas)=names(correlations)=c("GLM","GEE","GAMLSS","LME4","GAMM","VineCopula")
  colnames(loglik_table)=c("LogLik","DF","AIC","BIC")

  ########## VARIOGRAM SCORES
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
    
      
      #Calculate variogram scores
      #vs2[model]    = vs_sample(y=y,dat=sim_model_out[[model]],p=2)
      #vs2_wt[model] = vs_sample(y=y,dat=sim_model_out[[model]],p=2,w_vs=w_vs)


      #########################################

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
conv_check_final <- conv_check_sum / num_outer_sims #Correct convergence rate for each model

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
    mean(model_times$GLM),
    mean(model_times$GEE),
    mean(model_times$RENOSIG),
    mean(model_times$LME4),
    mean(model_times$GAMM),
    mean(model_times$VineCopula)
  ),
  Convergence_Rate = c(
    mean(conv_rates$GLM),
    mean(conv_rates$GEE),
    mean(conv_rates$RENOSIG),
    mean(conv_rates$LME4),
    mean(conv_rates$GAMM),
    NA_real_
  )
)


######PLOTTING

library(ggplot2)
library(gridExtra)

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

skew_out=0; dist_name_out="Normal"
# Display all three plots together
p=grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3
  , top = paste("Multivariate Normal (T=5) Coefficient Estimates with 95% CI | rho=",rho," | sims=",num_outer_sims,sep=""))

ggsave(p, filename = paste("Charts/Multivariate_Normal_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,".png",sep=""), width = 9, height = 4)

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
write.csv(loglik_count_table_extended, file = paste("Charts/ChartData/Multivariate_Normal_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_LogLik.csv",sep=""))
write.csv(coefficients_table_extended, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_Coef.csv",sep=""))
write.csv(ses_table_extended, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_SE.csv",sep=""))
write.csv(true_sim, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_TrueSim.csv",sep=""))
write.csv(timing_summary, file = paste("Charts/ChartData/Multivariate_Normal_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_Timing.csv",sep=""), row.names=FALSE)
write.csv(vs2_all, file = paste("Charts/ChartData/Multivariate_",dist_name_out,"_T5_rho",rho,"_sims",num_outer_sims,"_skew",skew_out,"_VariogramScores.csv",sep=""), row.names=FALSE)

print("Simulation complete.")
