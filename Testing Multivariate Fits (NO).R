######################################### NORMAL

library(VineCopula)
n=1000;d=5

simRVine=function(n,rho=0,data_in=NULL){

    # define 5-dimensional R-vine tree structure matrix
    Matrix <- c(5, 2, 3, 1, 4,
                0, 2, 3, 4, 1,
                0, 0, 3, 4, 1,
                0, 0, 0, 4, 1,
                0, 0, 0, 0, 1)
    Matrix <- matrix(Matrix, 5, 5)
    # define R-vine pair-copula family matrix
    family <- c(0, 1, 1, 1, 1,
                0, 0, 1, 1, 1,
                0, 0, 0, 1, 1,
                0, 0, 0, 0, 1,
                0, 0, 0, 0, 0)
    family <- matrix(family, 5, 5)
    # define R-vine pair-copula parameter matrix
    
    if(rho==0){
        #automatically select parameters for normal copula
         par=RVineCopSelect(data_in, familyset=1,Matrix)$par
    }   else {
        par <- c(0, rho^4, rho^3, rho^2, rho,
            0, 0, rho^3, rho^2, rho,
            0, 0, 0, rho^2, rho,
            0, 0, 0, 0, rho,
            0, 0, 0, 0, 0)
    }

    par <- matrix(par, 5, 5)
    # define second R-vine pair-copula parameter matrix
    par2 <- matrix(0, 5, 5)

    ## define RVineMatrix object
    RVM <- RVineMatrix(Matrix = Matrix, family = family,
                    par = par, par2 = par2,
                    names = c("V1", "V2", "V3", "V4", "V5"))

    ## see the object's content or a summary
    #str(RVM)
    #summary(RVM)

    ### inspect the model using plots
    #if (FALSE) plot(RVM)  # tree structure
    #contour(RVM)  # contour plots of all pair-copulas

    ## simulate from the vine copula model
    out=RVineSim(n, RVM)
    return(out)

}

simFit=function(n,rho=0,sims=100,data_in=NULL,mean_in=0,sigma_in=1,x1_in=0,x2_in=0,coef_in){
    coef_out=matrix(0,nrow=sims,ncol=3)
    for(i in 1:sims) {
        out=simRVine(n=n,rho=rho,data_in=data_in)
        out_adj=qnorm(out,mean=mean_in,sd=sigma_in)+matrix(rep(x1*coef_in[1],d),ncol=d)+matrix(rep(x2*coef_in[2],d),ncol=d)

        random_variable=c(out_adj[,1],out_adj[,2],out_adj[,3],out_adj[,4],out_adj[,5])
        patient=as.factor(rep(1:n,d))
        x1_long=rep(x1,d)
        x2_long=as.factor(rep(x2,d))
        data_long=data.frame(random_variable,patient,x1_long,x2_long)

        data_long=data_long[order(data_long$patient),]

        coef_out[i,]= glm(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian())$coefficients
        #capture GLM outputs into a matrix
    }
    result=cbind(colMeans(coef_out),apply(coef_out,2,sd))
    colnames(result)=c("Mean","SD")
    rownames(result)=c("Intercept","X1","X2")
    return(result)
}

#final_high_cor=simFit(n=1000,rho=0.75,sims=100)
#final_low_cor=simFit(n=1000,rho=0.25,sims=100)

rho=0.5; num_outer_sims=5
sims=list()
  coef_sum = ses_sum = coef_count = ses_count = NULL
for (j in 1:num_outer_sims) {

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

  model_glm=           glm(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian())
  model_gee=        glmgee(        random_variable~1+x1_long+as.factor(x2_long), data=data_long, family=gaussian(), id=patient, corstr="Unstructured")
  model_re_nosig <- gamlss(formula=random_variable~1+x1_long+as.factor(x2_long)+random(as.factor(patient)), data=as.data.frame(out_adj), family=NO())
  model_lme4 <-       lmer(formula=random_variable~1+x1_long+as.factor(x2_long)+ (1|patient), data=as.data.frame(out_adj))
  model_gamm =        gamm(formula=random_variable~1+x1_long+as.factor(x2_long), random=list(patient=~1), data=as.data.frame(out_adj), family=gaussian)

  residuals_matrix <- cbind(
    model_glm$residuals[as.character(1:1000)],
    model_glm$residuals[as.character(1001:2000)],
    model_glm$residuals[as.character(2001:3000)],
    model_glm$residuals[as.character(3001:4000)],
    model_glm$residuals[as.character(4001:5000)]
  )

  coef_in=summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2]

  final_vinecop=simFit(n=1000,rho=0,sims=100,data_in=pnorm(residuals_matrix),mean_in=coef_in[1,1],sigma_in=sigma(model_glm),x1_in=x1,x2_in=x2,coef_in=c(coef_in[2,1],coef_in[3,1]))

  results_table=list()

  results_table[[1]]=summary(model_glm)$coeff[,1:2]
  results_table[[2]]=summary(model_gee)$coefficients[1:(nrow(summary(model_gee)$coefficients)-2),1:2]
  results_table[[3]]=cbind(summary(model_re_nosig)[1:3],summary(model_re_nosig)[5:7])
  #esults_table[[4]]=cbind(summary(model_re_np)[1:4],summary(model_re_np)[8:11])
  results_table[[5]]=summary(model_lme4)$coefficients[,c(1,2)]
  results_table[[6]]=cbind(summary(model_gamm$lme)$coefficients[[1]],sqrt(diag(model_gamm$lme$varFix)))

  rownames(results_table[[1]])=rownames(results_table[[2]])=rownames(results_table[[3]])=rownames(results_table[[5]])=rownames(results_table[[6]])=c("(Intercept)","x1_long","as.factor(x2_long)1")
  colnames(results_table[[1]])=colnames(results_table[[2]])=colnames(results_table[[3]])=colnames(results_table[[5]])=colnames(results_table[[6]])=c("Estimate","Std. Error")
  names(results_table)=c("GLM","GEE","RE-NOSIG","RE-NP","LME4","GAMM")

  #    dfs=c( ####Come back to DF
  #      n*2-df.residual(model_glm)
  #      , n*2-model_gee$df.residual
  #      , model_re_nosig$df.fit
  #      , lme_EDF
  #      , lme_EDF #GAMM doesn't provide dfs so using GEE as proxy
  #    )
  #    logLiks=c(
  #      logLik(model_glm)
  #      , model_gee$logLik
  #      , logLik(model_re_nosig)
  #      , logLik(model_lme4)
  #      , logLik(model_gamm$lme)
  #    )

      #Creating summary tables
      coefficients_table= rbind(
        results_table[[1]][,1]
        , results_table[[2]][,1]
        , results_table[[3]][,1]
        , results_table[[5]][,1]
        , results_table[[6]][,1]
        , final_vinecop[,1]
      )
      ses_table= rbind(
        results_table[[1]][,2]
        , results_table[[2]][,2]
        , results_table[[3]][,2]
        , results_table[[5]][,2]
        , results_table[[6]][,2]
        , final_vinecop[,2]
      )
  #    loglik_table= rbind(
  #      c(logLiks[1],dfs[1])
  #      , c(logLiks[2],dfs[2])
  #      , c(logLiks[3],dfs[3])
  #      , c(logLiks[5],dfs[5])
  #      , c(logLiks[6],dfs[6])
  #    )

  sims[[j]]=coefficients_table

  rownames(coefficients_table)=rownames(ses_table)=c("GLM","GEE","RE-NOSIG","LME4","GAMM","VineCopula")

    if (is.null(coef_sum)) {
    coef_sum <- ifelse(is.na(coefficients_table), 0, coefficients_table)
    ses_sum <- ifelse(is.na(ses_table), 0, ses_table)
    coef_count <- !is.na(coefficients_table)
    ses_count <- !is.na(ses_table)
  } else {
    coef_sum <- coef_sum + ifelse(is.na(coefficients_table), 0, coefficients_table)
    ses_sum <- ses_sum + ifelse(is.na(ses_table), 0, ses_table)
    coef_count <- coef_count + !is.na(coefficients_table)
    ses_count <- ses_count + !is.na(ses_table)
  }

}

coefficients_table_extended <- coef_sum / coef_count
ses_table_extended <- ses_sum / ses_count
coefficients_table_extended[coef_count == 0] <- NA_real_
ses_table_extended[ses_count == 0] <- NA_real_


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
  y_min <- min(plot_data_list[[i]]$Lower)
  y_max <- max(plot_data_list[[i]]$Upper)
  y_range <- y_max - y_min
  y_margin <- y_range * 0.05  # 5% margin
  
  ggplot(plot_data_list[[i]], aes(x = Model, y = Coefficient)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
    labs(
      title = colnames(coefficients_table_extended)[i],
      x = "Model",
      y = "Coefficient Estimate"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5) +
    ylim(y_min - y_margin, y_max + y_margin)
})

# Display all three plots together
grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3
  , top = paste("Multivariate Normal (T=5) Coefficient Estimates with 95% CI | rho=",rho,sep=""))

# Alternatively, display them individually:
# plots[[1]]  # (Intercept)
# plots[[2]]  # x1_long
# plots[[3]]  # as.factor(x2_long)1
