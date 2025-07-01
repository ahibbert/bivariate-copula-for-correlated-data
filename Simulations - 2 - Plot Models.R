###########CHOOSE DISTRIBUTION AND PARAMETERS######################################

require(latex2exp)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
source("link_functions.R")
plotVersusTrue <- function (limits,inputs,true,x_input,xlab,ylab,scaled=FALSE,type="ALL",plotTrue=TRUE,remove_outliers=FALSE) {
  
  #inputs=t1intercepts; true=mu1;x_input=tau; xlab=NA; ylab=NA; scaled=FALSE; limits=NA
  
  if (scaled==TRUE) {
    inputs = inputs / true - 1
    true = rep(0,nrow(inputs))
    inputs[!is.finite(inputs)]=0
  }
  
  inputs=as.data.frame(inputs)
  smooth_fit=inputs*0
  
  for (col_name in colnames(inputs)) {
    y=inputs[,col_name]
    weights=rep(1,length(y))
    if((col_name=="summary_re_nosig" | col_name=="summary_re_np") & remove_outliers==TRUE) {
      weights=as.numeric(y<mean(y,trim=0.2)+abs(mean(y,trim=0.2)*10) & y>mean(y,trim=0.2)-abs(mean(y,trim=0.2)*10))
    }
    
    if(all(is.na(y))) {smooth_fit[,col_name]=y} else {
      lo <- loess(y~x_input,weights=weights)  
      smooth_fit[,col_name]=predict(lo)
    }
  }
  
  if(plotTrue==TRUE) {
    true_smooth=predict(loess(true~x_input))  
  }
  
  plot<-ggplot() + labs(x = xlab, y=ylab) +
    {if(!(is.na(limits[1])||is.na(limits[2]))){xlim(limits[1],limits[2])}} +
    {if(!(is.na(limits[3])||is.na(limits[4]))){ylim(limits[3],limits[4])}} +
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_glm'], color="GLM"),linetype = 3) + 
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_gee'], color="GEE"),linetype = 3,position=position_jitter()) + #
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_lme4'], color="LME4"),linetype = 5) + 
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_gamm'], color="GAMM"),linetype = 5) +
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_re_nosig'], color="GAMLSS (4)"),linetype = 5) +
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_re_np'], color="GAMLSS NP (5)"),linetype = 5) +
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_cop'], color="GJRM (C)"),linetype = 4) +
    geom_line(data=inputs, aes(x=x_input, y=smooth_fit[,'summary_cop_n'], color="GJRM (N)"),linetype = 4) +
    {if(plotTrue==TRUE){geom_line(data=inputs, aes(x=x_input, y=true_smooth, color="True"))}} +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMM","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 8, name = "Dark2"),"#000000"))
  return(plot)
}

#bt_mode=FALSE;files_in=c("Data/results_combined_B1_B2_NO_1000_2024-11-27.RData","Data/results_combined_B1_B2_PO_1000_2024-12-04.RData","Data/results_combined_B1_B2_GA_1000_2024-11-28.RData","Data/results_combined_B1_B2_LO_1000_2024-11-27.RData")
#bt_mode=TRUE; files_in=c("Data/results_combined_B1_Bt_NO_1000_2024-12-05.RData","Data/results_combined_B1_Bt_PO_1000_2024-12-05.RData","Data/results_combined_B1_Bt_GA_1000_2024-11-28.RData","Data/results_combined_B1_Bt_LO_1000_2024-11-26.RData")

files=files_in
multiplot=TRUE; plotcount=plot_count_lik=0
for (filename in files)
{
  #filename=files[2]
  load(filename); print(filename)
  dist=substr(filename,29,30); print(dist)
  ###################### 1. DATA SETUP##################################################
  print("DATA SETUP STAGE")
  
  if ( (multiplot==FALSE) || (multiplot==TRUE && (!exists("first_run"))) ) {
    plot_count_lik<-plotcount<-0; bias_plots <- lik_plots_skew <- skew_error_plots <- skew_bias_plots <- error_plots <- lik_plots<-list(); first_run=FALSE;bias_all<-list();
  }
  
  check_length=c()
  for (i in 1:length(results_combined)) {
    check_length[i]<-nrow(results_combined[[i]])  
  }
  
  if(!(min(check_length)==max(check_length))) {print("WARNING: Some simulations are incomplete")}
  
  results <- results_combined[check_length==max(check_length)]
  
  #Take out parameters
  parameters=matrix(0,nrow=length(results),ncol=6)
  for (z in 1:length(results)) {
    parameters[z,]=results[[z]][nrow(results[[z]]),1:6]
  }
  colnames(parameters)<-c("n","a","b","c","mu1","mu2")
  
  ###Wrapper for easy plotting
  
  set.seed(1000)
  options(scipen=999)
  
  tau<-emp_cor<-mu1<-mu2<-skew<-vector(length=length(results))
  t1intercepts<-t1error<-t2intercepts<-t2error<-aic<-bic<-adj_bic<-aic_4<-loglik<-matrix(ncol=(nrow(results[[1]])-2),nrow=length(results))
  
  for (i in 1:length(results)) {
    t1intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_1"])
    t2intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_2"])
    t1error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b1"])
    t2error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b2"])
    
    loglik[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"LogLik"])
    
    aic[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"AIC"])
    aic_4[i,]=-2*t(results[[i]][1:(nrow(results[[i]])-2),"LogLik"])+4*t(results[[i]][1:(nrow(results[[i]])-2),"EDF"])
    
    adj_bic[i,]=-2*t(results[[i]][1:(nrow(results[[i]])-2),"LogLik"])+log(1000)*t(results[[i]][1:(nrow(results[[i]])-2),"EDF"])
    
    bic[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"BIC"])
    
    tau[i]=results[[i]][(nrow(results[[i]])-1),6]/100
    emp_cor[i]=results[[i]][(nrow(results[[i]])-1),7]/100
    
    colnames(t1intercepts)<-colnames(t2intercepts)<-colnames(aic_4)<-colnames(t1error)<-colnames(t2error)<-colnames(loglik)<-colnames(aic)<-colnames(adj_bic)<-colnames(bic)<-colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))
  }
  
  loglik[,"summary_gamm"]=aic[,"summary_gamm"]=bic[,"summary_gamm"]=aic_4[,"summary_gamm"]=adj_bic[,"summary_gamm"]=NA #Temporary until we fix this in calculations
  
  
  #Theoretical errors
  if(dist=="GA") {
    #Parameters
    mu1=parameters[,"a"]*parameters[,"mu1"]
    mu2=parameters[,"a"]*parameters[,"mu2"]
    if(bt_mode==TRUE) {
      mu2=log(mu2/mu1)
    }
    #Errors
    load(file="Data/numDerivResults_20231127.rds")
    load("Data/se_mles_20231127_n100sims20_ALL.rds") #se_mles
    
    se_sim=se_mles[,c(1,2,3)]/sqrt(10)
    se_nd=numDerivResults[,c(1,2,5)]
    
    se_final=se_nd
    for (i in 1:3) {
      is.nan(se_final[,i])
      se_final[is.nan(se_final[,i]),i]=se_sim[is.nan(se_final[,i]),i]
    }
    
    trueSE<-se_final
    
    skew=2/sqrt(parameters[,"a"])
  }
  if(dist=="LO") {
    #Parameters
    mu1=parameters[,"mu1"]
    mu2=parameters[,"mu2"]
    
    
    p=mu1;q=1-mu1
    skew1=(q-p)/sqrt(p*q)
    p=mu2;q=1-mu2
    skew2=(q-p)/sqrt(p*q)
    skew=(skew2+skew1)/2
    
    if(bt_mode==TRUE) {
      mu2=logit(mu2)-logit(mu1)
    }
    
    load(file="Data/lo_mle")
    trueSE=sqrt(merge(parameters,lo_mle,by=c("a","b","c","mu1","mu2"))[,c("var_B1","var_B2","var_Bt")])
    colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")
    
  }
  if(dist=="NO") {
    #Parameters
    mu1=parameters[,"mu1"]
    mu2=parameters[,"mu2"]
    
    if(bt_mode==TRUE) {
      mu2=mu2-mu1
    }
    skew=mu1*0
    #Errors
    trueSE<-t(rbind((parameters[,"a"]/sqrt(parameters[,"n"]))
                    ,(parameters[,"b"]/sqrt(parameters[,"n"]))
                    ,sqrt((parameters[,"a"]^2)+(parameters[,"b"]^2)-2*parameters[,"a"]*parameters[,"b"]*parameters[,"c"])/sqrt(parameters[,"n"])))
    colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")
  }
  if(dist=="PO") {
    load("Data/nb_mle")
    trueSE<-matrix(ncol=3,nrow=length(results))
    
    for (i in 1:length(results)) {
      trueSE[i,]<-c(results[[i]]["actuals",c("se_b1","se_b2","LogLik")]) 
      mu1[i]=c(results[[i]]["actuals",c("b_1")]) 
      mu2[i]=c(results[[i]]["actuals",c("b_2")]) 
      skew[i]=results[[i]][(nrow(results[[i]])-1),8]/10000
    }
    colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")
    
    trueSE[,3]=sqrt(merge(parameters,nb_mle,by=c("a","b","c","mu1","mu2"))[,"var_Bt"])
    
    if(bt_mode==TRUE) {mu2=log(mu2/mu1)}
      
  }
  
  if(bt_mode==TRUE) {
    t2intercepts[,grepl("cop",colnames(t1intercepts))]= t2intercepts[,grepl("cop",colnames(t1intercepts))] - t1intercepts[,grepl("cop",colnames(t1intercepts))]   
  }
  
  
  ###################### 2.PLOT SETUP######################
  
  print("--BIAS AND ERROR PLOTS STAGE")
  
  library(latex2exp)
  limits_bias = c(.1,0.7,-1,if(dist=="LO"&bt_mode==TRUE) {2} else {1}); xlab=TeX("Kendall's \\tau")
  limits_error <- c(limits_bias[1:2],0,if(dist=="PO"){.25} else if (dist=="LO") {.4} else{.125})
  limits_bias_skew = c(min(skew),max(skew),-1,1); xlabskew=TeX("Skewness")
  limits_error_skew <- c(limits_bias_skew[1:2],0,if(dist=="PO"){.25}else{.125})
  
  limits_bias=limits_error=c(.2,.8,NA,NA)
  #if(dist=="NO") {tau=parameters[,"c"]; xlab="Pearson Correlation"}
  
  plotcount=plotcount+1
  bias_plots[[plotcount]]<- plotVersusTrue(limits_bias
                 ,if(dist=="NO"){t1intercepts}else if(dist=="LO") {logit_inv(t1intercepts)} else{exp(t1intercepts)}
                 ,mu1
                 ,tau
                 ,xlab
                 ,ylab=TeX("$(\\hat{\\mu_1}/\\mu_1)-1$")
                 , scaled=TRUE)
  error_plots[[plotcount]]<- plotVersusTrue(limits_error
                                            ,t1error
                                            ,trueSE[,"mu1_se"]
                                            ,tau
                                            ,xlab
                                            ,ylab=TeX("$SE(\\hat{\\beta_{1}})$")
                                            ,scaled=FALSE
                                            ,remove_outliers=TRUE)
  
  skew_bias_plots[[plotcount]]<- plotVersusTrue(limits_bias_skew
                                           ,if(dist=="NO"){t1intercepts}else if(dist=="LO") {logit_inv(t1intercepts)} else{exp(t1intercepts)}
                                           ,mu1
                                           ,skew
                                           ,xlabskew
                                           ,ylab=TeX("$(\\hat{\\mu_1}/\\mu_1)-1$")
                                           , scaled=TRUE)
  skew_error_plots[[plotcount]]<- plotVersusTrue(limits_error_skew
                                            ,t1error
                                            ,trueSE[,"mu1_se"]
                                            ,skew
                                            ,xlabskew
                                            ,ylab=TeX("$SE(\\hat{\\beta_{1}})$")
                                            ,scaled=FALSE
                                            ,remove_outliers=TRUE)
  
  plotcount=plotcount+1
  bias_plots[[plotcount]]<- plotVersusTrue(limits_bias
                 ,if(dist=="NO"){cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t2intercepts)])}
                  else{cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t1intercepts)])}
                 ,if(dist=="LO"&bt_mode==FALSE) {logit(mu2)} else {mu2}
                 ,tau
                 ,xlab
                 ,ylab=TeX("$(\\hat{\\beta_t}/\\beta_t)-1$")
                 ,scaled=TRUE)
  error_plots[[plotcount]]<- plotVersusTrue(limits_error
                                            ,t2error
                                            ,if(bt_mode==TRUE){trueSE[,"mu2_se_Bt"]}else{trueSE[,"mu2_se_B2"]}
                                            ,tau
                                            ,xlab
                                            ,ylab=TeX("$SE(\\hat{\\beta_{t}})$")
                                            ,scaled=FALSE
                                            ,remove_outliers=TRUE)
  skew_error_plots[[plotcount]]<- plotVersusTrue(limits_error_skew
                                            ,t2error
                                            ,if(bt_mode==TRUE){trueSE[,"mu2_se_Bt"]}else{trueSE[,"mu2_se_B2"]}
                                            ,skew
                                            ,xlabskew
                                            ,ylab=TeX("$SE(\\hat{\\beta_{t}})$")
                                            ,scaled=FALSE
                                            ,remove_outliers=TRUE)
  
  skew_bias_plots[[plotcount]]<- plotVersusTrue(limits_bias_skew
                                           ,if(dist=="NO"){cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t2intercepts)])}
                                           else{cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t1intercepts)])}
                                           ,if(dist=="LO"&bt_mode==FALSE) {logit(mu2)} else {mu2}
                                           ,mu2
                                           ,skew
                                           ,xlabskew
                                           ,ylab=TeX("$(\\hat{\\mu_2}/\\mu_2)-1$")
                                           ,scaled=TRUE)
  
  print("LIKELIHOOD PLOTS")
  ########Likelihoods
  limits_lik <- c(.2,.8,NA,NA)
  limits_lik_skew <- c(.2,.8,NA,NA)

  #if (dist=="PO") {limits_lik <- c(limits_bias[1:2],-5000,1000);limits_lik_skew<-c(limits_bias_skew[1:2],-5000,1000)}
  
  plot_count_lik=plot_count_lik+1
  lik_plots[[plot_count_lik]]<- plotVersusTrue(limits_lik
                                               ,loglik
                                               ,NA
                                               ,tau
                                               ,xlab
                                               ,ylab="LogLik"
                                               ,scaled=FALSE
                                               ,plotTrue = FALSE
                                               ,remove_outliers = TRUE)
  lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(limits_lik_skew
                                               ,loglik
                                               ,NA
                                               ,skew
                                               ,xlabskew
                                               ,ylab="LogLik"
                                               ,scaled=FALSE
                                               ,plotTrue = FALSE
                                               ,remove_outliers = TRUE)
  plot_count_lik=plot_count_lik+1
  lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                               ,aic
                                               ,NA
                                               ,tau
                                               ,xlab
                                               ,ylab="AIC"
                                               ,scaled=FALSE
                                               ,plotTrue = FALSE
                                               ,remove_outliers = TRUE)
  lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                    ,aic
                                                    ,NA
                                                    ,skew
                                                    ,xlabskew
                                                    ,ylab="AIC"
                                                    ,scaled=FALSE
                                                    ,plotTrue = FALSE
                                                    ,remove_outliers = TRUE)
  
  plot_count_lik=plot_count_lik+1
  lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                               ,aic_4#bic
                                               ,NA
                                               ,tau
                                               ,xlab
                                               ,ylab="GAIC (4)"
                                               ,scaled=FALSE
                                               ,plotTrue = FALSE
                                               ,remove_outliers = TRUE)
  
  lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                    ,aic_4#bic
                                                    ,NA
                                                    ,skew
                                                    ,xlabskew
                                                    ,ylab="GAIC (4)"
                                                    ,scaled=FALSE
                                                    ,plotTrue = FALSE
                                                    ,remove_outliers = TRUE)
  plot_count_lik=plot_count_lik+1
  lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                               ,adj_bic#bic
                                               ,NA
                                               ,tau
                                               ,xlab
                                               ,ylab="BIC"
                                               ,scaled=FALSE
                                               ,plotTrue = FALSE
                                               ,remove_outliers = TRUE)
  
  lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                    ,adj_bic#bic
                                                    ,NA
                                                    ,skew
                                                    ,xlabskew
                                                    ,ylab="BIC"
                                                    ,scaled=FALSE
                                                    ,plotTrue = FALSE
                                                    ,remove_outliers = TRUE)
  
}


###Bias / Error
if(bt_mode==TRUE) {
  ggarrange(bias_plots[[2]],error_plots[[2]]
            ,bias_plots[[4]],error_plots[[4]]
            ,bias_plots[[6]],error_plots[[6]]
            ,bias_plots[[8]],error_plots[[8]]
            ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("NO","NO","NB","NB","GA","GA","LO","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Tau
  ggsave(file=paste("Charts/simulation_bias_plus_error_Bt",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=11,dpi=900)
  
} else {
  ggarrange(bias_plots[[1]],error_plots[[1]]
            ,bias_plots[[3]],error_plots[[3]]
            ,bias_plots[[5]],error_plots[[5]]
            ,bias_plots[[7]],error_plots[[7]]
            ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("NO","NO","NB","NB","GA","GA","LO","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Tau
  ggsave(file=paste("Charts/simulation_bias_plus_error_B2",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=11,dpi=900)
  
  ggarrange(plotlist=skew_bias_plots[c(3,5,7)],ncol=3,nrow=1,labels=c("NB","GA","LO"),common.legend = TRUE,hjust=0.1,font.label = list(size = 12)) + bgcolor("white") + border(color = "white") # Likelihoods x Tau
  ggsave(file=paste("Charts/simulation_bias_skew",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=4,dpi=900)
 
  ggarrange(plotlist=lik_plots[c(1,5,9,13,4,8,12,16)]
            ,common.legend=TRUE, ncol=4, nrow=2, labels=c("NO","NB","GA","LO","NO","NB","GA","LO"),hjust=0.1,font.label = list(size = 12)) + bgcolor("white") + border(color = "white") # Likelihoods x Tau
  ggsave(file=paste("Charts/simulation_loglik_BIC_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=6,dpi=900)
  
  ggarrange(plotlist=lik_plots[c(2,6,10,14,3,7,11,15)]
            ,common.legend=TRUE, ncol=4,nrow=2, labels=c("NO","NB","GA","LO","NO","NB","GA","LO"),hjust=0.1,font.label = list(size = 12)) + bgcolor("white") + border(color = "white") # Likelihoods x Tau
  ggsave(file=paste("Charts/simulation_AIC_GAIC_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=6,dpi=900)
}

#############Bias v skew table

#bias_all[[plotcount/2]]<-cbind(trunc(tau*10,1)*10,trunc(skew),(if(dist=="NO"){t1intercepts}else{exp(t1intercepts)}/mu1)-1)

#bias_all_combined<-rbind(bias_all[[1]],bias_all[[2]])
#bias_all<-bias_all_combined

#bias_table_list<-list()
#for (i in 1:3) {
#  dataset<- bias_all[[plotcount/2]][,c(1,2,4+i)]
#  colnames(dataset)<-c("tau","skew","bias")
#  summary<-aggregate(dataset[,"bias"] ~ dataset[,"tau"] + dataset[,"skew"], data = dataset, mean, na.rm = TRUE)
#  colnames(summary)<-c("tau","skew","bias")
#  bias_table_list[[i]]<-round(xtabs(summary[,"bias"] ~ summary[,"tau"] + summary[,"skew"]),2)
#}

