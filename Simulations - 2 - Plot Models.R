###########CHOOSE DISTRIBUTION AND PARAMETERS######################################

require(latex2exp)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
source("link_functions.R")
plotVersusTrue <- function (limits,inputs,true,tau,xlab,ylab,scaled=FALSE,type="ALL",plotTrue=TRUE) {
  
  if (scaled==TRUE) {
    inputs = inputs / true-1
    true = true/true-1
  }
  
  inputs=as.data.frame(inputs)
  
  plot<-ggplot() + labs(x = xlab, y=ylab) +
    {if(!(is.na(limits[1])||is.na(limits[2]))){xlim(limits[1],limits[2])}} +
    {if(!(is.na(limits[3])||is.na(limits[4]))){ylim(limits[3],limits[4])}} +
    geom_smooth(data=inputs, aes(x=tau, y=summary_glm, color="GLM"),linetype = "dashed",se=FALSE) + 
    geom_smooth(data=inputs, aes(x=tau, y=summary_gee, color="GEE"),linetype = "dashed",se=FALSE,position=position_jitter()) + #
    geom_smooth(data=inputs, aes(x=tau, y=summary_lme4, color="LME4"),linetype = "dashed",se=FALSE) + 
    geom_smooth(data=inputs, aes(x=tau, y=summary_re_nosig, color="GAMLSS (4)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=inputs, aes(x=tau, y=summary_re_np, color="GAMLSS NP (5)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=inputs, aes(x=tau, y=summary_cop, color="GJRM (C)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=inputs, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),linetype = "dashed",se=FALSE) +
    {if(plotTrue==TRUE){geom_smooth(data=inputs, aes(x=tau, y=true, color="True"), span=1,se=FALSE)}} +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"#000000"))
  return(plot)
}

#B1/B2 results
#bt_mode=FALSE
#load("Data/results_combined_nointNO_1000_2024-03-15.RData"); dist="NO" ####NO
#load("Data/results_combined_B1_B2_PO_1000_2024-08-27.RData");dist="PO" #load("results_combined_noint_rangeupPO_1000_2024-03-19.RData"); dist="PO"
#load("Data/results_combined_nointGA_1000_2024-03-15.RData"); dist="GA" ####GA
#load("Data/results_combined_B1_B2_LO_1000_2024-10-30.RData");dist="LO"

#B0/Bt results
bt_mode=TRUE
#load("Data/results_combined_B0_BtNO_1000_2024-07-23.RData");dist="NO"
#load("Data/results_combined_B0_Bt_PO_1000_2024-08-28.RData");dist="PO" #load("results_combined_B0_BtPO_1000_2024-07-24.RData");dist="PO"
#load("Data/results_combined_B0_BtGA_1000_2024-07-23.RData");dist="GA"
load("Data/results_combined_B1_Bt_LO_1000_2024-10-30.RData");dist="LO"

multiplot=TRUE

#################################1. DATA SETUP##################################################

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
  trueSE<-numDerivResults[,c(1,2,5)]
  skew=2/sqrt(parameters[,"a"])
}

if(dist=="LO") {
  #Parameters
  mu1=parameters[,"mu1"]
  mu2=parameters[,"mu2"]
  
  if(bt_mode==TRUE) {
    mu2=logit(mu2)-logit(mu1)
  }
  
  p=mu1;q=1-mu1
  skew1=(q-p)/sqrt(p*q)
  p=mu2;q=1-mu2
  skew2=(q-p)/sqrt(p*q)
  skew=(skew2+skew1)/2
  
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
  #Errors
  trueSE<-t(rbind((parameters[,"a"]*sqrt(1-(parameters[,"c"]^2)))/sqrt(parameters[,"n"])
                  ,(parameters[,"b"]*sqrt(1-(parameters[,"c"]^2)))/sqrt(parameters[,"n"])
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


###################### PLOT SETUP######################

library(latex2exp)
limits_bias = c(.1,0.7,-1,if(dist=="LO"&bt_mode==TRUE) {2} else {1}); xlab=TeX("Kendall's \\tau")
limits_error <- c(limits_bias[1:2],0,if(dist=="PO"){.25} else if (dist=="LO") {.4} else{.125})
limits_bias_skew = c(min(skew),max(skew),-1,1); xlabskew=TeX("Skewness")
limits_error_skew <- c(limits_bias_skew[1:2],0,if(dist=="PO"){.25}else{.125})
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
                                          ,scaled=FALSE)

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
                                          ,scaled=FALSE)

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
                                          ,scaled=FALSE)
skew_error_plots[[plotcount]]<- plotVersusTrue(limits_error_skew
                                          ,t2error
                                          ,if(bt_mode==TRUE){trueSE[,"mu2_se_Bt"]}else{trueSE[,"mu2_se_B2"]}
                                          ,skew
                                          ,xlabskew
                                          ,ylab=TeX("$SE(\\hat{\\beta_{t}})$")
                                          ,scaled=FALSE)

skew_bias_plots[[plotcount]]<- plotVersusTrue(limits_bias_skew
                                         ,if(dist=="NO"){cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t2intercepts)])}
                                         else{cbind((t2intercepts)[,1:5],(t2intercepts)[,6:ncol(t1intercepts)])}
                                         ,if(dist=="LO"&bt_mode==FALSE) {logit(mu2)} else {mu2}
                                         ,mu2
                                         ,skew
                                         ,xlabskew
                                         ,ylab=TeX("$(\\hat{\\mu_2}/\\mu_2)-1$")
                                         ,scaled=TRUE)

########Likelihoods
limits_lik <- c(limits_bias[1:2],NA,NA)
limits_lik_skew <- c(limits_bias_skew[1:2],NA,NA)
if (dist=="PO") {limits_lik <- c(limits_bias[1:2],-5000,1000);limits_lik_skew<-c(limits_bias_skew[1:2],-5000,1000)}

plot_count_lik=plot_count_lik+1
lik_plots[[plot_count_lik]]<- plotVersusTrue(limits_lik
                                             ,loglik
                                             ,NA
                                             ,tau
                                             ,xlab
                                             ,ylab="LogLik"
                                             ,scaled=FALSE
                                             ,plotTrue = FALSE)
lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(limits_lik_skew
                                             ,loglik
                                             ,NA
                                             ,skew
                                             ,xlabskew
                                             ,ylab="LogLik"
                                             ,scaled=FALSE
                                             ,plotTrue = FALSE)
plot_count_lik=plot_count_lik+1
lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                             ,aic
                                             ,NA
                                             ,tau
                                             ,xlab
                                             ,ylab="AIC"
                                             ,scaled=FALSE
                                             ,plotTrue = FALSE)
lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                  ,aic
                                                  ,NA
                                                  ,skew
                                                  ,xlabskew
                                                  ,ylab="AIC"
                                                  ,scaled=FALSE
                                                  ,plotTrue = FALSE)

plot_count_lik=plot_count_lik+1
lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                             ,aic_4#bic
                                             ,NA
                                             ,tau
                                             ,xlab
                                             ,ylab="GAIC (4)"
                                             ,scaled=FALSE
                                             ,plotTrue = FALSE)

lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                  ,aic_4#bic
                                                  ,NA
                                                  ,skew
                                                  ,xlabskew
                                                  ,ylab="GAIC (4)"
                                                  ,scaled=FALSE
                                                  ,plotTrue = FALSE)
plot_count_lik=plot_count_lik+1
lik_plots[[plot_count_lik]]<- plotVersusTrue(c(limits_lik[1:2],limits_lik[c(4,3)]*-2)
                                             ,adj_bic#bic
                                             ,NA
                                             ,tau
                                             ,xlab
                                             ,ylab="BIC"
                                             ,scaled=FALSE
                                             ,plotTrue = FALSE)

lik_plots_skew[[plot_count_lik]]<- plotVersusTrue(c(limits_lik_skew[1:2],limits_lik_skew[c(4,3)]*-2)
                                                  ,adj_bic#bic
                                                  ,NA
                                                  ,skew
                                                  ,xlabskew
                                                  ,ylab="BIC"
                                                  ,scaled=FALSE
                                                  ,plotTrue = FALSE)

###################PLOT FUNCTION#############################

#ggarrange(plotlist=bias_plots       ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("N","N","NB","NB","G","G"),hjust=-.01) + bgcolor("white") + border(color = "white") # Bias x Tau
#ggsave(file=paste("simulation_bias_B0_Bt",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=9,dpi=900)
#ggarrange(plotlist=c(bias_plots[c(1)],error_plots[c(1)],bias_plots[c(3)],error_plots[c(3)],bias_plots[c(5)],error_plots[c(5)],bias_plots[c(7)],error_plots[c(7)])       ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("N","N","NB","NB","G","G","LO","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Tau
#ggsave(file=paste("Charts/simulation_bias_plus_error_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=11,dpi=900)

ggarrange(plotlist=c(bias_plots[c(2)],error_plots[c(2)],bias_plots[c(4)],error_plots[c(4)],bias_plots[c(6)],error_plots[c(6)],bias_plots[c(8)],error_plots[c(8)])       ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("N","N","NB","NB","G","G","LO","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Tau
ggsave(file=paste("Charts/simulation_bias_plus_error_Bt",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=11,dpi=900)



#ggarrange(plotlist=c(skew_bias_plots[c(1)],skew_bias_plots[c(3)],skew_bias_plots[c(5)])  ,common.legend=TRUE, ncol=3, nrow=1,      labels=c("NB","G","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Skew#ggsave(file=paste("simulation_bias_skew_AIO_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=3.5,dpi=900)
#ggsave(file=paste("Charts/simulation_bias_by_skew_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=4,dpi=900)
#ggarrange(plotlist=error_plots      ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("N","N","NB","NB","G","G"),hjust=-0.01) + bgcolor("white") + border(color = "white") # Error x Tau
#ggsave(file=paste("simulation_error_B0_Bt",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=9,dpi=900)
#ggsave(file=paste("simulation_error_AIO_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=9,dpi=900)
#ggarrange(plotlist=skew_error_plots ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("P","P","G","G")) + bgcolor("white") + border(color = "white") # Error x Skew
#ggsave(file=paste("simulation_error_skew_AIO_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=6,dpi=900)
#ggarrange(plotlist=lik_plots[c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]     ,common.legend=TRUE, ncol=4, nrow=plot_count_lik/4, labels=c("N","NB","G","LO","N","NB","G","LO","N","NB","G","LO","N","NB","G","LO"),hjust=0.1,font.label = list(size = 12)) + bgcolor("white") + border(color = "white") # Likelihoods x Tau
#ggsave(file=paste("Charts/simulation_loglik_AIO_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=15,dpi=900)
#ggarrange(plotlist=lik_plots_skew   ,common.legend=TRUE, ncol=4, nrow=plot_count_lik/4, labels=c("P","P","P","G","G","G")) + bgcolor("white") + border(color = "white") # Likelihoods x Skew
#ggsave(file=paste("simulation_loglik_AIO_skew_",parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=6,dpi=900)

#One at a time
#B1/B2 mode
ggarrange(plotlist=c(bias_plots[c(1)],error_plots[c(1)])       ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("LO","LO"),hjust=-.1) + bgcolor("white") + border(color = "white") # Bias x Tau
ggsave(file=paste("simulation_bias_plus_error_SINGLE_",dist,parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=3.45,dpi=900)

#B0/Bt mode
ggarrange(plotlist=bias_plots       ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("LO","LO"),hjust=-.01) + bgcolor("white") + border(color = "white") # Bias x Tau
ggsave(file=paste("Charts/simulation_bias_B0_Bt_SINGLE_",dist,parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=3.45,dpi=900)
ggarrange(plotlist=error_plots      ,common.legend=TRUE, ncol=2, nrow=plotcount/2,      labels=c("LO","LO"),hjust=-0.01) + bgcolor("white") + border(color = "white") # Error x Tau
ggsave(file=paste("Charts/simulation_error_B0_Bt_SINGLE_",dist,parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=8,height=3.45,dpi=900)

ggarrange(plotlist=lik_plots,common.legend=TRUE, ncol=4, nrow=1, labels=c("LO","LO","LO","LO"),hjust=0.1,font.label = list(size = 12)) + bgcolor("white") + border(color = "white") # Likelihoods x Tau
ggsave(file=paste("Charts/simulation_loglik_AIO_SINGLE_",dist,parameters[1,"n"],"_",Sys.Date(),".png",sep=""),last_plot(),width=9,height=3,dpi=900)


###guides(color=guide_legend(override.aes=list(fill=NA)))

#############Bias v skew table

bias_all[[plotcount/2]]<-cbind(trunc(tau*10,1)*10,trunc(skew),(if(dist=="NO"){t1intercepts}else{exp(t1intercepts)}/mu1)-1)

#bias_all_combined<-rbind(bias_all[[1]],bias_all[[2]])
#bias_all<-bias_all_combined

bias_table_list<-list()
for (i in 1:3) {
  dataset<- bias_all[[plotcount/2]][,c(1,2,4+i)]
  colnames(dataset)<-c("tau","skew","bias")
  summary<-aggregate(dataset[,"bias"] ~ dataset[,"tau"] + dataset[,"skew"], data = dataset, mean, na.rm = TRUE)
  colnames(summary)<-c("tau","skew","bias")
  bias_table_list[[i]]<-round(xtabs(summary[,"bias"] ~ summary[,"tau"] + summary[,"skew"]),2)
}

