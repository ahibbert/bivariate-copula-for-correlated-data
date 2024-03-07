###########CHOOSE DISTRIBUTION AND PARAMETERS######################################

#Load a simulation
#if(dist=="NO"){load("results_NO_n1000_1_2.rds")}
#if(dist=="GA"){load("results_combined_N_C0_n1000_geefix_mu1mu21012_GEEFIXV2.rds")}
#if(dist=="PO"){load("results_combined_PO_1000_2024-02-19")}

load("results_combined_GA_1000_2024-03-06.RData") ####GA
#load("results_combined_NO_1000_2024-03-06.RData") ####NO
#load("results_combined_PO_1000_2024-03-07.RData") ####NB/PO

results <- results_combined; dist="GA"

#Take out parameters
parameters=matrix(0,nrow=length(results_combined),ncol=6)
for (z in 1:length(results_combined)) {
  parameters[z,]=results_combined[[z]][nrow(results_combined[[z]]),1:6]
}
colnames(parameters)<-c("n","a","b","c","mu1","mu2")

#################################1. DATA SETUP##################################################

require(latex2exp)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)

set.seed(1000)
options(scipen=999)

#Theoretical errors
if(dist=="GA") {
  #Parameters
  mu1=parameters[,"a"]*parameters[,"mu1"]
  mu2=parameters[,"a"]*parameters[,"mu2"]
  #Errors
  load(file="numDerivResults_20231127.rds")
  trueSE<-numDerivResults[,c(1,2,5)]
}

if(dist=="NO") {
  #Parameters
  mu1=parameters[,"mu1"]
  mu2=parameters[,"mu2"]
  #Errors
  trueSE<-t(rbind((parameters[,"a"]*sqrt(1-(parameters[,"c"]^2)))/sqrt(n)
                  ,(parameters[,"b"]*sqrt(1-(parameters[,"c"]^2)))/sqrt(n)
                  ,sqrt((parameters[,"a"]^2)+(parameters[,"b"]^2)-2*parameters[,"a"]*parameters[,"b"]*parameters[,"c"])/sqrt(n)))
  colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")
}
if(dist=="PO") {
  #Parameters
  mu1=parameters[,"mu1"]*parameters[,"c"]
  mu2=parameters[,"mu2"]*parameters[,"c"]
  #Errors
  trueSE<-matrix(ncol=3,nrow=length(results))
  for (i in 1:length(results)) {
    trueSE[i,]<-results[[i]]["actuals",c("se_b1","se_b2","se_b2")]  
  }
  colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu2_se_Bt")
}
  
tau<-emp_cor<-vector(length=length(results))
t1intercepts<-t1error<-t2intercepts<-t2error<-aic<-bic<-loglik<-matrix(ncol=(nrow(results[[1]])-2),nrow=length(results))
for (i in 1:length(results)) {
  t1intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_1"])
  t2intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_2"])
  t1error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b1"])
  t2error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b2"])
  
  loglik[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"LogLik"])
  aic[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"AIC"])
  bic[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"BIC"])
  
  tau[i]=results[[i]][(nrow(results[[i]])-1),6]/100
  emp_cor[i]=results[[i]][(nrow(results[[i]])-1),7]/100
    
  colnames(t1intercepts)<-colnames(t2intercepts)<-colnames(t1error)<-colnames(t2error)<-colnames(loglik)<-colnames(aic)<-colnames(bic)<-colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))
}

###Wrapper for easy plotting
plotVersusTrue <- function (limits,inputs,true,tau,ylab,scaled=FALSE,type="ALL") {
  
  if (scaled==TRUE) {
    inputs = inputs / true-1
    true = true/true-1
  }
  
  inputs=as.data.frame(inputs)
  
  plot<-ggplot() + ylim(limits[3],limits[4]) + labs(x = TeX("Kendall's $\\tau$"), y=ylab) +
    {if(type=="ALL"||type=="non-GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_glm, color="GLM"),linetype = "dashed",se=FALSE)}} + 
    {if(type=="ALL"||type=="non-GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_gee, color="GEE"),linetype = "dashed",se=FALSE,position=position_jitter(w=0.005, h=0.0001))}} +
    {if(type=="ALL"||type=="non-GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_lme4, color="LME4"),linetype = "dashed",se=FALSE)}} + 
    {if(type=="ALL"||type=="non-GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_re_nosig, color="GAMLSS (4)"),linetype = "dashed",se=FALSE)}} +
    {if(type=="ALL"||type=="non-GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_re_np, color="GAMLSS NP (5)"),linetype = "dashed",se=FALSE)}} +
    {if(type=="ALL"||type=="GJRM") {geom_smooth(data=inputs, aes(x=tau, y=summary_cop, color="GJRM (C)"),linetype = "dashed",se=FALSE)}} +
    {if(type=="ALL"||type=="GJRM") { geom_smooth(data=inputs, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),linetype = "dashed",se=FALSE)}} +
    geom_smooth(data=inputs, aes(x=tau, y=true, color="True"), span=1,se=FALSE) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"#000000"))
  return(plot)
}

###################### BIAS CHARTS ######################

bias_1_plot<- plotVersusTrue(c(-1,1,-1,1)
               ,if(dist=="NO"){t1intercepts}else{exp(t1intercepts)}
               ,mu1
               ,tau
               ,ylab=TeX("$(\\hat{\\mu_1}/\\mu_1)-1$")
               , scaled=TRUE)
bias_2_plot<- plotVersusTrue(c(-1,1,-1,1)
               ,if(dist=="NO"){cbind((t2intercepts+t1intercepts)[,1:5],(t2intercepts)[,6:ncol(t2intercepts)])}
                else{cbind(exp(t2intercepts+t1intercepts)[,1:5],exp(t2intercepts)[,6:ncol(t1intercepts)])}
               ,mu2
               ,tau
               ,ylab=TeX("$(\\hat{\\mu_2}/\\mu_2)-1$")
               , scaled=TRUE)
bias_3_plot<- plotVersusTrue(c(-1,1,-1,1)
              ,if(dist=="NO"){cbind((t2intercepts)[,1:5],(t2intercepts-t1intercepts)[,6:ncol(t1intercepts)])}
               else{cbind(exp(t2intercepts)[,1:5],exp(t2intercepts-t1intercepts)[,6:ncol(t1intercepts)])}
              ,if(dist=="NO"){(mu2-mu1)}else{(mu2/mu1)}
              ,tau
              ,ylab=TeX("$\\left(\\frac{\\hat{\\mu_2}}{\\hat{\\mu_1}}\\div\\frac{\\mu_2}{\\mu_1}\\right)-1$")
              ,scaled=TRUE)

plot.new()
  ggarrange(bias_1_plot,bias_2_plot, bias_3_plot,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")  + guides(color=guide_legend(override.aes=list(fill=NA)))

ggsave(file=paste("simulation_bias_AIO_",dist,"_",n,"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=3,dpi=900)
  
###################### ERROR CHARTS #####################

error_1_plot<- plotVersusTrue(c(-1,1,0,.08)
                             ,t1error
                             ,trueSE[,"mu1_se"]
                             ,tau
                             ,ylab=TeX("$SE(\\hat{\\beta_{1}})$")
                             ,scaled=FALSE)
error_2_plot<- plotVersusTrue(c(-1,1,0,.08)
                              ,t1error
                              ,trueSE[,"mu2_se_B2"]
                              ,tau
                              ,ylab=TeX("$SE(\\hat{\\beta_{2}})$")
                              ,scaled=FALSE
                              ,type="GJRM")
error_2_plot_bt<- plotVersusTrue(c(-1,1,0,.08)
                              ,t1error
                              ,trueSE[,"mu2_se_Bt"]
                              ,tau
                              ,ylab=TeX("$SE(\\hat{\\beta_{t}})$")
                              ,scaled=FALSE
                              ,type="non-GJRM")

ggarrange(error_1_plot,error_2_plot,error_2_plot_bt,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")
  
ggsave(file=paste("simulation_error_AIO_",dist,"_",n,"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=3,dpi=900)








  ##########Individual plots###############
  par(mfrow=c(2,3))
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_glm"]      ,xlab="True SE",ylab="Estimated SE", main="GLM")      ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_gee"]      ,xlab="True SE",ylab="Estimated SE", main="GEE")      ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_re_nosig"] ,xlab="True SE",ylab="Estimated SE", main="GLMM (4)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_re_np"]       ,xlab="True SE",ylab="Estimated SE", main="GLMM (5)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_cop"]      ,xlab="True SE",ylab="Estimated SE", main="GJRM (C)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu1_se"],y=t1error[,"summary_cop_n"]    ,xlab="True SE",ylab="Estimated SE", main="GJRM (N)") ; abline(a=0, b=1,col="red")
  
  par(mfrow=c(2,3))
  plot(x=trueSE[,"mu2_se_Bt"],y=t2error[,"summary_glm"]      ,xlab="True SE",ylab="Estimated SE", main="GLM")      ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu2_se_Bt"],y=t2error[,"summary_gee"]      ,xlab="True SE",ylab="Estimated SE", main="GEE")      ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu2_se_Bt"],y=t2error[,"summary_re_nosig"] ,xlab="True SE",ylab="Estimated SE", main="GLMM (4)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu2_se_Bt"],y=t2error[,"summary_re"]       ,xlab="True SE",ylab="Estimated SE", main="GLMM (5)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu2_se_B2"],y=t2error[,"summary_cop"]      ,xlab="True SE",ylab="Estimated SE", main="GJRM (C)") ; abline(a=0, b=1,col="red")
  plot(x=trueSE[,"mu2_se_B2"],y=t2error[,"summary_cop_n"]    ,xlab="True SE",ylab="Estimated SE", main="GJRM (N)") ; abline(a=0, b=1,col="red")
  
  ##############TESTING
  plot.new()  
  par(mfrow=c(3,2))
  
  scatter.smooth(t1error[,1]~tau,main="Error GLM",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  scatter.smooth(t1error[,2]~tau,main="Error GEE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  
  scatter.smooth(t1error[,3]~tau,main="Error RE no sig",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  scatter.smooth(t1error[t1error[,4]>.0001,4]~tau[t1error[,4]>.0001],main="Error RE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  
  scatter.smooth(t1error[,5]~tau,main="Error Copula Clayton",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  scatter.smooth(t1error[,6]~tau,main="Error Copula Normal",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  lines(lowess(t1error[,1]~tau),col="red")
  
  plot(t1intercepts[,1]/log(t1intercepts[,7]/mu1)-1~tau,main="Bias GLM",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  plot(t1intercepts[,2]/log(t1intercepts[,7]/mu1)-1~tau,main="Bias GEE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  
  plot(t1intercepts[,3]/log(t1intercepts[,7]/mu1)-1~tau,main="Bias RE no sig",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  plot(t1intercepts[t1error[,4]>.0001,4]/log(t1intercepts[t1error[,4]>.0001,7]/mu1)-1~tau[t1error[,4]>.0001],main="Bias RE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  
  plot(t1intercepts[,5]/log(t1intercepts[,7]/mu1)-1~tau,main="Bias Copula Clayton",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  plot(t1intercepts[,6]/log(t1intercepts[,7]/mu1)-1~tau,main="Bias Copula Normal",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  abline(h=0,col="red")
  
  ##############BOXPLOTS
  plot.new()
  par(mfrow=c(3,2))
  plot(t1intercepts[,1]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GLM"); abline(h=0,col="blue")
  plot(t1intercepts[,2]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GEE"); abline(h=0,col="blue")
  plot(t1intercepts[,3]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GLMM no sig"); abline(h=0,col="blue")
  plot(t1intercepts[,4]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GLMM"); abline(h=0,col="blue")
  plot(t1intercepts[,5]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GJRM (Clayton)"); abline(h=0,col="blue")
  plot(t1intercepts[,6]/log(t1intercepts[,7]/mu1)-1~as.factor(round(tau*10)/10),ylim=c(-5,5),xlab="Tau", ylab="Bias",main="Bias GJRM (Gaussian)"); abline(h=0,col="blue")
  
  plot.new()
  par(mfrow=c(3,2))
  plot(t1error[,1]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLM"); abline(h=0,col="blue")
  plot(t1error[,2]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GEE"); abline(h=0,col="blue")
  plot(t1error[,3]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLMM no sig"); abline(h=0,col="blue")
  plot(t1error[t1error[,4]>0.001,4]~as.factor(round(tau[t1error[,4]>0.001]*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLMM"); abline(h=0,col="blue")
  plot(t1error[,5]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GJRM (Clayton)"); abline(h=0,col="blue")
  plot(t1error[,6]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GJRM (Gaussian)"); abline(h=0,col="blue")
  
  plot.new()
  par(mfrow=c(3,2))
  plot(t2intercepts[,1]/log((t1intercepts[,7]/mu1)/(t1intercepts[,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GLM",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  plot(t2intercepts[,2]/log((t1intercepts[,7]/mu1)/(t1intercepts[,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GEE",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  plot(t2intercepts[,3]/log((t1intercepts[,7]/mu1)/(t1intercepts[,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GLMM no sig",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  plot(t2intercepts[t2error[,4]>0.001,4]/log((t1intercepts[t2error[,4]>0.001,7]/mu1)/(t1intercepts[t2error[,4]>0.001,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GLMM",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  plot(t2intercepts[,5]/log((t1intercepts[,7]/mu1)/(t1intercepts[,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GJRM (Copula)",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  plot(t2intercepts[,6]/log((t1intercepts[,7]/mu1)/(t1intercepts[,7]/mu2))+1~as.factor(round(tau*10)/10),main="Bias GJRM (Gaussian)",ylim=c(-4,4),xlab="Tau", ylab="Bias"); abline(h=0,col="blue")
  
  plot.new()
  par(mfrow=c(3,2))
  plot(t2error[,1]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLM"); abline(h=0,col="blue")
  plot(t2error[,2]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GEE"); abline(h=0,col="blue")
  plot(t2error[,3]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLMM no sig"); abline(h=0,col="blue")
  plot(t2error[t2error[,4]>0.001,4]~as.factor(round(tau[t2error[,4]>0.001]*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GLMM"); abline(h=0,col="blue")
  plot(t2error[,5]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GJRM (Clayton)"); abline(h=0,col="blue")
  plot(t2error[,6]~as.factor(round(tau*10)/10),ylim=c(0,.3),xlab="Tau", ylab="Standard Error",main="Standard Error GJRM (Gaussian)"); abline(h=0,col="blue")
  
  