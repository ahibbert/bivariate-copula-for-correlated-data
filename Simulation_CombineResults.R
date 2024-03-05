###########CHOOSE DISTRIBUTION AND PARAMETERS######################################

#Combine simulations

#if(dist=="NO"){load("results_NO_n1000_1_2.rds")}
#if(dist=="GA"){load("results_combined_N_C0_n1000_geefix_mu1mu21012_GEEFIXV2.rds")}
#if(dist=="PO"){load("results_combined_PO_1000_2024-02-19")}

results=results_combined

###HOPE THIS ONE WORKS
parameters=matrix(0,nrow=length(results),ncol=8)
for (i in length(results)) {
  parameters[z,]=results[[z]][nrow(results[[z]]),]
}







#a=.1+.1*1:20; b=.1+.1*1:20; mu1=10; mu2=12; n=1000;dist="GA"
#a=.5*1:5; b=.5*1:5;c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9); mu1=1; mu2=2; n=1000;dist="NO"
#a=NA; b=NA;c=c(.25,.5,1,5); mu1=c(.25,.5,1,5); mu2=c(.25,.5,1,5); n=1000;dist="PO"

parameters=matrix(0,nrow=length(a)*length(b)*length(mu1)*length(mu2)*length(c),ncol=6)
z=1
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    for (k in 1:length(mu1)) {
      for (l in 1:length(mu2)) {
        for (m in 1:length(c))  {
          parameters[z,]= c(n,a[i],b[j],c[m],mu1[k],mu2[l])
          z=z+1
        }
      }
    }
  }
}
colnames(parameters)<-c("n","a","b","c","mu1","mu2")






#################################1. DATA SETUP##################################################

require(latex2exp)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)

set.seed(1000)
options(scipen=999)

if(dist=="NO"){tau=parameters[,"c"]}
if(dist=="GA"){}#UPDATE value of tau
if(dist=="PO"){pcol=((1+(1/parameters[,"mu1"]))*(1+(1/parameters[,"mu2"])))^-0.5;tau=pcol}

#Theoretical errors
if(dist=="GA") {
  #Parameters
  mu1=parameters["a"]/parameters[,"mu1"]
  mu2=parameters["a"]/parameters[,"mu2"]
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
t1intercepts=matrix(ncol=(nrow(results[[1]])-2),nrow=length(results))
t1error=t1intercepts;t2error=t1intercepts;t2intercepts=t1intercepts
for (i in 1:length(results)) {
  t1intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_1"])
  t2intercepts[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"b_2"])
  t1error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b1"])
  t2error[i,]=t(results[[i]][1:(nrow(results[[i]])-2),"se_b2"])
  
  tau[i]=results[[i]][(nrow(results[[i]])-1),5]
  emp_cor[i]=results[[i]][(nrow(results[[i]])-1),6]
    
  colnames(t1intercepts)=colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))
  colnames(t2intercepts)=colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))
  colnames(t1error)=colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))
  colnames(t2error)=colnames(t(results[[i]][1:(nrow(results[[i]])-2),"b_1"]))

}

###################### BIAS CHARTS ######################

  ylim_max<-1

  if(dist=="GA" || dist=="PO") {data_input<-as.data.frame(cbind(exp(t1intercepts),mu1,mu2,tau))}
  if(dist=="NO") {data_input<-as.data.frame(cbind((t1intercepts),mu1,mu2,tau))}

  bias_1_plot<-ggplot() + ylim(-1*ylim_max,1*ylim_max) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\hat{\\mu_1}/\\mu_1)-1$")) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/(mu1))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(mu1))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_lme4)/(mu1))-1, color="LME4"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(mu1))-1, color="GAMLSS (4)"),level=.99,position=position_jitter(w=0.005, h=0.0001)) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_np)/(mu1))-1, color="GAMLSS NP (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(mu1))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(mu1))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"black"))
  
  if(dist=="GA"||dist=="PO") {data_input<-as.data.frame(cbind(exp(t2intercepts+t1intercepts)[,1:5],exp(t2intercepts)[,6:ncol(t1intercepts)],mu1,mu2,tau))}
  if(dist=="NO") {data_input<-as.data.frame(cbind((t2intercepts+t1intercepts)[,1:5],(t2intercepts)[,6:ncol(t2intercepts)],mu1,mu2,tau))}

  bias_2_plot<-ggplot() + ylim(-1*ylim_max,1*ylim_max) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\hat{\\mu_2}/\\mu_2)-1$")) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/(mu2))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(mu2))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_lme4)/(mu2))-1, color="LME4"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(mu2))-1, color="GAMLSS (4)"),level=.99,position=position_jitter(w=0.005, h=0.0001)) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_np)/(mu2))-1, color="GAMLSS NP (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(mu2))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(mu2))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"black"))
  
  if(dist=="GA"||dist=="PO") {data_input<-as.data.frame(cbind(exp(t2intercepts)[,1:5],exp(t2intercepts-t1intercepts)[,6:ncol(t1intercepts)],mu1,mu2/mu1,tau))}
  if(dist=="NO") {data_input<-as.data.frame(cbind((t2intercepts)[,1:5],(t2intercepts-t1intercepts)[,6:ncol(t1intercepts)],mu1,mu2-mu1,tau))}
  bias_3_plot<-ggplot()  + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\frac{\\hat{\\mu_2}}{\\hat{\\mu_1}}\\div\\frac{\\mu_2}{\\mu_1})-1$")) + ylim(-1*ylim_max,1*ylim_max)+
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/((V17)))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(V17))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(V17))-1, color="GAMLSS (4)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_np)/(V17))-1, color="GAMLSS NP (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(V17))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(V17))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +

    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"black"))
    
plot.new()
  ggarrange(bias_1_plot,bias_2_plot, bias_3_plot,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")  + guides(color=guide_legend(override.aes=list(fill=NA)))

#ggsave(file="simulation_bias_charts_all_in_one.png",last_plot(),width=12,height=3,dpi=900)
ggsave(file=paste("simulation_bias_AIO_",dist,"_",n,"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=3,dpi=900)
  
###################### ERROR CHARTS #####################

  if(dist=="GA"||dist=="PO") {data_input<-as.data.frame(cbind(t1error,tau,trueSE))}
  if(dist=="NO") {data_input<-as.data.frame(cbind(t1error,tau,trueSE))} ###Add theoretircal errors
  
#  data_input<-as.data.frame(cbind(t1error,tau,numDerivResults[,c(1,2)]))
  error_1_plot<- ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{1}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),linetype = "dashed",se=FALSE) + 
    geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),linetype = "dashed",se=FALSE,position=position_jitter(w=0.005, h=0.0001)) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_lme4, color="LME4"),linetype = "dashed",se=FALSE) + 
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GAMLSS (4)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_np, color="GAMLSS NP (5)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=mu1_se, color="True"), span=1,se=FALSE) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"#000000"))

  #                    ) +
  #scale_linetype_manual(name="linetype", c('GLM'=1, 'GEE'=2, 'GLMM (4)'=3,'GLMM (5)'=4, 'GJRM (C)'=5, 'GJRM (N)'=6))
  #theme(legend.position = "right", legend.title=element_text(size=20),
  #      legend.text=element_text(size=14))
  
  #data_input<-as.data.frame(cbind(t2error,tau,numDerivResults[,c(1,2)]))
  
  if(dist=="GA"||dist=="PO") {data_input<-as.data.frame(cbind(t2error,tau,trueSE))}
  if(dist=="NO") {data_input<-as.data.frame(cbind(t2error,tau,trueSE))} ###Add theoretircal errors
  
  error_2_plot<-
    ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{2}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=mu2_se_B2, color="True"), span=1,level=.95,se=FALSE) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),level=.99) + 
    #geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GLMM (4)"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_re, color="GLMM (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),level=.95,se=FALSE,linetype = "dashed") +
    
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GLMM (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"#000000"))
    
  error_2_plot_bt<-
    ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{t}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=mu2_se_Bt, color="True"), span=1,level=.95,se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),level=.95,se=FALSE,linetype = "dashed",position=position_jitter(w=0.005, h=0.0001)) + 
    geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_lme4, color="LME4"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GAMLSS (4)"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_np, color="GAMLSS NP (5)"),level=.95,se=FALSE,linetype = "dashed") +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),level=.99) +
    
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","LME4","GAMLSS (4)","GAMLSS NP (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 7, name = "Dark2"),"#000000"))
  
  ggarrange(error_1_plot,error_2_plot,error_2_plot_bt,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")
  
  ggsave(file=paste("simulation_error_AIO_",dist,"_",n,"_",Sys.Date(),".png",sep=""),last_plot(),width=12,height=3,dpi=900)
  
  #ggsave(file="simulation_error_charts_all_in_one.png",last_plot(),width=12,height=3,dpi=900)
  #ggsave(file="simulation_error_charts_all_in_one_MVTNORMAL.png",last_plot(),width=12,height=3,dpi=900)
  
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
  
  