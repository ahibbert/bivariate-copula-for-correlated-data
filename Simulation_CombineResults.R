########01 DATA SETUP (DON'T EDIT) ###############

#Combine simulations
#load("results_NO_n1000_1_2.rds")
load("results_combined_N_C0_n1000_geefix_mu1mu21012_GEEFIXV2.rds")

require(latex2exp)
require(ggplot2)
require(ggpubr)
library(RColorBrewer)
set.seed(1000)
options(scipen=999)

a=.1+.1*1:20; b=.1+.1*1:20; mu1=10; mu2=12; n=1000;dist="GA"
#a=.5*1:5; b=.5*1:5;c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9); mu1=1; mu2=2; n=1000;dist="NO"

  #write.csv(cbind(parameters,tau,marginal_skew_1,marginal_skew_2, t1intercepts, t2intercepts,t1error,t2error),file="SimulationResults.csv")

  tau=results[[1]][7,3]
  for (i in 2:length(results)) {tau=rbind(tau,results[[i]][7,3]) }
  
  pcol=results[[1]][7,4]
  for (i in 2:length(results)) {pcol=rbind(pcol,results[[i]][7,4]) }
  
  parameters=results[[1]][8,1:4]
  for (i in 2:length(results)) {parameters=rbind(parameters,results[[i]][8,1:4]) }
  
  if(dist=="GA") {
    #For Gamma
    mu1=parameters[,1]/parameters[,4]
    mu2=parameters[,1]/parameters[,3]
  }
  if(dist=="NO") {
    mu1=parameters[,3]
    mu2=parameters[,4]
  }
  
  #Theoretical errors
  if(dist=="GA") {
    load(file="numDerivResults_20231127.rds")
    trueSE<-numDerivResults[,c(1,2,5)]
  }
  if(dist=="NO") {
    trueSE<-rbind(parameters[,1]/sqrt(n),parameters[,2]/sqrt(n),sqrt((parameters[,1]^2)+(parameters[,2]^2))/sqrt(n))
    colnames(trueSE)<-c("mu1_se","mu2_se_B2","mu_se_Bt")
  }
  
  t1error=results[[1]][1:6,3]
  for (i in 2:length(results)) {t1error=rbind(t1error,results[[i]][1:6,3]) }
  t2error=results[[1]][1:6,4]
  for (i in 2:length(results)) {t2error=rbind(t2error,results[[i]][1:6,4]) }
  
  #########BIAS
  t1intercepts=results[[1]][c(1:8),1]
  for (i in 2:length(results)) {t1intercepts=rbind(t1intercepts,results[[i]][c(1:8),1]) }
  t2intercepts=results[[1]][c(1:8),2]
  for (i in 2:length(results)) {t2intercepts=rbind(t2intercepts,results[[i]][c(1:8),2]) }

  ##########TABLE
  
  # par(mfrow=c(2,3))
  # i=1
  # a<-boxplot(t1intercepts[,i]/t1intercepts[,7]-1~as.factor(round(tau*5)/5))
  # z=cbind(i,a$stats,a$n)
  # for (i in 2:6) {
  #   a<-boxplot(t1intercepts[,i]/t1intercepts[,7]-1~as.factor(round(tau*5)/5))
  #   z=rbind(z,cbind(i,a$stats,a$n))
  # }
  # 
  # time1biastable <- z[c(3,3+5,3+5*2,3+5*3,3+5*4,3+5*5),]
  # 
  # i=1
  # a<-boxplot(t1error[,1]~as.factor(round(tau*5)/5))
  # z=cbind(i,a$stats,a$n)
  # for (i in 2:6) {
  #   a<-boxplot(t1error[,i]~as.factor(round(tau*5)/5))
  #   z=rbind(z,cbind(i,a$stats,a$n))
  # }
  # 
  # time1errortable <- z[c(3,3+5,3+5*2,3+5*3,3+5*4,3+5*5),]
  # 
  # i=1
  # a<-boxplot(t2intercepts[,i]/t2intercepts[,7]-1~as.factor(round(tau*5)/5))
  # z=cbind(i,a$stats,a$n)
  # for (i in 2:6) {
  #   a<-boxplot(t2intercepts[,i]/t2intercepts[,7]-1~as.factor(round(tau*5)/5))
  #   z=rbind(z,cbind(i,a$stats,a$n))
  # }
  # 
  # time2biastable <- z[c(3,3+5,3+5*2,3+5*3,3+5*4,3+5*5),]
  # 
  # i=1
  # a<-boxplot(t2error[,1]~as.factor(round(tau*5)/5))
  # z=cbind(i,a$stats,a$n)
  # for (i in 2:6) {
  #   a<-boxplot(t2error[,i]~as.factor(round(tau*5)/5))
  #   z=rbind(z,cbind(i,a$stats,a$n))
  # }
  # 
  # time2errortable <- z[c(3,3+5,3+5*2,3+5*3,3+5*4,3+5*5),]
  # 
  # summaryresultstable <- rbind(time1biastable,time1errortable,time2biastable,time2errortable)
  # colnames(summaryresultstable) <- c("model","0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0","n")
  # 
  #write.csv(summaryresultstable,file="simulation_full_results_table_n1000_geefix.csv")
  #write.csv(cbind(t1intercepts,tau,marginal_skew_1,marginal_skew_2),file="simulation_skewness_tables_n1000_geefix.csv")
  

###################### BIAS CHARTS ######################
  
  if(dist=="GA") {data_input<-as.data.frame(cbind(exp(t1intercepts[,1:6]),mu1,tau))}
  if(dist=="NO") {data_input<-as.data.frame(cbind((t1intercepts[,1:6]),mu1,tau))}
  
  colnames(data_input)[7:8] <- c("actuals","tau")
  bias_1_plot<-ggplot() + ylim(-1,1) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\hat{\\mu_1}/\\mu_1)-1$")) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/(actuals))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(actuals))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(actuals))-1, color="GLMM (4)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re)/(actuals))-1, color="GLMM (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(actuals))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(actuals))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"black"))

  data_input<-as.data.frame(cbind(exp(t2intercepts[,1:6]+t1intercepts[,1:6]),mu2,tau))
  colnames(data_input)[7:8] <- c("actuals","tau")
  bias_2_plot<-ggplot() + ylim(-1,1) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\hat{\\mu_2}/\\mu_2)-1$")) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/(actuals))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(actuals))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(actuals))-1, color="GLMM (4)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re)/(actuals))-1, color="GLMM (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(actuals))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(actuals))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"black"))
  
  data_input<-as.data.frame(cbind(exp(t2intercepts[,1:6]),mu2/mu1,tau))
  colnames(data_input)[7:8] <- c("actuals","tau")
  bias_3_plot<-ggplot()  + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$(\\frac{\\hat{\\mu_2}}{\\hat{\\mu_1}}\\div\\frac{\\mu_2}{\\mu_1})-1$")) + ylim(-1,1)+
    geom_smooth(data=data_input, aes(x=tau, y=((summary_glm)/(actuals))-1, color="GLM"),level=.99) + 
    geom_smooth(data=data_input, aes(x=tau, y=((summary_gee)/(actuals))-1, color="GEE"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re_nosig)/(actuals))-1, color="GLMM (4)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_re)/(actuals))-1, color="GLMM (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop)/(actuals))-1, color="GJRM (C)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=((summary_cop_n)/(actuals))-1, color="GJRM (N)"),level=.99) +
    geom_hline(aes(yintercept=0,color="True"), size=1) +
    guides(color=guide_legend(override.aes=list(fill=NA))) +

    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"black"))
    
plot.new()
  ggarrange(bias_1_plot,bias_2_plot, bias_3_plot,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")  + guides(color=guide_legend(override.aes=list(fill=NA)))

  
#ggsave(file="simulation_bias_charts_all_in_one.png",last_plot(),width=12,height=3,dpi=900)
  
###################### ERROR CHARTS #####################

  if(dist=="GA") {data_input<-as.data.frame(cbind(t1error,tau,trueSE))}
  if(dist=="NO") {data_input<-as.data.frame(cbind(t1error,tau,trueSE))} ###Add theoretircal errors
  
#  data_input<-as.data.frame(cbind(t1error,tau,numDerivResults[,c(1,2)]))
  error_1_plot<- ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{1}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),linetype = "dashed",se=FALSE) + 
    geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),linetype = "dashed",se=FALSE,position=position_jitter(w=0.005, h=0.0001)) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GLMM (4)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re, color="GLMM (5)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),linetype = "dashed",se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=mu1_se, color="True"), span=1,se=FALSE) +
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"#000000"))

  #                    ) +
  #scale_linetype_manual(name="linetype", c('GLM'=1, 'GEE'=2, 'GLMM (4)'=3,'GLMM (5)'=4, 'GJRM (C)'=5, 'GJRM (N)'=6))
  #theme(legend.position = "right", legend.title=element_text(size=20),
  #      legend.text=element_text(size=14))
  
  #data_input<-as.data.frame(cbind(t2error,tau,numDerivResults[,c(1,2)]))
  error_2_plot<-
    ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{2}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=mu2_se_B2, color="True"), span=1,level=.95,se=FALSE) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),level=.99) + 
    #geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GLMM (4)"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_re, color="GLMM (5)"),level=.99) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),level=.95,se=FALSE,linetype = "dashed") +
    
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"#000000"))
    
  error_2_plot_bt<-
    ggplot() + ylim(0,.08) + xlim(0.15,.75) + labs(x = TeX("Kendall's $\\tau$"), y=TeX("$SE(\\hat{\\beta_{t}})$")) +
    geom_smooth(data=data_input, aes(x=tau, y=mu2_se_Bt, color="True"), span=1,level=.95,se=FALSE) +
    geom_smooth(data=data_input, aes(x=tau, y=summary_glm, color="GLM"),level=.95,se=FALSE,linetype = "dashed") + 
    geom_smooth(data=data_input, aes(x=tau, y=summary_gee, color="GEE"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re_nosig, color="GLMM (4)"),level=.95,se=FALSE,linetype = "dashed") +
    geom_smooth(data=data_input, aes(x=tau, y=summary_re, color="GLMM (5)"),level=.95,se=FALSE,linetype = "dashed") +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_cop, color="GJRM (C)"),level=.99) +
    #geom_smooth(data=data_input, aes(x=tau, y=summary_cop_n, color="GJRM (N)"),level=.99) +
    
    scale_colour_manual(name="Model", breaks=c("GLM","GEE","GLMM (4)","GLMM (5)","GJRM (C)","GJRM (N)","True")
                        , values=c(brewer.pal(n = 6, name = "Dark2"),"#000000"))
  
  ggarrange(error_1_plot,error_2_plot,error_2_plot_bt,common.legend=TRUE,nrow=1, ncol=3, legend="right",labels="AUTO") + #,labels=c("(a)","(b)","(c)","(d)"), font.label = list(size=12,face="plain"
    bgcolor("white")+border(color = "white")
  
  #ggsave(file="simulation_error_charts_all_in_one.png",last_plot(),width=12,height=3,dpi=900)
  
  
#########Time 1 BIAS AND ERROR AGAINST TAU (OLD - separate plots) ##############
  
    #Bias time 1 (paramater)
    data_input<-as.data.frame(cbind(t1intercepts[,1:6],log(mu2),tau))
    colnames(data_input)[7:8] <- c("actuals","tau")
    a<-ggplot(data=data_input, aes(x=tau, y=summary_glm/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glm")
    b<-ggplot(data=data_input, aes(x=tau, y=summary_gee/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gee")
    c<-ggplot(data=data_input, aes(x=tau, y=summary_re_nosig/actuals-1))  + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glmm (4)")
    d<-ggplot(data=data_input, aes(x=tau, y=summary_re/actuals-1))        + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glmm (5)")
    e<-ggplot(data=data_input, aes(x=tau, y=summary_cop/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gjrm (C)")
    f<-ggplot(data=data_input, aes(x=tau, y=summary_cop_n/actuals-1))     + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gjrm (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="bias_time_1_par.jpeg",last_plot(),width=10,height=6,dpi=300)
    
    #Bias time 2 (paramater)
    data_input<-as.data.frame(cbind(t2intercepts[,1:6],log(mu1)-log(mu2),tau))
    colnames(data_input)[7:8] <- c("actuals","tau")
    a<-ggplot(data=data_input, aes(x=tau, y=summary_glm/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glm")
    b<-ggplot(data=data_input, aes(x=tau, y=summary_gee/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gee")
    c<-ggplot(data=data_input, aes(x=tau, y=summary_re_nosig/actuals-1))  + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glmm (4)")
    d<-ggplot(data=data_input, aes(x=tau, y=summary_re/actuals-1))        + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="glmm (5)")
    e<-ggplot(data=data_input, aes(x=tau, y=summary_cop/actuals-1))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gjrm (C)")
    f<-ggplot(data=data_input, aes(x=tau, y=summary_cop_n/actuals-1))     + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(-4,4) + labs(x = "tau", y="bias", title="gjrm (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="bias_time_2_par.jpeg",last_plot(),width=10,height=6,dpi=300)
    
    
    #Bias time 1 (MEAN)
    data_input<-as.data.frame(cbind(t1intercepts[,1:6],mu2,tau))
    colnames(data_input)[7:8] <- c("actuals","tau")
    a<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_glm)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLM")
    b<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_gee)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GEE")
    c<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_re_nosig)/(actuals))-1))  + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLMM (4)")
    d<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_re)/(actuals))-1))        + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLMM (5)")
    e<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_cop)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GJRM (C)")
    f<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_cop_n)/(actuals))-1))     + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GJRM (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="bias_time_1_mean.jpeg",last_plot(),width=10,height=6,dpi=300)
    
    #Bias time 2 (MEAN)
    data_input<-as.data.frame(cbind(t2intercepts[,1:6]+t1intercepts[,1:6],mu1,tau))
    colnames(data_input)[7:8] <- c("actuals","tau")
    a<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_glm)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLM")
    b<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_gee)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GEE")
    c<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_re_nosig)/(actuals))-1))  + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLMM (4)")
    d<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_re)/(actuals))-1))        + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GLMM (5)")
    e<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_cop)/(actuals))-1))       + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GJRM (C)")
    f<-ggplot(data=data_input, aes(x=tau, y=(exp(summary_cop_n)/(actuals))-1))     + geom_point(size=0.5,color="gray") + geom_smooth(level=0.99) + ylim(-1,1) + labs(x = "tau", y="bias", title="GJRM (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="bias_time_2_mean.jpeg",last_plot(),width=10,height=6,dpi=300)
    
    #Error time 1
    a<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_glm))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLM")
    b<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_gee))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GEE")
    c<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_re_nosig))  + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLMM (4)")
    d<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_re))        + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLMM (5)")
    e<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_cop))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GJRM (C)")
    f<-ggplot(data=as.data.frame(cbind(t1error,tau)), aes(x=tau, y=summary_cop_n))     + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GJRM (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="error_time_1.jpeg",last_plot(),width=10,height=6,dpi=300)
    

                          
    #Error time 2
    a<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_glm))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLM")
    b<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_gee))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GEE")
    c<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_re_nosig))  + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLMM (4)")
    d<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_re))        + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GLMM (5)")
    e<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_cop))       + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GJRM (C)")
    f<-ggplot(data=as.data.frame(cbind(t2error,tau)), aes(x=tau, y=summary_cop_n))     + geom_point(size=0.5,color="gray") + geom_smooth() + ylim(0,0.1) + labs(x = "tau", y="standard error", title="GJRM (N)")
    ggarrange(a,b,c,d,e,f)
    #ggsave(file="error_time_2.jpeg",last_plot(),width=10,height=6,dpi=300)
    
    ###Combined Error time 1 / time 2
    
   
      #                    , values=c('GLM'='red', 'GEE'='pink', 'GLMM (4)'='blue','GLMM (5)'='purple', 'GJRM (C)'='orange', 'GJRM (N)'='yellow')
      #                    ) +
      #scale_linetype_manual(name="linetype", c('GLM'=1, 'GEE'=2, 'GLMM (4)'=3,'GLMM (5)'=4, 'GJRM (C)'=5, 'GJRM (N)'=6))
      #theme(legend.position = "right", legend.title=element_text(size=20),
      #      legend.text=element_text(size=14))
      
    ggarrange(error_1_plot,error_2_plot,common.legend = TRUE)
    
########## REFERENCE
  #     
  #   plot(t1intercepts[,1]/t1intercepts[,7]-1~tau,main="Bias GLM",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   plot(t1intercepts[,2]/t1intercepts[,7]-1~tau,main="Bias GEE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   scatter.smooth(t1error[,1]~tau,main="Error GLM",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #   scatter.smooth(t1error[,2]~tau,main="Error GEE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #   
  #   plot(t1intercepts[,3]/t1intercepts[,7]-1~tau,main="Bias RE no sig",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   plot(t1intercepts[t1error[,4]>.0001,4]/t1intercepts[t1error[,4]>.0001,7]-1~tau[t1error[,4]>.0001],main="Bias RE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   scatter.smooth(t1error[,3]~tau,main="Error RE no sig",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #   scatter.smooth(t1error[t1error[,4]>.0001,4]~tau[t1error[,4]>.0001],main="Error RE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #   
  #   plot(t1intercepts[,5]/t1intercepts[,7]-1~tau,main="Bias Copula Clayton",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   plot(t1intercepts[,6]/t1intercepts[,7]/mu1-1~tau, main="Bias Copula Normal",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  #   abline(h=0,col="red")
  #   scatter.smooth(t1error[,5]~tau,main="Error Copula Clayton",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #   scatter.smooth(t1error[,6]~tau,main="Error Copula Normal",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  #   lines(lowess(t1error[,1]~tau),col="red")
  #  
  #  #########Time 2 BIAS AND ERROR AGAINST TAU
  # plot.new()  
  # 
  # par(mfrow=c(3,4))
  # plot(t2intercepts[,1]/t2intercepts[,7]-1~tau,main="Bias GLM",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # plot(t2intercepts[,2]/t2intercepts[,7]-1~tau,main="Bias GEE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # scatter.smooth(t2error[,1]~tau,main="Error GLM",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # scatter.smooth(t2error[,2]~tau,main="Error GEE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # 
  # plot(t2intercepts[,3]/t2intercepts[,7]-1~tau,main="Bias RE no sig",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # plot(t2intercepts[t2error[,4]>.0001,4]/t2intercepts[t2error[,4]>.0001,7]-1~tau[t2error[,4]>.0001],main="Bias RE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # scatter.smooth(t2error[,3]~tau,main="Error RE no sig",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # scatter.smooth(t2error[t2error[,4]>.0001,4]~tau[t2error[,4]>.0001],main="Error RE",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # 
  # plot(t2intercepts[,5]/t2intercepts[,7]-1~tau,main="Bias Copula Clayton",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # plot(t2intercepts[,6]/t2intercepts[,7]-1~tau,main="Bias Copula Normal",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # abline(h=0,col="red")
  # scatter.smooth(t2error[,5]~tau,main="Error Copula Clayton",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # scatter.smooth(t2error[,6]~tau,main="Error Copula Normal",ylim=c(0,.3),xlab="Tau", ylab="Error",col="gray")
  # lines(lowess(t2error[,1]~tau),col="red")
  # 

  
###############################################REFERENCE##################################
  
  # 
  # 
  #   
  # #########Time 1 Error AGAINST TAU
  # 
  # plot.new()  
  # par(mfrow=c(3,4))
  # 
  # ######Time 2 BIAS AGAINST TAU (should just be zero) 
  # plot(t2intercepts[,1]~tau,main="Bias GLM",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2intercepts[,2]~tau,main="Bias GEE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2error[,1]~tau,main="Error GLM",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # plot(t2error[,2]~tau,main="Error GEE",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # 
  # plot(t2intercepts[,3]~tau,main="Bias RE no sig",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2intercepts[,4]~tau,main="Bias RE",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2error[,3]~tau,main="Error RE no sig",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # plot(t2error[,4]~tau,main="Error RE",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # 
  # plot(t2intercepts[,5]~tau,main="Bias Copula Clayton",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2intercepts[,6]~tau,main="Bias Copula Normal",ylim=c(-4,4),xlab="Tau", ylab="Bias")
  # plot(t2error[,5]~tau,main="Error Copula Clayton",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # plot(t2error[,6]~tau,main="Error Copula Normal",ylim=c(0,.4),xlab="Tau", ylab="Error",col="blue")
  # 
  # ######Time 2 Error AGAINST TAU 
  
#PLOT for b=1:5
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
  
  