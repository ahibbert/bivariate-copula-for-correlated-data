#Lipids Data
  #require(sas7bdat)
  #lipid <- read.sas7bdat("LipidsData.sas7bdat")
  #lipids_merged<-(merge(lipid[lipid$MONTH==0,],lipid[lipid$MONTH==24,],by="PATIENT"))
  lipids_merged<-readRDS("lipids_merged.rds")
  gamma_c_mu1<-lipids_merged$TRG.x
  gamma_c_mu2<-lipids_merged$TRG.y
  
#Stock prices over 10 years

  #ASX2018<-read.table("20180102.txt", header=FALSE, sep=",")
  #ASX1998<-read.table("19980102.txt", header=FALSE, sep=",")
  #ASX98_18<-merge(ASX1998,ASX2018,by="V1")
  ASX98_18<-readRDS("ASX98_18.rds")
  gamma_c_mu1<-ASX98_18$V6.x
  gamma_c_mu2<-ASX98_18$V6.y

#Avocado prices
  #avo<-read.table("avocado prices.csv", header=T, sep=",")
  avo<-readRDS("avo.rds")
  gamma_c_mu1<-avo[avo$Date=="4/01/2015","AveragePrice"]
  gamma_c_mu2<-avo[avo$Date=="25/03/2018","AveragePrice"]

  
############# NEW PLOTTING
  
  library(moments)
  skewness(gamma_c_mu1)
  skewness(gamma_c_mu2)
  
  require(MASS)  
  require(copula)
  require(gamlss)
  require(dglm)
  require(ggplot2)
  require(ggpubr)
  u<-0
  v<-0
  
  fit <- dglm(gamma_c_mu1~1, family=Gamma(link="log"))
  mu=exp(fit$coefficients)
  shape=exp(-1*fit$dispersion.fit$coefficients)
  scale <- mu/shape
  
  u<-pgamma(gamma_c_mu1,shape=shape,scale=scale)
  
  fit <- dglm(gamma_c_mu2~1, family=Gamma(link="log"))
  mu=exp(fit$coefficients)
  shape=exp(-1*fit$dispersion.fit$coefficients)
  scale <- mu/shape
  
  v<-pgamma(gamma_c_mu2,shape=shape,scale=scale)

  #ASX# bins=30; binsc=15; limx=50; limy=0.15
  #lipids 
  bins=20; binsc=15; limx=3; limy=2
  # Plot density as points
  a<-ggplot(data = as.data.frame(gamma_c_mu1)) +
    geom_histogram(data = as.data.frame(gamma_c_mu1), aes(x=gamma_c_mu1, y=..density..),bins=bins) +
    geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu1, y=dgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])), color="blue", size = .75) +
    ylim(0,limy) +
    xlim(0,limx) +
    labs(x="time 1 margin",title="time 1 margin") +
    theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
  b<-ggplot(data = as.data.frame(gamma_c_mu2)) +
    geom_histogram(data = as.data.frame(gamma_c_mu2), aes(x=gamma_c_mu2, y=..density..),bins=bins) +
    geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu2, y=dgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])), color="blue", size = .75) +
    ylim(0,limy) +
    labs(x="time 2 margin",title="time 2 margin") +
    xlim(0,limx) +
    theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
  c<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
    geom_point(size=0.5,color="black") + 
    labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2") +
    geom_smooth(method="loess", level=.95) +
    xlim(0,limx) +
    ylim(0,limx) 
  d<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
    #geom_point(size=0.25,color="black") + 
    geom_density_2d(contour_var="density",bins=binsc,color="black") + 
    scale_fill_brewer() +
    labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2 (unif. transform)",fill="density")
  #e<-ggplot(data=as.data.frame(fittedClayton),aes(x=V1,y=V2)) + 
  #  #geom_point(size=0.25,color="black") + 
  #  geom_density_2d(contour_var="density",bins=20,color="black") +
  #  scale_fill_brewer() +
  #  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")
  #f<-ggplot(data=as.data.frame(fittedTDist),aes(x=V1,y=V2)) + 
  #  #geom_point(size=0.25,color="black") + 
  #  geom_density_2d( contour_var="density",bins=20,color="black") + 
  #  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted normal copula")
  ggarrange(a,b,d,nrow=1)
  #ggsave(file="applications_asx.jpeg",last_plot(),width=14,height=4,dpi=300)
  #ggarrange(d,e,f,common.legend = TRUE,nrow=1,legend="right")
  ##ggsave(file="example_bias_case_contour_plots.jpeg",last_plot(),width=14,height=4,dpi=300)
  
  
  
########################LOOK at shape
# 
#   his<-kde2d(u,v,h=.4)
#   par(mfrow=c(2,3))
#   
#   paste("Cor:",cor(gamma_c_mu1,gamma_c_mu2)," | Tau:",cor(gamma_c_mu1,gamma_c_mu2,method = "kendall"))
#         
#   scatter.smooth(gamma_c_mu1,gamma_c_mu2,xlab="Time 1",ylab="Time 2",main="Time 1 v Time 2")
#   scatter.smooth(u,v, main="Unif. Transform of Time 1 v Time 2 (Scatter)", xlab="Uniform Time 1",ylab = "Uniform Time 2")
#   persp(his,axes=T,scale=T,ticktype="detailed",theta=15, main="Unif. Transform of Time 1 v Time 2 (3D Density)",zlim=c(0,3),xlab="Uniform Time 1",ylab="Uniform Time 2",zlab="Density")
#   
#   fittedCopula<-rCopula(10000,copula=claytonCopula(coef(fitCopula(claytonCopula(),data=cbind(u,v))), dim = 2))
#   fittedNormal<-rCopula(10000,copula=normalCopula(coef(fitCopula(normalCopula(dim=2),data=cbind(u,v))), dim = 2))
#   his<-kde2d(fittedCopula[,1],fittedCopula[,2],h=.25)
#   hisn<-kde2d(fittedNormal[,1],fittedNormal[,2],h=.25)
#   require(VineCopula)
#   fits<-BiCopSelect(u,v,family=c(1,2,3,4,5,6),indeptest = TRUE)
#   persp(his,axes=T,scale=T,ticktype="detailed",theta=15, main="Unif. Transform of Time 1 v Time 2 (3D Density)",zlim=c(0,3),xlab="Uniform Time 1",ylab="Uniform Time 2",zlab="Density")
#   persp(his,axes=T,scale=T,ticktype="detailed",theta=15,main=paste("Fitted Clayton"),zlim=c(0,3))
#   persp(hisn,axes=T,scale=T,ticktype="detailed",theta=15,main=paste("Fitted Normal"),zlim=c(0,3))
#   fits$familyname
  
################ACTUAL MODELLING

set.seed(1000)
options(scipen=999)
#Setting up data in correct structure for random effect
n=length(gamma_c_mu1)
patient<-as.factor(seq(1:n))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")

#Loading required pacakges
require(gamlss); require(gee); require(geepack);

#Running each of the models for the application dataset
model_glm <- glm(random_variable~as.factor(time==1),data=dataset,family=Gamma(link = "log"),maxit=10000)
model_gee<-geese(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link="log"), mean.link = "log",corstr = "exchangeable",control=geese.control(trace=TRUE,maxit=10000))
model_re_nosig <- gamlss(random_variable~as.factor(time==1)+random(as.factor(patient))
                         ,data=dataset,family=GA(),method=RS())
model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                   , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
                   , method=CG(10000))  

#Estimate time 1,2
summary_glm<-c( summary(model_glm)$coeff[1]
                ,summary(model_glm)$coeff[2]
                ,summary(model_glm)$coeff[3]
                ,summary(model_glm)$coeff[4]
)


summary_gee<-c( summary(model_gee)$mean[1,1]
                ,summary(model_gee)$mean[2,1]
                ,summary(model_gee)$mean[1,2] #Robust SE - look into this
                ,summary(model_gee)$mean[2,2] #Robust SE - look into this
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
#Loading GJRM after gamlss to avoid overlap
require(GJRM)

#Setting up GJRM equations
eq.mu.1 <- gamma_c_mu1~1
eq.mu.2 <- gamma_c_mu2~1
fl <- list(eq.mu.1, eq.mu.2)

#Running GJRM for each of the copulas tested

model_copula<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
model_copula_n<-gjrm(fl, margins = c("GA" , "GA") , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_j<-gjrm(fl, margins = c("GA" , "GA") , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_g<-gjrm(fl, margins = c("GA" , "GA") , copula = "G0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_f<-gjrm(fl, margins = c("GA" , "GA") , copula = "F",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_amh<-gjrm(fl, margins = c("GA" , "GA") , copula = "AMH",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_fgm<-gjrm(fl, margins = c("GA" , "GA") , copula = "FGM",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_pl<-gjrm(fl, margins = c("GA" , "GA") , copula = "PL",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_h<-gjrm(fl, margins = c("GA" , "GA") , copula = "HO",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")

  aics=c(
  model_copula$logLik
  ,model_copula_n$logLik
  ,model_copula_j$logLik
  ,model_copula_g$logLik
  ,model_copula_f$logLik
  ,model_copula_amh$logLik
  ,model_copula_fgm$logLik
  ,model_copula_pl$logLik
  ,model_copula_h$logLik)
  
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
  summary_cop_j<-c( model_copula_j$coefficients[1]
                    , model_copula_j$coefficients[2] - model_copula_j$coefficients[1]
                    , summary(model_copula_j)$tableP1[2] #SE for time 0
                    , summary(model_copula_j)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_g<-c( model_copula_g$coefficients[1]
                    , model_copula_g$coefficients[2] - model_copula_g$coefficients[1]
                    , summary(model_copula_g)$tableP1[2] #SE for time 0
                    , summary(model_copula_g)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_f<-c( model_copula_f$coefficients[1]
                    , model_copula_f$coefficients[2] - model_copula_f$coefficients[1]
                    , summary(model_copula_f)$tableP1[2] #SE for time 0
                    , summary(model_copula_f)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_amh<-c( model_copula_amh$coefficients[1]
                    , model_copula_amh$coefficients[2] - model_copula_amh$coefficients[1]
                    , summary(model_copula_amh)$tableP1[2] #SE for time 0
                    , summary(model_copula_amh)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_fgm<-c( model_copula_fgm$coefficients[1]
                    , model_copula_fgm$coefficients[2] - model_copula_fgm$coefficients[1]
                    , summary(model_copula_fgm)$tableP1[2] #SE for time 0
                    , summary(model_copula_fgm)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_pl<-c( model_copula_pl$coefficients[1]
                    , model_copula_pl$coefficients[2] - model_copula_pl$coefficients[1]
                    , summary(model_copula_pl)$tableP1[2] #SE for time 0
                    , summary(model_copula_pl)$tableP2[2] #SE for time 1
                    
  )
  summary_cop_h<-c( model_copula_h$coefficients[1]
                    , model_copula_h$coefficients[2] - model_copula_h$coefficients[1]
                    , summary(model_copula_h)$tableP1[2] #SE for time 0
                    , summary(model_copula_h)$tableP2[2] #SE for time 1
                    
  )
  
  application_model_summary <- rbind(summary_glm, summary_gee,summary_re_nosig,summary_re,summary_cop,summary_cop_n
        ,summary_cop_j
        ,summary_cop_g
        ,summary_cop_f
        ,summary_cop_amh
        ,summary_cop_fgm
        ,summary_cop_pl
        ,summary_cop_h
        )
  
  
  application_model_summary2<-cbind(application_model_summary,exp(application_model_summary[,1]),exp(application_model_summary[,2]+application_model_summary[,1]))
  colnames(application_model_summary2) <- c("bm1","bm2","se1","se2","ebm1","ebm1plus2")
  application_model_summary2
  
  #help("gjrm")
  #summary(model_copula$logLik)
  
  #model_glm$aic
  #model_re_nosig$aic
  #5*2-2*model_copula$logLik
  