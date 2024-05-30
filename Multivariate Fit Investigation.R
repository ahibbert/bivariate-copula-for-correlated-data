source("common_functions.R")

# a. Simulation parameters
set.seed(1000);options(scipen=999);

### RAND MULTIVARIATE ###

library(haven)
randFULL <- read_sas("Data/randhrs1992_2020v1.sas7bdat")
docnames<-(grep("R[0-9]+DOCTIM", colnames(randFULL),value=TRUE))
docnames_sorted<-docnames[order(as.numeric(regmatches(docnames,regexpr("[0-9]+",docnames))))]

rand_allyearsdocvisits<-randFULL[,c("HHIDPN","RABYEAR",docnames_sorted)]

#1775 individuals with no missing data
rand_mvt<-as.data.frame(na.omit(rand_allyearsdocvisits))
#save(rand_mvt,file="rand_mvt.rds")
load("rand_mvt.rds")

#plot(rand_mvt[,2:17])
#cor(rand_mvt[,2:17])
#cor(rand_mvt[,2:17],method="kendall")

library(gamlss)
fit<-fitDist(rand_mvt[,3],type = "counts")
fit$fits

patient<-as.factor(seq(1:nrow(rand_mvt)))
dataset<-matrix(data=NA,ncol=3,nrow=0)

for (t in 1:length(docnames)) {
  dataset_temp<-cbind(patient,t,rand_mvt[,docnames[t]])
  dataset <- rbind(dataset,dataset_temp)
}

colnames(dataset)<-c("patient","time","random_variable")
dataset<-dataset[order(dataset[,"patient"]),]

#####################################Copula fits

plot.new()
par(mfrow=c(3,5))
for (i in 1:length(docnames)) {
  histDist(dataset[dataset[,"time"]==i,"random_variable"],family=ZISICHEL,main=paste("Time",i),xlim=c(0,100),ylim=c(0,.2))
}
plot(1:15,colMeans(rand_mvt)[3:17],main="Mean no. doctor visits in year",xlab="Time",ylab="Mean doctor visits")

par(mfrow=c(1,1))
boxplot(dataset[,"random_variable"]~dataset[,"time"],ylim=c(0,20),ylab="doctor visits",xlab="time")
boxplot(rand_mvt$R15DOCTIM~rand_mvt$RABYEAR,ylim=c(0,40))
  
gamlss_fits <- list()
family_spec=ZISICHEL
fit_no<-matrix(NA,ncol=0,nrow=nrow(dataset[dataset[,"time"]==1,]))
fit_unif<-matrix(NA,ncol=0,nrow=nrow(dataset[dataset[,"time"]==1,]))
for (i in 1:length(docnames)) {
  gamlss_fits[[i]]<-gamlss(dataset[dataset[,"time"]==i,"random_variable"]~1,family=family_spec,method=RS(100))
  fit_no<-cbind(fit_no,(gamlss_fits[[i]]$residuals))
  fit_unif<-cbind(fit_unif,pNO(gamlss_fits[[i]]$residuals))
}
colnames(fit_no)<-colnames(fit_unif)<-docnames_sorted

library(ggplot2)
library(ggpubr)

#pairs(fit_no)
#pairs(fit_unif)

copFits<-list()
plots<-list()
library(VineCopula)
for (i in 1:(length(docnames)-1)) {
  copFits[[i]]<-BiCopSelect(fit_unif[,i],fit_unif[,i+1])  
  data<-as.data.frame(fit_unif[,c(i,i+1)])
  colnames(data)<-c("V1","V2")
  plots[[i]]<-ggplot(data) +
    geom_density_2d(aes(x=V1,y=V2))
}


ggarrange(plotlist=plots)

library(VineCopula)
vinefit<-RVineStructureSelect(fit_unif[,c(11,12,13,14,15)])
summary(vinefit)
contour(vinefit)

vinefit$pair.logLik

#library(network);plot(vinefit)


#########Random effect fit

gamlss_glm_fit_nbi<-gamlss(random_variable~time                                                 ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=NBI,method=RS(100))
gamlss_glm_fit_zis<-gamlss(random_variable~time                                                 ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_glm_fit_zis_timecat<-gamlss(random_variable~as.factor(time)                              ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_re_fit_zis<-gamlss(random_variable~time+random(as.factor(patient))                       ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_re_fit_zis_timecat<-gamlss(random_variable~as.factor(time)+random(as.factor(patient))    ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))

term.plot(gamlss_re_fit_zis_timecat,ylim="free")

results<-as.data.frame(rbind(
c(logLik(gamlss_glm_fit_nbi), "GLM NBI + time", gamlss_glm_fit_nbi$df.fit)
,c(logLik(gamlss_glm_fit_zis), "GLM ZIS  + time", gamlss_glm_fit_zis$df.fit)
,c(logLik(gamlss_glm_fit_zis_timecat), "GLM ZIS + factor(time)",gamlss_glm_fit_zis_timecat$df.fit)
,c(logLik(gamlss_re_fit_zis), "GLMM ZIS + time",gamlss_re_fit_zis$df.fit)
,c(logLik(gamlss_re_fit_zis_timecat), "GLMM ZIS + factor(time)",gamlss_re_fit_zis_timecat$df.fit)
,c(logLik(gamlss_fits[[11]])+logLik(gamlss_fits[[12]])+logLik(gamlss_fits[[13]])+logLik(gamlss_fits[[14]])+logLik(gamlss_fits[[15]])+vinefit$logLik,"GJRM ~ factor(time)",4*5+10*2)
))

colnames(results)<-c("LogLik","Model","EDF")
results<-results[c("Model","LogLik","EDF")]
results$LogLik<-round(as.numeric(results$LogLik))
results$EDF<-round(as.numeric(results$EDF))

results$AIC2 <- results$LogLik*-2+results$EDF*2
results$AIC4 <- results$LogLik*-2+results$EDF*4
results$BIC <- round(results$LogLik*-2+results$EDF*log(nrow(rand_mvt)))

#sum of marginal fits




