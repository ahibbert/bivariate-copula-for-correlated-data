source("common_functions.R")

library(gamlss)
library(e1071)

require(sas7bdat)
library(haven)

#hars <- read.sas7bdat("Data/h20n_r.sas7bdat")
#rand <- read.sas7bdat("Data/randhrs1992_2020v1.sas7bdat")

rand <- read_sas("Data/randhrs1992_2020v1.sas7bdat",col_select=c("HHIDPN","R14IWENDY","R15IWENDY"
                                                                             ,"R14HSPTIM" #Number of hospital stays last 2 years
                                                                             ,"R15HSPTIM" #Number of hospital stays last 2 years
                                                                             ,"R14HSPNIT" #Nights in hospital last 2 years
                                                                             ,"R15HSPNIT" #Nights in hospital last 2 years
                                                                             ,"R14DOCTIM" #Doctor visits last 2 years
                                                                             ,"R15DOCTIM" #Doctor visits last 2 years
                                                                             ,"R14OOPMD" #Out of pocket medical exp. last 2 years
                                                                             ,"R15OOPMD" #Out of pocket medical exp. last 2 years
))

rand <- read_sas("Data/randhrs1992_2020v1.sas7bdat",col_select=c("HHIDPN","R14IWENDY","R15IWENDY"
                                                                 ,"R14DOCTIM" #Doctor visits last 2 years
                                                                 ,"R15DOCTIM" #Doctor visits last 2 years
))

rand_hosp_visits<-as.data.frame(rand[!(is.na(rand$R14HSPTIM))&!(is.na(rand$R15HSPTIM))&rand$R14HSPTIM>0&rand$R15HSPTIM>0,c("HHIDPN","R14IWENDY","R15IWENDY","R14HSPTIM","R15HSPTIM")])
rand_hosp_nights<-as.data.frame(rand[!(is.na(rand$R14HSPNIT))&!(is.na(rand$R15HSPNIT))&rand$R14HSPNIT>0&rand$R15HSPNIT>0,c("HHIDPN","R14IWENDY","R15IWENDY","R14HSPNIT","R15HSPNIT")])
rand_doc_visits <-as.data.frame(rand[!(is.na(rand$R14DOCTIM))&!(is.na(rand$R15DOCTIM))&rand$R14DOCTIM>0&rand$R15DOCTIM>0, c("HHIDPN","R14IWENDY","R15IWENDY","R14DOCTIM","R15DOCTIM")])
rand_oop_exp    <-as.data.frame(rand[!(is.na(rand$R14OOPMD))&!(is.na(rand$R15OOPMD))&rand$R14OOPMD>0&rand$R15OOPMD>0,        c("HHIDPN","R14IWENDY","R15IWENDY","R14OOPMD","R15OOPMD")])

###Updated to not exclude zero***
rand_doc_visits <-as.data.frame(rand[!(is.na(rand$R14DOCTIM))&!(is.na(rand$R15DOCTIM))&rand$R14DOCTIM>=0&rand$R15DOCTIM>=0, c("HHIDPN","R14IWENDY","R15IWENDY","R14DOCTIM","R15DOCTIM")])


#2018 & 2020/21

listdata <- list()
listdata[[1]]<-rand_hosp_visits
listdata[[2]]<-rand_hosp_nights
listdata[[3]]<-rand_doc_visits
listdata[[4]]<-rand_oop_exp

par(mfrow=c(4,3))
for (i in 1:4) {
  i=3
  dataset<-listdata[[i]]
  hist(dataset[,4])
  hist(dataset[,5])
  plot(dataset[,4],dataset[,5])
  print(c(cor(dataset[,4],dataset[,5]),cor(dataset[,4],dataset[,5],method="kendall"),skewness(dataset[,4],dataset[,5])))
  
  rand <- read_sas("Data/randhrs1992_2020v1.sas7bdat",col_select=c("HHIDPN","R14IWENDY","R15IWENDY"
                                                                   ,"R14DOCTIM" #Doctor visits last 2 years
                                                                   ,"R15DOCTIM" #Doctor visits last 2 years
  ))
  gamma_c_mu1<-rand[,4]
  gamma_c_mu2<-rand[,5]
  dist="PO"
  
  patient<-as.factor(seq(1:length(gamma_c_mu1)))
  dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
  colnames(dataset)<-c("patient","random_variable","time")
  dataset<-dataset[order(dataset$patient),]
  if(i==4){dist="GA"}else(dist="PO")
  
  fits1<-fitDist(dataset[dataset$time==0,"random_variable"],type="counts")
  fits2<-fitDist(dataset[dataset$time==1,"random_variable"],type="counts")
  
  
  results<-fitBivModels(data=dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=FALSE)
  if(dist=="NO"){clean_results<-results}else{clean_results<-cbind(results,round(exp(results[,c(1,2)]),4));
  clean_results<-cbind(clean_results[,c(9,10)],clean_results[,1:8]);colnames(clean_results)<-c("mu_1","mu_2",colnames(clean_results[,c(3:10)]))}
  print(clean_results)
  
}


library(gamlss)
library(lme4)
library(gee)
library(gamlss.mx)

model_glm <- gamlss(formula=random_variable~-1+as.factor(time==1), data=dataset, family=GA()) 
model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")

model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
model_re_np <- gamlssNP(formula=random_variable~-1+as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                        , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)



model_re_nosig <- gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=ZAPIG())


require(GJRM)

margin_dist="GA"
#Setting up GJRM equations
eq.mu.1 <- formula(gamma_c_mu1~1)
eq.mu.2 <- formula(gamma_c_mu2~1)
fl <- list(eq.mu.1, eq.mu.2)
model_copula<-    gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")





rand[!is.na(rand$R14HSPTIM),]



##Wave is year so we want 15 and 14

rand[,c("HHIDPN","R14IWENDY","R15IWENDY"
        ,"R14HSPTIM" #Number of hospital stays last 2 years
        ,"R15HSPTIM" #Number of hospital stays last 2 years
        ,"R14HSPNIT" #Nights in hospital last 2 years
        ,"R15HSPNIT" #Nights in hospital last 2 years
        ,"R14DOCTOR" #Doctor visits last 2 years
        ,"R15DOCTOR" #Doctor visits last 2 years
        ,"R14OOPMD" #Out of pocket medical exp. last 2 years
        ,"R15OOPMD" #Out of pocket medical exp. last 2 years
        )]
rand$R1IWENDY


#Second number is wave
R15HSPTIM #Number of hospital stays 12m
R15HSPNIT #Nights in hospital 12m
R15DOCTIM #Doctor visits last 12m
R15OOPMD #Out of pocket medical expenditure


save(data, file = "data.Rdata")

summary(hars)

hist(hars$RN106[!is.na(hars$RN106)&hars$RN106<9999998&hars$RN106>0],breaks=100)

head(hars)




#################asthma and BMI##############
bmi_data_read<-read.csv(file="Data/BMI_IOS_SCD_Asthma.csv")
#bmi_data<-bmi_data_read[(bmi_data_read$Observation_number==1|bmi_data_read$Observation_number==2)&!is.na(bmi_data_read$Fres_PP),]

a<-bmi_data_read[bmi_data_read$Observation_number==1,]
b<- bmi_data_read[bmi_data_read$Observation_number==2,]
m<- merge(a,b,by="Subject.ID",all.x=TRUE)
m_nona<-m[!is.na(m$Observation_number.y),]#&m$X5Hz_PP.x>0&m$X5Hz_PP.y>0


m_nona<-m[!is.na(m$Observation_number.y)&m$X5Hz_PP.x>0&m$X5Hz_PP.y>0,]#&m$X5Hz_PP.x>0&m$X5Hz_PP.y>0
gamma_c_mu1<-m_nona$X5Hz_PP.x
gamma_c_mu2<-m_nona$X5Hz_PP.y


m_nona<-m[!is.na(m$Observation_number.y)&!is.na(m$Fres_PP.x)&!is.na(m$Fres_PP.y),]#&m$X5Hz_PP.x>0&m$X5Hz_PP.y>0
gamma_c_mu1<-m_nona$Fres_PP.x
gamma_c_mu2<-m_nona$Fres_PP.y


############Dementia

data_read<-read.csv(file="Data/dementia.csv")

a<-data_read[data_read$Visit=="x1",]
b<- data_read[data_read$Visit=="x2",]
m<- na.omit(merge(a,b,by="Subject_ID",all.x=TRUE))

gamma_c_mu1<-m$Group.x
gamma_c_mu2<-m$MR_Delay.y
######

library(sas7bdat)

child_health<-read.sas7bdat("Data/hlth_12.sas7bdat")
head(child_health)

child_health_cost_only<-child_health[,c("IDind","wave","M29","M30","M35","M36","M38","M39","M50","M44")]

summary(child_health_cost_only)##Start with 2015 and 2014
unique(child_health_cost_only[,"wave"]) #2011 2015 1989 1991 1993 2000 2004 2006 2009 1997

par(mfrow=c(3,3))
for (i in c("M29","M30","M35","M36","M38","M39","M50","M44")) {
  outcome=child_health_cost_only[child_health_cost_only[,i]>0&!is.na(child_health_cost_only[,i]),i]
  hist(outcome,breaks=100,main=i)  
}

plot.new() # 1989 1991 1993 1997 2000 2004 2006 2009 2011 2015
library(e1071)
source("common_functions.R")
par(mfrow=c(3,3))
for (i in c("M29","M30","M38","M39","M50")) {
  #i="M30"
  years=c(1993,1991)
  cost_2015_2014<-merge(child_health_cost_only[child_health_cost_only[,"wave"]==years[1],c("IDind",i)],child_health_cost_only[child_health_cost_only[,"wave"]==years[2],c("IDind",i)],by="IDind")
  cost_2015_2014_clean<-cost_2015_2014[cost_2015_2014[,2]>0&cost_2015_2014[,3]>0&is.finite(cost_2015_2014[,2])&is.finite(cost_2015_2014[,3])&!(cost_2015_2014[,3] %in% c(999,9999,99999,999999))&!(cost_2015_2014[,2] %in% c(999,9999,99999,999999)),]
  print(c(i,nrow(cost_2015_2014_clean),cor(cost_2015_2014_clean[,2],cost_2015_2014_clean[,3],method="kendall"),skewness(cost_2015_2014_clean[,2])))
  if(nrow(cost_2015_2014_clean)>10) {
    plot(cost_2015_2014_clean[,2],cost_2015_2014_clean[,3])
    
    gamma_c_mu1<-cost_2015_2014_clean[,3]
    gamma_c_mu2<-cost_2015_2014_clean[,2]
    
    patient<-as.factor(seq(1:length(gamma_c_mu1)))
    dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
    colnames(dataset)<-c("patient","random_variable","time")
    dataset<-dataset[order(dataset$patient),]
    dist="GA"
  
    results<-fitBivModels(data=dataset,dist,include="non-GJRM",a,b,c,mu1,mu2,calc_actuals=FALSE)
    if(dist=="NO"){clean_results<-results}else{clean_results<-cbind(results,round(exp(results[,c(1,2)]),4));
    clean_results<-cbind(clean_results[,c(9,10)],clean_results[,1:8]);colnames(clean_results)<-c("mu_1","mu_2",colnames(clean_results[,c(3:10)]))}
    print(clean_results)
  }
}

###Only selecting where value is >0 and not na for both time points


par(mfrow=c(3,3))
length(child_health_cost_only[child_health_cost_only[,3]>0&!is.na(child_health_cost_only[,3]),3])
length(child_health_cost_only[child_health_cost_only[,4]>0&!is.na(child_health_cost_only[,4]),4])
length(child_health_cost_only[child_health_cost_only[,5]>0&!is.na(child_health_cost_only[,5]),5])
length(child_health_cost_only[child_health_cost_only[,6]>0&!is.na(child_health_cost_only[,6]),6])
length(child_health_cost_only[child_health_cost_only[,7]>0&!is.na(child_health_cost_only[,7]),7])
length(child_health_cost_only[child_health_cost_only[,8]>0&!is.na(child_health_cost_only[,8]),8])
length(child_health_cost_only[child_health_cost_only[,9]>0&!is.na(child_health_cost_only[,9]),9])



library(e1071)
skewness(child_health_cost_only[child_health_cost_only[,3]>0&!is.na(child_health_cost_only[,3]),3])
skewness(child_health_cost_only[child_health_cost_only[,4]>0&!is.na(child_health_cost_only[,4]),4])
skewness(child_health_cost_only[child_health_cost_only[,5]>0&!is.na(child_health_cost_only[,5]),5])
skewness(child_health_cost_only[child_health_cost_only[,6]>0&!is.na(child_health_cost_only[,6]),6])
skewness(child_health_cost_only[child_health_cost_only[,7]>0&!is.na(child_health_cost_only[,7]),7])
skewness(child_health_cost_only[child_health_cost_only[,8]>0&!is.na(child_health_cost_only[,8]),8])
skewness(child_health_cost_only[child_health_cost_only[,9]>0&!is.na(child_health_cost_only[,9]),9])


############Diabetes hospitalisations

data_read<-read.csv(file="Data/diabetic_data.csv")

#data_read<-data_read[data_read$number_inpatient>0,]
admissions<-data.frame(table(data_read$patient_nbr))
two_plus_admissions<-data.frame(admissions[admissions$Freq>1, "Var1"])
colnames(two_plus_admissions)<-"patient_nbr"
m<-data.frame(merge(data_read,two_plus_admissions,by="patient_nbr"))

m$admission_nbr <- with(m, ave(seq_along(patient_nbr), patient_nbr, FUN = seq_along))

time1<-m[m$admission_nbr==1&(m$age=="[0-10)"|m$age=="[10-20)"),] 
time2<-m[m$admission_nbr==2,]

m_merge<-merge(time1,time2,by="patient_nbr")

gamma_c_mu1<-m_merge$time_in_hospital.x
gamma_c_mu2<-m_merge$time_in_hospital.y

length(gamma_c_mu1)

skewness(gamma_c_mu1)
skewness(gamma_c_mu2)
cor(gamma_c_mu1,gamma_c_mu2,method = "kendall")
dist="PO"

patient<-as.factor(seq(1:length(gamma_c_mu1)))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")
dataset<-dataset[order(dataset$patient),]

plotDist(dataset,dist)

results<-fitBivModels(data=dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=FALSE)
if(dist=="NO"){clean_results<-results}else{clean_results<-cbind(results,round(exp(results[,c(1,2)]),4));
clean_results<-cbind(clean_results[,c(9,10)],clean_results[,1:8]);colnames(clean_results)<-c("mu_1","mu_2",colnames(clean_results[,c(3:10)]))}
clean_results

library(gamlss)
library(lme4)
library(gee)
library(gamlss.mx)

model_glm <- gamlss(formula=random_variable~as.factor(time==1), data=dataset, family=GA()) 
model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")

model_gee<-gee(random_variable~-1+as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")

model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
model_lme4 <- glmer(formula=random_variable~-1+as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))
model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                        , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)

m_nona[m$Observation_number.x]





skewness(bmi_data$R5Hz_PP[bmi_data$Observation_number==1])
skewness(bmi_data$R20Hz_PP[bmi_data$Observation_number==1])
skewness(bmi_data$X5Hz_PP[bmi_data$Observation_number==1])
skewness(bmi_data$Fres_PP[!is.na(bmi_data$Fres_PP)])


par(mfrow=c(2,2))
hist(bmi_data$R5Hz_PP[bmi_data$Observation_number==1])
hist(bmi_data$R20Hz_PP[bmi_data$Observation_number==1])
hist(bmi_data$X5Hz_PP[bmi_data$Observation_number==1],breaks=100)
hist(bmi_data$X5Hz_PP[bmi_data$Observation_number==2],breaks=100)
hist(bmi_data$Fres_PP[bmi_data$Observation_number==1])


plot(bmi_data$X5Hz_PP[bmi_data$Observation_number==1],bmi_data$X5Hz_PP[bmi_data$Observation_number==2])


par(mfrow=c(2,2))
fitDist(bmi_data$R5Hz_PP)$fits
fitDist(bmi_data$R20Hz_PP)$fits
fitDist(bmi_data$X5Hz_PP)$fits
fitDist(bmi_data$Fres_PP)$fits

histDist(bmi_data$R5Hz_PP,family="GA")
histDist(bmi_data$R20Hz_PP,family="GA")
histDist(bmi_data$X5Hz_PP,family="GA")
histDist(bmi_data$Fres_PP,family="GA")


glm(formula=Fres_PP~-1+as.factor(Observation_number==1), data=bmi_data, family=Gamma(link="log"))
glmer(formula=Fres_PP~-1+as.factor(Observation_number==1) + (1|Subject.ID), data=bmi_data, family=Gamma(link="log"))


glm<-gamlss(list(Fres_PP~(Observation_number)),data=bmi_data,family=GA())
summary(glm)
glmm<-gamlss(list(Fres_PP~(Observation_number)+random(as.factor(Subject.ID))),data=bmi_data,family=GA())
summary(glmm)




par(mfrow=c(3,3))
hist(bmi_data$Fres_PP[bmi_data$Observation_number==1])
hist(bmi_data$Fres_PP[bmi_data$Observation_number==2])
hist(bmi_data$Fres_PP[bmi_data$Observation_number==3])
hist(bmi_data$Fres_PP[bmi_data$Observation_number==4])
hist(bmi_data$Fres_PP[bmi_data$Observation_number==5])
hist(bmi_data$Fres_PP[bmi_data$Observation_number==6])


plot(bmi_data$Fres_PP[bmi_data$Observation_number==1],bmi_data$Fres_PP[bmi_data$Observation_number==2])

help(gamlss)

#gamlss(formula=random_variable~-1+as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA())
