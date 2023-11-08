library(MASS)
library(psych)
library(copula)
library(VineCopula)
library(gamlss)
library(moments)
require(ggpubr)
require(ggplot2)
require(dglm)

#####Try and fit copula to the most extreme case

set.seed(1000)
options(scipen=999)
a=0.25; b=1.75 #.25, .75 the most extreme
#a=0.75; b=1.25 #First MLE calc
mu1=1; mu2=2; n=1000 #100,500,1000,5000,n=10000


##########NADARAJAH & GUPTA SIM
w<-rbeta(n,a,b) #Mean .5
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1) #Mean 6 * .5 = 3
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2) #Mean 12 * .5 = 6

##########is this an RE?
#gamma_c_mu1<-rgamma(n,shape=3,scale=2) #Mean 6 * .5 = 3
#gamma_c_mu2<-gamma_c_mu1 + rgamma(n,shape=4,scale=5) #Mean 12 * .5 = 6


#cor(gamma_c_mu1,gamma_c_mu2,method="kendall")
#skewness(gamma_c_mu1)
#skewness(gamma_c_mu2)

patient<-as.factor(seq(1:n))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")

u<-0
v<-0

u<-pgamma(gamma_c_mu1,shape=a,scale=1/mu1)
v<-pgamma(gamma_c_mu2,shape=a,scale=1/mu2)

#u<-pgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])
#v<-pgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])

#u<-pgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])
#v<-pgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])


#################BEST COPULA FIT
require(GJRM)

#Setting up GJRM equations
eq.mu.1 <- gamma_c_mu1~1
eq.mu.2 <- gamma_c_mu2~1
fl <- list(eq.mu.1, eq.mu.2)

#Running GJRM for each of the copulas tested
#"N", "C0", "GAL0", "J0", "G0", "F", "AMH", "FGM", "T", "PL", "HO"

#"N", "C0", "C90", "C180", "C270", "GAL0", "GAL90", "GAL180", "GAL270", "J0", "J90", "J180", "J270", "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO"

#BiCopSelect(u,v)

model_copula<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
model_copula_n<-gjrm(fl, margins = c("GA" , "GA") , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_j<-gjrm(fl, margins = c("GA" , "GA") , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_g<-gjrm(fl, margins = c("GA" , "GA") , copula = "G0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_f<-gjrm(fl, margins = c("GA" , "GA") , copula = "F",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_amh<-gjrm(fl, margins = c("GA" , "GA") , copula = "AMH",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_fgm<-gjrm(fl, margins = c("GA" , "GA") , copula = "FGM",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_pl<-gjrm(fl, margins = c("GA" , "GA") , copula = "PL",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
#model_copula_h<-gjrm(fl, margins = c("GA" , "GA") , copula = "HO",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")

#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "T",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik ############7575.195
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C180",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik 
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik ############ 7575.195
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J180",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik ############7575.195
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G180",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL=90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL180",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik

#Weird copulas
#"C0C90", "C0C270", "C180C90", "C180C270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270", "G0G90", "G0G270", "G180G90", "G180G270", "J0J90", "J0J270", "J180J90" and "J180J270"
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0C90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0C270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C180C90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "C180C270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL0GAL90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL0GAL270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL180GAL90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "GAL180GAL270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G0G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G0G270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G180G90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "G180G270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J0J90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J0J270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J180J90",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik
#model_copula_c180<-gjrm(fl, margins = c("GA" , "GA") , copula = "J180J270",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B"); model_copula_c180$logLik

#aics=c(
#  model_copula$logLik
#  ,model_copula_n$logLik
#  ,model_copula_j$logLik
#  ,model_copula_g$logLik
#  ,model_copula_f$logLik
#  ,model_copula_amh$logLik
#  ,model_copula_fgm$logLik
#  ,model_copula_pl$logLik
#  ,model_copula_h$logLik)


####UKNOWN FIT u,v

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

fittedClayton=rCopula(n, archmCopula(family="clayton",BiCopSelect(u,v,family=3)$par,dim=2))
fittedTDist=rCopula(n, tCopula(BiCopSelect(u,v,family=2)$par,dim=2,df=BiCopSelect(u,v,family=2)$par2))

# Plot density as points
z<-ggplot(data = as.data.frame(gamma_c_mu1)) +
  geom_histogram(data = as.data.frame(gamma_c_mu1), aes(x=gamma_c_mu1, y=..density..),bins=75) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu1, y=dgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])), color="blue", size = .75) +
  #ylim(0,10) +
  #xlim(0,2) +
  labs(x="time 1 margin",title="time 1 margin") +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
x<-ggplot(data = as.data.frame(gamma_c_mu2)) +
  geom_histogram(data = as.data.frame(gamma_c_mu2), aes(x=gamma_c_mu2, y=..density..),bins=50) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu2, y=dgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])), color="blue", size = .75) +
  #ylim(0,10) +
  labs(x="time 2 margin",title="time 2 margin") +
  #xlim(0,1) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
c<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
  geom_point(size=0.5,color="black") + 
  labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2") +
  #xlim(0,2) +
  #ylim(0,1) +
  geom_smooth(method="loess", level=.99) 
d<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=20,color="black") + 
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2 (unif. transform)",fill="density")
e<-ggplot(data=as.data.frame(fittedClayton),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=20,color="black") +
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")
f<-ggplot(data=as.data.frame(fittedTDist),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d( contour_var="density",bins=20,color="black") + 
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted normal copula")
ggarrange(z,x,c,nrow=1)
##ggsave(file="example_bias_case_margin_plots.jpeg",last_plot(),width=14,height=4,dpi=300)
ggarrange(d,e,f,common.legend = TRUE,nrow=1,legend="right")
##ggsave(file="example_bias_case_contour_plots.jpeg",last_plot(),width=14,height=4,dpi=300)

#plot(u,v,main="Uniform transform of both marginals",xlab="Time 1 Marginal Gamma (Uniform Transform)",ylab="Time 2 Marginal Gamma (Uniform Transform)")
#plot(fittedClayton,main="Simulation of Fitted Clayton Copula",xlab="Fitted Time 1 Marginal Gamma",ylab="Fitted Time 2 Marginal Gamma")
#plot(fittedTDist,main="Simulation of Fitted Normal Copula",xlab="Fitted Time 1 Marginal Gamma",ylab="Fitted Time 2 Marginal Gamma")

#par(mfrow=c(1,3))
#persp(kde2d(u,v,h=.4,n=65),main="Uniform transform of marginals",zlim=c(0,4))
#persp(kde2d(fittedClayton[,1],fittedClayton[,2],h=.4,n=65),main="Simulated fitted copula",zlim=c(0,4))
#persp(kde2d(fittedTDist[,1],fittedTDist[,2],h=.4,n=65),main="Simulated fitted normal",zlim=c(0,4))

require(gamlss)
require(gee)
require(lme4)
require(mgcv)
library(geepack)

model_glm <- glm(random_variable~as.factor(time==1),data=dataset,family=Gamma(link = "log"),maxit=10000)

model_gee<-geese(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link="log"), mean.link = "log", corstr = "exchangeable",control=geese.control(trace=TRUE,maxit=10000))#model_lme4<-glmer(random_variable~as.factor(time==1) + (1 | patient), data=dataset, family=Gamma(link="log")) #lme4
model_gamm<-gamm(random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=Gamma(link="log"),niterPQL=1000) #mgcv
model_re_nosig <- gamlss(random_variable~as.factor(time==1)+random(as.factor(patient)),data=dataset,family=GA(),method=CG(10000))
model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                   , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
                   , method=CG(10000))

summary_glm<-c( summary(model_glm)$coeff[1]
                ,summary(model_glm)$coeff[2] + summary(model_glm)$coeff[1]
                ,summary(model_glm)$coeff[3]
                ,summary(model_glm)$coeff[4]
)
summary_gee<-c( summary(model_gee)$mean[1,1]
                ,summary(model_gee)$mean[2,1] + summary(model_gee)$mean[1,1]
                ,summary(model_gee)$mean[1,2] #Robust SE - look into this
                ,summary(model_gee)$mean[2,2] #Robust SE - look into this
)

invisible(capture.output(
  summary_re_nosig<-c( summary(model_re_nosig)[1]
                         ,summary(model_re_nosig)[2] + summary(model_re_nosig)[1]
                       ,summary(model_re_nosig)[4]
                       ,summary(model_re_nosig)[5]
  )
))
invisible(capture.output(
  summary_re<-c( summary(model_re)[1]
                 ,summary(model_re)[2] + summary(model_re)[1]
                 ,summary(model_re)[5]
                 ,summary(model_re)[6]
  )
))
actuals<-c( log(a*(1/mu1))
            , log(a*(1/mu2))
            , 0#model_copula$tau
            , 0
)

plot()
summary(model_re)
#Time 1
plot(model_re_nosig$mu.fv[1:2000], model_re_nosig$y[1:2000],xlim=c(0,.5),ylim=c(0,.5))
lines(c(0, 5), c(0, 5), type='l',col="red")
plot(model_re$mu.fv[1:2000], model_re$y[1:2000],xlim=c(0,.5),ylim=c(0,.5))
lines(c(0, 5), c(0, 5), type='l',col="red")
#Time 2
plot(model_re_nosig$mu.fv[2001:4000], model_re_nosig$y[2001:4000],xlim=c(0,.5),ylim=c(0,.5))
lines(c(0, 5), c(0, 5), type='l',col="red")
plot(model_re$mu.fv[2001:4000], model_re$y[2001:4000],xlim=c(0,.5),ylim=c(0,.5))
lines(c(0, 5), c(0, 5), type='l',col="red")





plot(model_re$mu.fv[1:2000]/model_re$mu.fv[2001:4000], model_re$y[1:2000]/model_re$y[2001:4000])

plot(model_re$mu.fv[1:2000], model_re$y,xlim=c(0,.5),ylim=c(0,.5))
lines(c(0, 5), c(0, 5), type='l',col="red")



plot(model_glm$fitted.values,model_glm$data$random_variable)
plot(model_glm$fitted.values,model_glm$data$random_variable)


e<-ggplot(data=as.data.frame(cbind(model_re_nosig$mu.fv, model_re_nosig$y)),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=20,color="black") +
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")

f<-ggplot(data=as.data.frame(cbind(model_re$mu.fv, model_re$y)),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=20,color="black") +
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")

ggarrange(e,f)

require(GJRM)
eq.mu.1 <- gamma_c_mu1~1
eq.mu.2 <- gamma_c_mu2~1
fl <- list(eq.mu.1, eq.mu.2)
model_copula<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_n<-gjrm(fl, margins = c("GA" , "GA") , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")

#AIC for copula

model_glm$aic
summary(model_glm)
2*4-2*model_copula_n$logLik
2*6-2*model_copula$logLik
summary(model_copula)
model_re_nosig$aic

2*6-2*logLik(model_re_nosig)
summary(model_re_nosig)
model_re$aic


#time<-c(start_time[2:5]-start_time[1:4],start_time[7:8]-start_time[6:7])
#time
#plot(time)

summary_cop<-c( model_copula$coefficients[1]
                , model_copula$coefficients[2]
                , summary(model_copula)$tableP1[2] #SE for time 0
                , summary(model_copula)$tableP2[2] #SE for time 1
)
summary_cop_n<-c( model_copula_n$coefficients[1]
                  , model_copula_n$coefficients[2] 
                  , summary(model_copula_n)$tableP1[2] #SE for time 0
                  , summary(model_copula_n)$tableP2[2] #SE for time 1
                  
)

rbind(summary_glm, summary_gee,summary_re_nosig,summary_re,summary_cop,summary_cop_n,actuals)
exp(rbind(summary_glm, summary_gee,summary_re_nosig,summary_re,summary_cop,summary_cop_n,actuals))

#Results for bias case (EXP)
#(Intercept) (Intercept)                  
#summary_glm       0.27790949 0.133037154 1.045690 1.065221
#summary_gee       0.27790948 0.133037154 1.046683 1.065204
#summary_re_nosig  0.01844942 0.008929757 1.011789 1.016713
#summary_re        0.01619289 0.014708712 1.000906 1.023602
#summary_cop       0.28406920 0.135318812 1.042085 1.041990
#summary_cop_n     0.41422316 0.197520118 1.054372 1.054211
#actuals           0.25000000 0.125000000 1.000000 1.000000


#Results for RE model
#(Intercept) (Intercept)                        
#summary_glm         1.802780   1.4608502 0.011221056 0.015868969
#summary_gee         1.802780   1.4608502 0.012989056 0.015865001
#summary_re_nosig    1.718899   1.5016149 0.008028506 0.011354022
#summary_re          1.798033   1.3814264 0.008070156 0.008070156
#summary_cop         1.809899   1.4542516 0.012965674 0.009114965
#summary_cop_n       1.802569   1.4611255 0.012861624 0.009052048
#actuals            -1.386294  -0.6931472 0.000000000 0.000000000

# 
# #############
# 
# require(lme4)
# help(lme4)
# 
# summary(model_copula)
# 
# #########RE is the conditional model
# 
# model_re_nosig$mu.coefSmo
# 
# summary(model_lme4)$varcor
# summary(model_lme4)$coefficients
# 
# plot(model_lme4)
# 
# s=3.94481
# l=0.0944163 
# n=100
# 
# base=rnorm(n,0,3.32218^2)
# 
# mean=exp(-1.31)
# var=0.619
# shape=((mean)^2)/var
# scale=var/mean
# t1=base+rgamma(n,shape=shape,scale=scale)
# mean=exp(-1.31-0.7070887)
# var=0.619
# shape=((mean)^2)/var
# scale=var/mean
# t2=base+rgamma(n,shape=shape,scale=scale)
# 
# hist(base,main="Distribution of Estimated Random Effect")
# hist(t1); hist(t2)
# 
# 
# rbind (c(var(gamma_c_mu1),cov(gamma_c_mu1,gamma_c_mu2)),
#        c(cov(gamma_c_mu1,gamma_c_mu2),var(gamma_c_mu2)))

#########################################
# 
# set.seed(100)
# 
# mu1=1; mu2=2; n=100;a=10;s=3
# base=rnorm(n,mean=0,sd=s)
# 
# gamma_c_mu1<-base+rgamma(n,shape=a,scale=1/mu1) #Mean 6 * .5 = 3
# gamma_c_mu2<-base+rgamma(n,shape=a,scale=1/mu2) #Mean 12 * .5 = 6
# 
# cor(gamma_c_mu1,gamma_c_mu2,method="kendall")
# 
# patient<-as.factor(seq(1:n))
# dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
# colnames(dataset)<-c("patient","random_variable","time")
# 
# par(mfrow=c(2,2))
# #u<-pgamma(gamma_c_mu1,shape=a,scale=1/mu1)
# #v<-pgamma(gamma_c_mu2,shape=a,scale=1/mu2)
# u<-pgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])
# v<-pgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])
# his<-kde2d(u,v,h=.4)
# plot(u,v)
# persp(his,axes=T,scale=T,ticktype="detailed",theta=15)
# 
# require(gamlss); require(gee)
# 
# model_glm <- glm(random_variable~as.factor(time==1),data=dataset,family=Gamma(link = "log"),maxit=10000)
# model_re_nosig <- gamlss(random_variable~as.factor(time==1)+random(as.factor(patient)),data=dataset,family=GA(),method=RS(1000))
# model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
#                    , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
#                    , start.from	= model_re_nosig, method=CG(1000))
# 
# model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"),maxiter=10000,silent=TRUE)
# 
# 
# 
# 
# #####OLD
# 
# 
# head(as.data.frame(gamma_c_mu1))
# 
# ggplot(as.data.frame(gamma_c_mu1), aes(x=gamma_c_mu1)) + geom_histogram()
# 
# plot.new()
# par(mfrow=c(1,3))
# histDist(gamma_c_mu1, family=GA(),main="Time 1 Marginal Gamma (mu=0.25, sd=0.5)"
#          ,xlab="x",ylab="Density",nbins=200,xlim=c(0,.5)
#          ,col.main = "black"
#          ,col.lab = "black",
#          ,col.axis = "black"
#          ,line.col="blue"
#          ,col.hist="gray"
#          , border.hist	="black"
#          ,fg.hist="black")
# histDist(gamma_c_mu2, family=GA(),main="TIme 2 Marginal Gamma (mu=0.125, sd=0.25)"
#          ,xlab="x",ylab="Density",nbins=200,xlim=c(0,.25)
#          ,col.main = "black"
#          ,col.lab = "black",
#          ,col.axis = "black"
#          ,line.col="blue"
#          ,col.hist="gray"
#          , border.hist	="black"
#          ,fg.hist="black")
# plot(gamma_c_mu1,gamma_c_mu2,xlab="Time 1 Marginal Gamma",ylab="Time 2 Marginal Gamma",main="Time 2 v Time 1 Marginal Gamma")
# 
