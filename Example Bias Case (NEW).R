############# Overview of file
#############
#############

########### 0. DATA SETUP ##############

# a. Simulation parameters
set.seed(1000)
options(scipen=999)
a=1; b=1; mu1=10; mu2=12; n=1000 #100,500,1000,5000,n=10000

# b. Simualating Nadarajah and Gupta bivariate Gamma
w<-rbeta(n,a,b) #Mean .5
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=mu1) #Mean 6 * .5 = 3
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=mu2) #Mean 12 * .5 = 6

# c.Setting up as longitiduinal structured data
patient<-as.factor(seq(1:n))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")
dataset<-dataset[order(dataset$patient),]

################# 1. Investigating the dependence structure and best copula fit #############
require(GJRM)
library(MASS)
library(psych)
library(copula)
library(VineCopula)
library(gamlss)
library(moments)
require(ggpubr)
require(ggplot2)
require(dglm)

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


################# 2. GLM, GEE, GLMM fits to the data #########

### Running all models
require(gamlss)
require(gee)
require(lme4)
require(mgcv)
library(geepack)

model_glm <- glm(random_variable~as.factor(time==1),data=dataset,family=Gamma(link = "log"),maxit=10000)
model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link="log"), corstr = "exchangeable")#model_lme4<-glmer(random_variable~as.factor(time==1) + (1 | patient), data=dataset, family=Gamma(link="log")) #lme4
model_gamm<-gamm(random_variable~as.factor(time==1), random=list(patient=~1), data=dataset, family=Gamma(link="log"),niterPQL=1000) #mgcv
model_re_nosig <- gamlss(random_variable~as.factor(time==1)+random(as.factor(patient)),data=dataset,family=GA(),method=CG(10000))
model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                   , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
                   , method=CG(10000))

model_re_nosig_np<-gamlss(random_variable~as.factor(time==1)+random(as.factor(patient)),data=dataset,family=GA(),method=CG(10000),mixture="np")
help(gamlss.random)
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

### Investigating



########### 3. GJRM fits #########

require(GJRM)
eq.mu.1 <- gamma_c_mu1~1
eq.mu.2 <- gamma_c_mu2~1
fl <- list(eq.mu.1, eq.mu.2)
model_copula<-gjrm(fl, margins = c("GA" , "GA") , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
model_copula_n<-gjrm(fl, margins = c("GA" , "GA") , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")



#AIC for copula

#model_glm$aic
#summary(model_glm)
#2*4-2*model_copula_n$logLik
#2*6-2*model_copula$logLik
#summary(model_copula)
#model_re_nosig$aic

#2*6-2*logLik(model_re_nosig)
#summary(model_re_nosig)
#model_re$aic


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

########### 4. Combining results #########

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


########## 5. Investigating GLMM fits

library(gamlss)
library(mgcv)
set.seed(1000)

####Benchmark Gamma GLM
gamlss_fit_0<-gamlss(formula=random_variable~as.factor(time==1)
       , sigma.formula=~as.factor(time==1), data=dataset, family=GA()
       , method=CG(10000))
plot(gamlss_fit_0)

####Gamma + Random Effect
gamlss_fit_GA<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=GA()
                      #, method=CG(10000)
                      )
plot(gamlss_fit_GA)
####Gamma + Random Effect
#gamlss_fit_GA_s<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
#                      , sigma.formula=~as.factor(time==1)
#                      , data=dataset, family=GA(),method=CG(1000))
#plot(gamlss_fit_GA)

#####Look at the fitted random effect
#So - 1. the arndom effect doesn't look normal, 2. 
par(mfrow=c(2,2))
hist(log(gamma_c_mu1),xlim=c(-8,0),main="Margin 1 (log)",breaks =20)
hist(log(gamma_c_mu2),xlim=c(-8,0),main="Margin 2 (log)",breaks =20)
hist(gamlss_fit_GA$mu.coefSmo[[1]]$coef,xlim=c(-4,4),main="Random effect",breaks =20)
hist(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef,xlim=c(-8,0),main="Margin 2 minus Random effect",breaks =20)

####Looking at the 
plot.new()
par(mfrow=c(2,2))
histDist(gamma_c_mu1,main = "Original margin 1 (gamma)",family=GA(),ylim=c(0,14),nbins=100,xlim=c(0,.5))
histDist(gamma_c_mu2,main = "Original margin 2 (gamma)",family=GA(),ylim=c(0,14),nbins=50,xlim=c(0,.5))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),ylim=c(0,14),nbins=400,xlim=c(0,.5),main="Margin 1 minus random effect")
histDist(exp(log(gamma_c_mu2)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),ylim=c(0,14),nbins=600,xlim=c(0,.5),main="Margin 2 minus random effect")
###So adding RE makes it MORE skewed

#Not a substantial improvement
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),nbins=300,xlim=c(0,.5),main="GA",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GB2(),nbins=300,xlim=c(0,.5),main="GB2",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=BCTo(),nbins=300,xlim=c(0,.5),main="BCTo",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GP(),nbins=300,xlim=c(0,.5),main="GP",ylim=c(0,14))

histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GA(),nbins=300,xlim=c(0,.5),main="GA",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GB2(),nbins=300,xlim=c(0,.5),main="LOGNO",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=BCTo(),nbins=300,xlim=c(0,.5),main="LOGNO2",ylim=c(0,14))
histDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),family=GP(),nbins=300,xlim=c(0,.5),main="GP",ylim=c(0,14))


fitDist(exp(log(gamma_c_mu2)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),type="realplus")$fits
#GB2      BCTo        GP   PARETO2  PARETO2o     BCPEo     BCCGo        GG    LOGNO2     LOGNO      WEI3      WEI2       WEI 
#-2555.342 -2554.560 -2553.936 -2553.936 -2553.936 -2553.029 -2550.251 -2549.281 -2528.517 -2528.517 -2457.936 -2457.936 -2457.936 
#GIG        GA       EXP        IG    IGAMMA 
#-2439.929 -2382.074 -2285.360 -2054.721 -1926.681 

fitDist(exp(log(gamma_c_mu1)-gamlss_fit_GA$mu.coefSmo[[1]]$coef),type="realplus")$fits
#BCTo      BCPEo    PARETO2         GP   PARETO2o        GB2      BCCGo         GG     LOGNO2      LOGNO       WEI2       WEI3 
#-2244.0461 -2244.0391 -2243.6023 -2243.6023 -2243.6023 -2243.0053 -2238.9137 -2237.5020 -2215.2062 -2215.2062 -2150.5406 -2150.5406 
#WEI        GIG         GA        EXP     IGAMMA         IG 
#-2150.5406 -2076.3217 -2067.8069 -1910.3904  -967.5800  -790.0647 


#Won't coverge
#gamlss_fit_GB2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
#                      , data=dataset, family=GB2()
#                      , method=CG(1000)
#                      )
#
plot(gamlss_fit_GB2)

#Plots look ridiculous but global deviance is better 
gamlss_fit_BCTo<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=BCTo()
                      , method=RS(100)
                      )

plot(gamlss_fit_BCTo)

##This is the only 2-parameter distribution with a better AIC fit
gamlss_fit_GP<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                        , data=dataset, family=GP()
                        #,method=RS(1000)
                        , method = CG(10000)
                        #, method=mixed(10,1000)
)

plot(gamlss_fit_GP) ###Fit looks ridiculous

gamlss_fit_PARETO2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                      , data=dataset, family=PARETO2()
                      #,method=RS(1000)
                      #, method = CG(10000)
                      #, method=mixed(10,1000)
)

plot(gamlss_fit_PARETO2) ###Fit looks ridiculous

gamlss_fit_LOGNO2<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                           , data=dataset, family=LOGNO2()
                           #,method=RS(1000)
                           #, method = CG(10000)
                           #, method=mixed(10,1000)
)

plot(gamlss_fit_LOGNO2) ###Fit looks ridiculous

gamlss_fit_LOGNO<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                          , data=dataset, family=LOGNO()
                          #,method=RS(1000)
                          #, method = CG(10000)
                          #, method=mixed(10,1000)
)

plot(gamlss_fit_LOGNO) ###Fit looks ridiculous

gamlss_fit_WEI<-gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient))
                         , data=dataset, family=WEI2()
                         #,method=RS(1000)
                         #, method = CG(10000)
                         #, method=mixed(10,1000)
)

plot(gamlss_fit_WEI) ###Fit looks ridiculous
summary(gamlss_fit_WEI)


gamlss_fit_GA$aic
gamlss_fit_BCTo$aic
gamlss_fit_GP$aic
gamlss_fit_PARETO2$aic
gamlss_fit_LOGNO2$aic
gamlss_fit_LOGNO$aic
gamlss_fit_WEI$aic


