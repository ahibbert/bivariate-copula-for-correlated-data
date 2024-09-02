library(MASS)
library(psych)
library(copula)
library(VineCopula)
library(gamlss)
library(moments)
require(ggpubr)
require(ggplot2)
require(dglm)

library(bigamma)

set.seed(1000)
options(scipen=999)

###Bivariate gamma of Mathai and Moschopoulos

n=2000; alpha_0=10; alpha_1= .5; alpha_2=.75; beta_0=1; beta_1=4;beta_2=5; ## >0

v_0=rgamma(n,shape=alpha_0,scale=beta_0)
v_1=rgamma(n,shape=alpha_1,scale=beta_1)
v_2=rgamma(n,shape=alpha_2,scale=beta_2)

x = (beta_0/beta_1)*v_0 + v_1
y = (beta_0/beta_2)*v_0 + v_2

par(mfrow=c(3,1))
hist(v_0)
hist(x)
hist(y)

###Bivariate SAT Gamma




gamma_c_mu1 <- x
gamma_c_mu2 <- y

#############REFERENCE#################
#a=0.25; b=1.75 #.25, .75 the most extreme
#mu1=1; mu2=2; n=2000 #100,500,1000,5000,10000

#w<-rbeta(n,a,b) #Mean .5
#gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1) #Mean 6 * .5 = 3
#gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2) #Mean 12 * .5 = 6

patient<-as.factor(seq(1:n))
dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
colnames(dataset)<-c("patient","random_variable","time")

u<-0
v<-0
#u<-pgamma(gamma_c_mu1,shape=a+b,scale=1/mu1)
#v<-pgamma(gamma_c_mu2,shape=a+b,scale=1/mu2)

u<-pgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])
v<-pgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])

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
xx<-ggplot(data = as.data.frame(gamma_c_mu2)) +
  geom_histogram(data = as.data.frame(gamma_c_mu2), aes(x=gamma_c_mu2, y=..density..),bins=50) +
  geom_line(aes(lty = 'fitted gamma',x=gamma_c_mu2, y=dgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])), color="blue", size = .75) +
  #ylim(0,10) +
  labs(x="time 2 margin",title="time 2 margin") +
  #xlim(0,1) +
  theme(legend.position = c(0.75, .92),legend.key = element_rect(fill = "transparent"),legend.title = element_blank(), legend.background = element_blank())
c<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
  geom_point(size=0.5,color="black") + 
  labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2") +
  geom_smooth(method="loess", level=.99) #+
  #xlim(0,2) +
  #ylim(0,1) 
d<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=15,color="black") + 
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="margin 1 v margin 2 (unif. transform)",fill="density")
e<-ggplot(data=as.data.frame(fittedClayton),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d(contour_var="density",bins=15,color="black") +
  scale_fill_brewer() +
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted clayton copula")
f<-ggplot(data=as.data.frame(fittedTDist),aes(x=V1,y=V2)) + 
  #geom_point(size=0.25,color="black") + 
  geom_density_2d( contour_var="density",bins=15,color="black") + 
  labs(x = "time 1 margin", y="time 2 margin", title="simulated fitted normal copula")
ggarrange(z,xx,c,nrow=1)
##ggsave(file="example_bias_case_margin_plots.jpeg",last_plot(),width=14,height=4,dpi=300)
ggarrange(d,e,f,common.legend = TRUE,nrow=1,legend="right")
##ggsave(file="example_bias_case_contour_plots.jpeg",last_plot(),width=14,height=4,dpi=300)
