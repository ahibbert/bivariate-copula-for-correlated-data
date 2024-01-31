library(MASS)
library(psych)
library(copula)
library(VineCopula)
library(gamlss)
library(moments)
require(ggpubr)
require(ggplot2)
require(dglm)
library(latex2exp)

plotSimBiGamma <- function(n,a,b,mu1,mu2,cbins,xlim,ylim,type)  {
  
  w<-rbeta(n,a,b) #Mean .5
  gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1) #Mean 6 * .5 = 3
  gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2) #Mean 12 * .5 = 6
  
  #cor(gamma_c_mu1,gamma_c_mu2,method="kendall")
  #skewness(gamma_c_mu1)
  #skewness(gamma_c_mu2)
  
  #patient<-as.factor(seq(1:n))
  #dataset<-as.data.frame(rbind(cbind(patient,gamma_c_mu1,0),cbind(patient,gamma_c_mu2,1)))
  #colnames(dataset)<-c("patient","random_variable","time")
  
  tau=cor(gamma_c_mu1,gamma_c_mu2,method="kendall")
  
  u<-0
  v<-0
  #u<-pgamma(gamma_c_mu1,shape=a+b,scale=1/mu1)
  #v<-pgamma(gamma_c_mu2,shape=a+b,scale=1/mu2)
  
  u<-pgamma(gamma_c_mu1,shape=fitdistr(gamma_c_mu1,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu1,"gamma")$estimate[2])
  v<-pgamma(gamma_c_mu2,shape=fitdistr(gamma_c_mu2,"gamma")$estimate[1],rate=fitdistr(gamma_c_mu2,"gamma")$estimate[2])
  
  if(type=="biv") {
    plot1<-ggplot(data=as.data.frame(cbind(gamma_c_mu1,gamma_c_mu2)),aes(x=gamma_c_mu1,y=gamma_c_mu2)) + 
      #geom_point(size=0.1,color="gray") + 
      geom_density_2d(contour_var="density",bins=cbins,color="black") + 
      scale_fill_brewer()  +
      #theme(panel.background = element_blank()) +
      xlim(0,xlim) +
      ylim(0,ylim) +
      labs(x = TeX("$Y_1$"), y=TeX("$Y_2$")
           #, title=TeX(paste("$\\mu_1=",(a/mu1),"\\ \\beta=",b,"\\ \\tau=",round(tau*10^2)/10^2))
           ,fill="density")
  }
  else
  {
    plot1<-ggplot(data=as.data.frame(cbind(u,v)),aes(x=u,y=v)) +
      #geom_point(size=0.25,color="black") + 
      geom_density_2d(contour_var="density",bins=cbins,color="black") + 
      scale_fill_brewer() +
      #theme(panel.background = element_blank()) +
      labs(x = TeX("$Y_1$ (Uniform)"), y=TeX("$Y_2$ (Uniform)")#, title=TeX(paste("$\\mu_1=",(a/mu1),"\\ \\beta=",b,"\\ \\tau=",round(tau*10^2)/10^2))
           ,fill="density")
  }
  #plot1<-persp(kde2d(gamma_c_mu1,gamma_c_mu2),main="Uniform transform of marginals",axes=T,scale=T,ticktype="detailed",theta=15,xlab="Y1",ylab="Y2",zlab="Density") #zlim=c(0,4) #,h=.4,n=65
  
  return(plot1)
}

set.seed(100)
options(scipen=999)

plot_list = list()

i=1; plot_list[[i]]<-plotSimBiGamma(a=2, b=.2,mu1=10, mu2=12, n=10000,cbins=12,xlim=.5,ylim=.5,type="biv");
i=2; plot_list[[i]]<-plotSimBiGamma(a=.5, b=.2,mu1=10, mu2=12, n=10000,cbins=12,xlim=.04,ylim=.04,type="biv");
i=3; plot_list[[i]]<-plotSimBiGamma(a=2, b=2,mu1=10, mu2=12, n=10000,cbins=12,xlim=.5,ylim=.5,type="biv");
i=4; plot_list[[i]]<-plotSimBiGamma(a=.5, b=2,mu1=10, mu2=12, n=10000,cbins=12,xlim=.04,ylim=.04,type="biv");

i=5; plot_list[[i]]<-plotSimBiGamma(a=2, b=.2,mu1=10, mu2=12, n=10000,cbins=11,xlim=.6,ylim=.3,type="unif");
i=6; plot_list[[i]]<-plotSimBiGamma(a=.5, b=.2,mu1=10, mu2=12, n=10000,cbins=11,xlim=.1,ylim=.05,type="unif");
i=7; plot_list[[i]]<-plotSimBiGamma(a=2, b=2,mu1=10, mu2=12, n=10000,cbins=11,xlim=.6,ylim=.3,type="unif");
i=8; plot_list[[i]]<-plotSimBiGamma(a=.5, b=2,mu1=10, mu2=12, n=10000,cbins=11,xlim=.01,ylim=.005,type="unif");


ggarrange(plotlist = plot_list[c(1,2,5,6,3,4,7,8)],nrow=2,ncol=4,
                        common.legend = TRUE, legend="bottom")

# library(grid)
# pl <- plot_list[c(1,2,5,6,3,4,7,8)]
 
 #N <- length(pl)
# nr <- 2
# nc <- 4
# 
 #combine <- rbind(tableGrob(t(c("\u03B2=0.2","\u03B2=2")), theme = ttheme_minimal(), rows = ""), 
 #                cbind(tableGrob(c("     \u03BC=0.2","     \u03BC=0.05","     \u03BC=0.05","     \u03BC=0.2"), theme = ttheme_minimal()), 
#                       arrangeGrob(grobs = pl),  size = "last"), size = "last")
# grid.newpage()
# grid.draw(combine)
 
library(gridExtra)
# Create list of plots
set.seed(0)
pl = plot_list[c(1,3,2,4,5,7,6,8)]
# Create row and column titles
col.titles = c("     \u03C3\u00B2=0.5","     \u03C3\u00B2=2","     \u03C3\u00B2=0.5","     \u03C3\u00B2=2")
row.titles = c("\u03B2=0.2","\u03B2=2")

grid.newpage()

nr=2
# Add row titles
pl[1:nr] = lapply(1:nr, function(i) arrangeGrob(pl[[i]], left=row.titles[i]))

# Add column titles and lay out plots
plot=grid.arrange(grobs=lapply(c(1,3,5,7), function(i) {
  arrangeGrob(grobs=pl[i:(i+1)], top=col.titles[(i+1)/2], ncol=1)
}), ncol=4)

#ggsave(file="bivariate_distribution_contours.png",plot,width=12,height=6,dpi=300) 



