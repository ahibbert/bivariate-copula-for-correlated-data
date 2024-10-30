###numerical derivatives for the poisson

source("link_functions.R")

LO_dist <- function(par, n=1) {
  #Compound multiple poisson of Stein & Juritz, 1987
  a=par[1];b=par[2];c=a=par[3];mu1=par[4];mu2=par[5]
  
  require(MASS)
  a=1;b=1
  normData<-mvrnorm(n,mu=c(0,0),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
  margin_1<-as.numeric(pnorm(normData[,1])<=mu1)
  margin_2<-as.numeric(pnorm(normData[,2])<=mu2)
  
  return(cbind(margin_1,margin_2))
}
a=NA; b=NA;c=c(.1,.25,.5,.75,.9); mu1=c(.1,.25,.5,.75,.9); mu2=c(.1,.25,.5,.75,.9); n=1000;dist="LO"

#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1
start=Sys.time()
iterations=length(a)*length(b)*length(mu1)*length(mu2)*length(c) #times 2 for types
dataset=list()
means=list()
mle=matrix(NA,ncol=6,nrow=0)
colnames=c("mean","var")
nb_par=matrix(NA,ncol=5,nrow=0)

for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    for (k in 1:length(mu1)) {
      for (l in 1:length(mu2)) {
        for (m in 1:length(c)) {
          start=Sys.time()
          dataset[[z]]=list()
          means[[z]]=matrix(NA,ncol=3,nrow=0)
          for (i in 1:1000) {
            dataset[[z]][[i]]=(LO_dist(par=c(a[i],b[j],c[m],mu1[k],mu2[l]),n=1000))
            means[[z]]=rbind((means[[z]]),c(logit(mean(dataset[[z]][[i]][,2])),logit(mean(dataset[[z]][[i]][,1])),logit(mean(dataset[[z]][[i]][,1]))-logit(mean(dataset[[z]][[i]][,2]))))
          }
          mle=rbind(mle,c(colMeans(means[[z]]),var(means[[z]][,1]),var(means[[z]][,2]),var(means[[z]][,3])))
          nb_par=rbind(nb_par,c(a[i],b[j],c[m],mu1[k],mu2[l]))
          print(c(z
                  ,z/(iterations) 
                  , (Sys.time()-start)
                  , (Sys.time()-start) / (z/iterations)))
          
          z=z+1
        }
      }
    }
  }
}

lo_mle<-cbind(nb_par,mle)
colnames(lo_mle)=c("a","b","c","mu1","mu2","mean_B1","mean_B2","mean_Bt","var_B1","var_B2","var_Bt")

save(lo_mle,file="Data/lo_mle")