###numerical derivatives for the poisson

poisson_dist <- function(par, n=1) {
  #Compound multiple poisson of Stein & Juritz, 1987
  a=par[1];b=par[2];c=a=par[3];mu1=par[4];mu2=par[5]
  
  mixing_dist<-rgamma(n,shape=c,scale=b)
  
  margin_1=vector(length = n) 
  margin_2=vector(length = n) 
  for (i in 1:n) {
    margin_1[i]=rpois(1,mu1*mixing_dist[i])
  }
  for (i in 1:n) {
    margin_2[i]=rpois(1,mu2*mixing_dist[i])
  }
  return(cbind(margin_1,margin_2))
}

a=NA; b=c(.2,.5,1,2,5);c=c(.2,.5,1,2,5); mu1=c(.5,1,2,5); mu2=c(.5,1,2,5); n=1000;dist="PO"
#a=NA; b=c(.2,.5,1,2,5);c=c(.2,.5,1,2,5); mu1=c(.5,1,2,5); mu2=c(.5,1,2,5); n=1000;dist="PO"

#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1
start=Sys.time()
iterations=length(a)*length(b)*length(mu1)*length(mu2)*length(c)*2 #times 2 for types
dataset=list()
means=list()
mle=matrix(NA,ncol=2,nrow=0)
colnames=c("mean","var")
nb_par=matrix(NA,ncol=5,nrow=0)

for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    for (k in 1:length(mu1)) {
      for (l in 1:length(mu2)) {
        for (m in 1:length(c)) {
          dataset[[z]]=list()
          means[[z]]=matrix(NA,ncol=1,nrow=0)
          for (i in 1:1000) {
            dataset[[z]][[i]]=(poisson_dist(par=c(a[i],b[j],c[m],mu1[k],mu2[l]),n=1000))
            means[[z]]=rbind(means[[z]],log(mean(dataset[[z]][[i]][,2])/mean(dataset[[z]][[i]][,1])))
          }
          mle=rbind(mle,c(mean(means[[z]]),var(means[[z]])))
          nb_par=rbind(nb_par,c(a[i],b[j],c[m],mu1[k],mu2[l]))
        }
      }
    }
  }
}

nb_mle<-cbind(nb_par,mle)
colnames(nb_mle)=c("a","b","c","mu1","mu2","mean_Bt","var_Bt")

save(nb_mle,file="nb_mle")