############## 0. Required functions ################### 

options(scipen=999)
source("common_functions.R")

############## 1. Run simulations for non-GJRM models########################
results<-list()
datasets<-list()
#a=.1+.1*1:20; b=.1+.1*1:20; c=NA; mu1=10; mu2=12; n=1000; dist="GA"
#a=.5*1:5; b=.5*1:5;c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9); mu1=1; mu2=2; n=1000;dist="NO"
a=NA; b=NA;c=c(.2,.5,1,2,5); mu1=c(.2,.5,1,2,5); mu2=c(.2,.5,1,2,5); n=1000;dist="PO"

#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1;
start=Sys.time()
iterations=length(a)*length(b)*length(mu1)*length(mu2)*length(c)*2 #times 2 for types

for (type in c("non-GJRM","GJRM")) {
  for (i in 1:length(a)) {
    for (j in 1:length(b)) {
      for (k in 1:length(mu1)) {
        for (l in 1:length(mu2)) {
          for (m in 1:length(c)) {
            set.seed(1000)
            datasets[[z]]<-generateBivDist(n,a[i],b[j],c[m],mu1[k],mu2[l],dist)
            tryCatch( {results[[z]] <- rbind(fitBivModels(datasets[[z]],dist,include=type,a[i],b[j],c[m],mu1[k],mu2[l]), c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA))}
                      , error = function(e) {results[[z]]<- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA), c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA))})
            print(c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA))
            print(c(z
                    ,z/(iterations) 
                    , (Sys.time()-start)
                    , (Sys.time()-start) / (z/iterations))
            )
            z = z + 1
          }
        }
      }
    }
  }
}

#only needed for NB case as lme4 and gamlss (4) have some failures occassionally
na_check=rep(TRUE,length(results)/2)
results_combined <- list()
for(i in 1:(length(results)/2)) {
  print(i); 
  results_combined[[i]]=tryCatch( {rbind(results[[i]][1:(nrow(results[[i]])-2),],results[[(length(results)/2)+i]])}, error = function(e) {c(NA,NA,NA,NA,NA,NA,NA,NA)} )
  na_check[i]=(is.null(nrow(results_combined[[i]])))
}
results_combined <- results_combined[na_check==FALSE]

save(results_combined,file=paste("results_combined_",dist,"_",n,"_",Sys.Date(),".RData",sep=""))

