############# Overview of file #####################

# Purpose of this file is to fit a series of different models to a simulated datasets and collate results
# Code iterates over parameter values c(a,b,c,mu1,mu2) for sample size n and saves the output results_combined to a file
# Results_combined includes parameter estimates for the mean at time 1 and 2 and their standard error + log likelihood
# Change 1) INPUTS to change simulations. Rest is automated. Define Bt_mode=TRUE for model with B1 and Bt, or model for just means if Bt_mode=FALSE
# Otherwise input vectors of parameters to be iterated over and the chosen distribution to simulate as dist="NO","PO","GA" or "LO"

#### 1) INPUTS ####

#Bt_mode=TRUE
#Final paper parameters
a=.5*1:5; b=.5*1:5;c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9); mu1=1; mu2=2; n=1000;dist="NO"
a=NA; b=c(.2,.5,1,2,5);c=c(.2,.5,1,2,5); mu1=c(.5,1,2,5); mu2=c(.5,1,2,5); n=1000;dist="PO"
a=.1+.1*1:20; b=.1+.1*1:20; c=NA; mu1=10; mu2=12; n=1000; dist="GA"
a=NA; b=NA;c=c(.1,.25,.5,.75,.9); mu1=c(.1,.25,.5,.75,.9); mu2=c(.1,.25,.5,.75,.9); n=1000;dist="LO"

#### 2) RUN SIMULATIONS ####

### THE REST IS AUTOMATIC: Saves into paste("Data/results_combined_B1_B2_",dist,"_",n,"_",Sys.Date(),".RData",sep="") ###

############## 1. Run simulations for non-GJRM models########################
options(scipen=999)
source("common_functions.R")

results<-list()
datasets<-list()

if(Bt_mode==TRUE) {
  model_struct_FUN=fitBivModels_Bt  
} else {
  model_struct_FUN=fitBivModels
}


#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1
start=Sys.time()
iterations=length(a)*length(b)*length(mu1)*length(mu2)*length(c)*2 #times 2 for types

for (type in c("non-GJRM","GJRM")) {
  d=1
  for (i in 1:length(a)) {
    for (j in 1:length(b)) {
      for (k in 1:length(mu1)) {
        for (l in 1:length(mu2)) {
          for (m in 1:length(c)) {
            set.seed(1000)
            # Generate a dataset for fit i
            datasets[[d]]<-generateBivDist(n,a[i],b[j],c[m],mu1[k],mu2[l],dist)
            # Fit all models
            tryCatch( {results[[z]] <- rbind(model_struct_FUN(datasets[[d]],dist,include=type,a[i],b[j],c[m],mu1[k],mu2[l]), c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA))}
                      , error = function(e) {results[[z]]<- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA), c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA))})
            # Print progress
            print(c(n,a[i],b[j],c[m],mu1[k],mu2[l],NA,NA,type))
            print(c(z
                    ,z/(iterations) 
                    , (Sys.time()-start)
                    , (Sys.time()-start) / (z/iterations))
            )
            z = z + 1; d=d+1
          }
        }
      }
    }
  }
}

# Only needed for NB case as lme4 and gamlss (4) have some failures occasionally
na_check=rep(TRUE,length(results)/2)
results_combined <- list()
for(i in 1:(length(results)/2)) {
  print(i); 
  results_combined[[i]]=tryCatch( {rbind(results[[i]][1:(nrow(results[[i]])-2),],results[[(length(results)/2)+i]])}, error = function(e) {c(NA,NA,NA,NA,NA,NA,NA,NA)} )
  na_check[i]=(is.null(nrow(results_combined[[i]])))
}
results_combined <- results_combined[na_check==FALSE]

##### SAVE LOCATION
if(Bt_mode==TRUE) {
  save(results_combined,file=paste("Data/results_combined_B1_Bt_",dist,"_",n,"_",Sys.Date(),".RData",sep=""))
} else {
  save(results_combined,file=paste("Data/results_combined_B1_B2_",dist,"_",n,"_",Sys.Date(),".RData",sep=""))
}

