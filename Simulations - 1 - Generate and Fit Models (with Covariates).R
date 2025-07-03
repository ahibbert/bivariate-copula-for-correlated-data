options(scipen=999)
source("common_functions.R")
results<-list()
datasets<-list()
#list2env(input_list[[1]], envir = .GlobalEnv)

model_struct_FUN=fitBivModels_Bt_withCov

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
            for (o in 1:length(x1)) {
              for (p in 1:length(x2)) {
                set.seed(1000)
                if(type == "non-GJRM") {datasets[[d]]<-generateBivDist_withCov(n,a[i],b[j],c[m],mu1[k],mu2[l],dist,x1[o],x2[p])}
                tryCatch( {results[[z]] <- model_struct_FUN(datasets[[d]],dist,include=type,a[i],b[j],c[m],mu1[k],mu2[l])}
                          , error = function(e) {results[[z]]<- NA})
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

##### SAVE LOCATION
if(Bt_mode==TRUE) {
  save(results_combined,file=paste("Data/results_combined_B1_Bt_",dist,"_",n,"_",Sys.Date(),".RData",sep=""))
} else {
  save(results_combined,file=paste("Data/results_combined_B1_B2_",dist,"_",n,"_",Sys.Date(),".RData",sep=""))
}

