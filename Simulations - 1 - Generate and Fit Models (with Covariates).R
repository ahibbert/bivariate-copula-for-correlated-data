# Set up a list of input configurations you wish to run
input_list <- list(
  list(Bt_mode=FALSE , a=1    , b=2        ,c=.75, mu1=1          , mu2=2                  , n=1000,dist="NO",x1=c(1,3),x2=1),
  list(Bt_mode=FALSE , a=.5*1:5    , b=.5*1:5        ,c=c(.1,.2,.3,.4,.5,.6,.7,.8,.9), mu1=1          , mu2=2                  , n=1000,dist="NO",x1=1,x2=1),
  list(Bt_mode=FALSE , a=NA        , b=c(.2,.5,1,2,5),c=c(.2,.5,1,2,5)               , mu1=c(.5,1,2,5), mu2=c(.5,1,2,5)        , n=1000,dist="PO",x1=1,x2=1),
  list(Bt_mode=FALSE , a=.1+.1*1:20, b=.1+.1*1:20    ,c=NA                           , mu1=10         , mu2=12                 , n=1000,dist="GA",x1=1,x2=1),
  list(Bt_mode=FALSE , a=NA        , b=NA            ,c=c(.1,.25,.5,.75,.9), mu1=c(.1,.25,.5,.75,.9)  , mu2=c(.1,.25,.5,.75,.9), n=1000,dist="LO",x1=1,x2=1)
)
list2env(input_list[[1]], envir = .GlobalEnv)

options(scipen=999)
source("common_functions.R")
results<-list()
datasets<-list()
#list2env(input_list[[1]], envir = .GlobalEnv)

model_struct_FUN=fitBivModels_Bt_withCov

#Code to iterate through various shapes of the bivariate distribution and fit the non-GJRM models
i=1; j=1; k=1; l=1; z=1
start=Sys.time()
iterations=length(a)*length(b)*length(mu1)*length(mu2)*length(c)*length(x1)*length(x2)*2 #times 2 for types

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
                results[[z]]$actuals=c(n,a[i],b[j],c[m],mu1[k],mu2[l],x1[o],x2[p])
                print(c(n,a[i],b[j],c[m],mu1[k],mu2[l],x1[o],x2[p]))
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
  results_combined[[i]]=list()
  for (j in 1:length(results[[i]])) {
    results_combined[[i]][[j]]=tryCatch( {rbind(results[[i]][[j]],results[[(length(results)/2)+i]][[j]])}, error = function(e) {NA} )  
  }
  
  na_check[i]=(is.null(nrow(results_combined[[i]])))
}
#results_combined <- results_combined[na_check==FALSE]

##### SAVE LOCATION
save(results_combined,file=paste("Data/results_combined",dist,"_",n,"_",Sys.Date(),".RData",sep=""))


