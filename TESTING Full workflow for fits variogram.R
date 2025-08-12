source("common_functions.R")
dist="NO";a=1;b=1;c=.75;mu1=1;mu2=2;x1=1;x2=1;n=100

#fits=fitBivModels_Bt_withCov(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE,cv=FALSE)

#dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
#fits=fitBivModels_Bt_withCov(      dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
#eval=evaluateModels(fits,vg_sims=100)

library(callr)
outer_sims=10
eval_outer=list()
for (i in 1:outer_sims) {
  print(i)
  eval <- callr::r(
    func = function(dist,a,b,c,mu1,mu2,x1,x2,n) {
      source("common_functions.R")
      dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
      #Set dataset as global as old S3 functions in gamlss and sometimes GJRM break without it
      assign("dataset", dataset, envir = .GlobalEnv)
      fits=fitBivModels_Bt_withCov(  dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
      eval=evaluateModels(fits,vg_sims=100)
      return(eval)
    },
    args = list(dist = dist, a = a, b = b, c = c,
                mu1 = mu1, mu2 = mu2, x1 = x1, x2 = x2, n = n)
  )
  eval_outer[[i]]=eval
}

score_items=list()
score_item_names=c("vs2_wt"   ,      "vs2",            "es","vs1"          ,  "vs2_wt_coronly")
for (i in 1:length(eval_outer)) {
  item_names=names(eval_outer[[i]])
  for(item in score_item_names) {
    if(i==1) {
      score_items[[item]]=eval_outer[[i]][[item]]  
    } else {
      score_items[[item]]=rbind(score_items[[item]],eval_outer[[i]][[item]])
    }
  }
}



#save(eval_outer,file="Data/eval_outer_NO_1_1_75_1_2_1_1_100.RData")
