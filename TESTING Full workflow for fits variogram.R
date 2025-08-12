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
score_item_names=c("vs2_wt"   ,      "vs2",            "es","vs1"          ,  "vs2_wt_coronly","logliks")
for (i in 1:length(eval_outer)) {
  #item_names=names(eval_outer[[i]])
  for(item in score_item_names) {
    new_item=if(item=="logliks"){eval_outer[[i]][[item]][,1]}else{eval_outer[[i]][[item]]}
     
    if(i==1) {
      score_items[[item]]=new_item
    } else {
      score_items[[item]]=rbind(score_items[[item]],new_item)
    }
  }
}

par_estimates=list()
par_item_names=c("coefficients","ses","sigmas","correlations")

for(i in 1:11) {
  par_estimates[[i]]=matrix(nrow=length(eval_outer),ncol=16)
}

for(i in 1:length(eval_outer)) {
  z=1
  for(item in par_item_names) {
    new_item=data.frame(eval_outer[[i]][[item]])
    for (j in 1:ncol((new_item))) {
      par_estimates[[z]][i,]=new_item[,j]
      z=z+1; 
    }
  }
}
names(par_estimates)=c("mu1","mu2","x1","x2","mu1_se","mu2_se","x1_se","x2_se","sigma1","sigma2","corr")
for(item in names(par_estimates)) {
  colnames(par_estimates[[item]])=c("glm", "gee", "re_nosig", "re_np", "lme4", "gamm",
                                    "cop", "cop_n", "cop_j", "cop_g", "cop_f",
                                    "cop_amh", "cop_fgm", "cop_pl", "cop_h", "cop_t")
}

