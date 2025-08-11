source("common_functions.R")
dist="NO";a=1;b=1;c=.75;mu1=1;mu2=2;x1=1;x2=1;n=1000

#fits=fitBivModels_Bt_withCov(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE,cv=FALSE)

dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
fits=fitBivModels_Bt_withCov(      dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
eval=evaluateModels(fits,model_list=c("glm","gee","cop_n"),vg_sims=100)

###OK Now we need to extract the info needed for the next function from fits.

library(callr)
eval <- callr::r(
  func = function(dist,a,b,c,mu1,mu2,x1,x2,n) {
    source("common_functions.R")
    dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
    #Set dataset as global as old S3 functions in gamlss and sometimes GJRM break without it
    assign("dataset", dataset, envir = .GlobalEnv)
    fits=fitBivModels_Bt_withCov(  dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
    eval=evaluateModels(fits,model_list=c("glm","gee","cop_n"),vg_sims=100)
    return(eval)
  },
  args = list(dist = dist, a = a, b = b, c = c,
              mu1 = mu1, mu2 = mu2, x1 = x1, x2 = x2, n = n)
)
eval
