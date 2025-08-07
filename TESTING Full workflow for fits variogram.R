source("common_functions.R")
dist="NO";a=1;b=1;c=.75;mu1=1;mu2=2;x1=1;x2=1;n=1000

#fits=fitBivModels_Bt_withCov(dataset,dist,include="ALL",a,b,c,mu1,mu2,calc_actuals=TRUE,cv=FALSE)

dataset=generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
fits=fitBivModels_Bt_withCov(      dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
eval=evaluateModels(fits,model_list=c("glm","gee","cop_n"))

###OK Now we need to extract the info needed for the next function from fits.



library(callr)
fits <- callr::r(
  func = function(dist,a,b,c,mu1,mu2,x1,x2,n) {
    source("common_functions.R")
    dataset <- generateBivDist_withCov(n=n,a=a,b=b,c=c,mu1=mu1,mu2=mu2,dist=dist,x1=x1,x2=x2)
    # Call the function to fit the models
    fitBivModels_Bt_withCov(  dataset=dataset,      dist=dist,      include = "ALL",      a=a,b=b,c=c,mu1=mu1,mu2=mu2,      calc_actuals = FALSE,      cv = FALSE    )
  },
  args = list(dist = dist, a = a, b = b, c = c,
              mu1 = mu1, mu2 = mu2, x1 = x1, x2 = x2, n = n)
)

### NEED ALL PARAMETERS - SIGMAS, CORRELATION PARAMETERS etc as we need to simulate it for random effect models
# saveRDS(fits, file = "fits.rds")
