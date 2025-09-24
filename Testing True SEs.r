source("common_functions.R")


    true_out=calcTrueCovariateValues(
  n=1000,
  a=1,
  b=.3,
  c=1,
  mu1=10,
  mu2=12,
  dist="GA",
  x1=1,
  x2=.01
    )
    true_out$par

ses10000_g=get_true_ses(
  n=10000,
  a=.2,
  b=.34,
  c=1,
  mu1=10,
  mu2=12,
  dist="GA",
  x1=1,
  x2=.01
)

ses100_g=get_true_ses(
  n=100,
  a=.2,
  b=.3,
  c=1,
  mu1=10,
  mu2=12,
  dist="GA",
  x1=1,
  x2=.01,
  sims=100
)

rbind(ses1000_g,ses100_g,ses10000_g)/sqrt(1000)

ses1000_1000=get_true_ses(  n=1000,  a=1,  b=1,  c=.1,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=1000)
ses1000_100=get_true_ses(  n=1000,  a=1,  b=1,  c=.05,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=100)
ses100_1000=get_true_ses(  n=100,  a=1,  b=1,  c=.1,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=1000)
ses100_100=get_true_ses(  n=100,  a=1,  b=1,  c=.05,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=100)

rbind(ses1000_1000*1000, ses1000_100*1000, ses100_1000*100, ses100_100*100)/sqrt(1000)


source("common_functions.R")
n=1000; a=1; b=1; c=.5; mu1=1; mu2=2; dist="NO"; x1=1; x2=.01
dataset = generateBivDist_withCov(n=n, a=a, b=b, c=c, mu1=mu1, mu2=mu2, dist=dist, x1=x1, x2=x2)
fits = fitBivModels_Bt_withCov(dataset=dataset, dist=dist, include="ALL",
                                         a=a, b=b, c=c, mu1=mu1, mu2=mu2,
                                         calc_actuals=FALSE, cv=FALSE,bt_mode = TRUE)
fits$coefficients
eval = evaluateModels(fits, vg_sims=100,bt_mode=TRUE)
eval

# Short script to analyze true_ses_matrix by distribution
dists <- substr(rownames(true_ses_matrix), 1, 2)
  
# Calculate means by distribution
result <- aggregate(true_ses_matrix, by = list(dists), FUN = mean, na.rm = TRUE)
names(result)[1] <- "distribution"

rownames(result) <- result$distribution
result=result[c(3,1,4,2),c(2:ncol(result))]
true=c(.1,.05,.05,.1)

result_adj=sqrt(result)
cbind(true,round(result_adj,4))