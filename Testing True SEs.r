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
  x2=.01,
  )
  true_out$par

ses10000_g=get_true_ses(
  n=10000,
  a=.2,
  b=.34,
  c=1,
  mu1=10,
  mu2=12,
  dist="PO",
  x1=1,
  x2=.01
)

ses1000=get_true_ses(
  n=1000,
  a=1,
  b=1,
  c=.79,
  mu1=1,
  mu2=2,
  dist="NO",
  x1=1,
  x2=.01,
  sims=100
)

ses100=get_true_ses(
  n=100,
  a=1,
  b=1,
  c=.79,
  mu1=1,
  mu2=2,
  dist="NO",
  x1=1,
  x2=.01,
  sims=100
)

rbind(ses1000_g,ses100_g,ses10000_g)/sqrt(1000)

ses1000_1000=get_true_ses(  n=1000,  a=1,  b=1,  c=.1,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=1000)
ses1000_100=get_true_ses(  n=1000,  a=1,  b=1,  c=.05,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=100)
ses100_1000=get_true_ses(  n=100,  a=1,  b=1,  c=.1,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=1000)
ses100_100=get_true_ses(  n=100,  a=1,  b=1,  c=.05,  mu1=1,  mu2=1,  dist="NO",  x1=1,  x2=.01,  sims=100)

rbind(ses1000_1000, ses1000_100, ses100_1000/10, ses100_100/10)


ses1000_1000=get_true_ses(  n=1000,  a=1,  b=2,  c=2,  mu1=2,  mu2=2,  dist="PO",  x1=1,  x2=.01,  sims=1000)
ses1000_100=get_true_ses(  n=1000,  a=1,  b=2,  c=2,  mu1=2,  mu2=2,  dist="PO",  x1=1,  x2=.01,  sims=100)
ses100_1000=get_true_ses(  n=100,  a=1,  b=2,  c=2,  mu1=2,  mu2=2,  dist="PO",  x1=1,  x2=.01,  sims=1000)
ses100_100=get_true_ses(  n=100,  a=1,  b=2,  c=2,  mu1=2,  mu2=2,  dist="PO",  x1=1,  x2=.01,  sims=100)

#Currently running CG -> can do SANN next | 
source("common_functions.R")
rm(ses_no, ses_ga, ses_po, ses_lo)
ses_no=get_true_ses(  n=1000,  a=1,  b=1,  c=.5,  mu1=1,  mu2=2,  dist="NO",  x1=1,  x2=.01,  sims=100,debug_mode=TRUE)
ses_ga=get_true_ses(  n=1000,  a=.2,  b=1,  c=2,  mu1=10,  mu2=12,  dist="GA",  x1=1,  x2=.01,  sims=100,debug_mode=TRUE)
ses_po=get_true_ses(  n=1000,  a=5,  b=5,  c=5,  mu1=5,  mu2=5,  dist="PO",  x1=1,  x2=.01,  sims=100,debug_mode=TRUE)
ses_lo=get_true_ses(  n=1000,  a=1,  b=2,  c=.5,  mu1=.4,  mu2=.6,  dist="LO",  x1=1,  x2=.01,  sims=100,debug_mode=TRUE)

out=rbind(ses_no$ses, ses_ga$ses, ses_po$ses, ses_lo$ses)
rownames(out)=c("NO","GA","PO","LO")
sqrt(out)

par(mfrow=c(4,2))
hist(ses_no[[3]][,1],main="mu1")
hist(ses_ga[[3]][,1])
hist(ses_po[[3]][,1],breaks=20)
hist(ses_lo[[3]][,1])
hist(ses_no[[3]][,3])
hist(ses_ga[[3]][,3])
hist(ses_po[[3]][,3])
hist(ses_lo[[3]][,3])





)

sqrt(rbind(ses1000_1000, ses1000_100, ses100_1000, ses100_100))


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


#######################TESTING PO
source("common_functions.R")

#po_out=cbind(true_params_matrix,skew_matrix)[true_params_matrix$dist=="PO",]
#po_out[nrow(po_out),]

ses1000_po_test=get_true_ses(
  n=1000,
  a=NA,
  b=.5,
  c=.5,
  mu1=1,
  mu2=1,
  dist="PO",
  x1=1,
  x2=.01001,
  debug_mode=TRUE
)
  true_out=calcTrueCovariateValues(
  n=1000,
  a=NA,
  b=.5,
  c=.5,
  mu1=1,
  mu2=1,
  dist="PO",
  x1=1,
  x2=.0101,
  )
  true_out$par


(sqrt(cov(ses1000_po_test[[3]])))/sqrt(1000)
test=generateBivDist_withCov(n=1000, a=NA, b=5, c=5, mu1=5, mu2=5, dist="PO", x1=1, x2=.0101)
log(mean(test$random_variable[test$time==0]))


par(mfrow=c(2,2))
hist(ses1000_po_test[[3]][,1],breaks=30,main="mu1")
hist(ses1000_po_test[[3]][,2],breaks=30,main="mu2")
hist(ses1000_po_test[[3]][,3],breaks=30,main="x1")
hist(ses1000_po_test[[3]][,4],breaks=30,main="x2")

mean(ses1000_po_test[[3]][,1])
median(mean(ses1000_po_test[[3]][,2]))