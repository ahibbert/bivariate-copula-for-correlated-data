    source("common_functions.R")


    true_out=calcTrueCovariateValues(
  n=1000,
  a=1,
  b=1,
  c=1,
  mu1=2,
  mu2=5,
  dist="NO",
  x1=1,
  x2=.01
    )
    true_out$par

ses=get_true_ses(
  n=100,
  a=1,
  b=1,
  c=1,
  mu1=1,
  mu2=2,
  dist="NO",
  x1=1,
  x2=.01
)

ses1000=get_true_ses(
  n=1000,
  a=1,
  b=1,
  c=1,
  mu1=1,
  mu2=2,
  dist="GA",
  x1=1,
  x2=.01
)
ses1000