
#Generate random variables from Nadarajah & Gupta's bivariate gamma
n=1000; mu1=1;mu2=2; a=.25; b=1.75;
w<-rbeta(n,a,b)
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu1)
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu2)
par=c(mu1,mu2,a,b)

#Implementating whittaker function

subWhittackerFunction <- function(t,l,m,p) {
  ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
}

####TEST SUBWHITTAKER FUNCTION
i=1;y1=gamma_c_mu1[14]; y2=gamma_c_mu2[14];
subWhittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),t=10)
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper = .Machine$double.xmax, rel.tol=.Machine$double.eps^.05)$value


whittackerFunction <- function(l,m,p) {
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * 
    integrate(subWhittackerFunction,l=l,m=m,p=p,lower=0,upper=Inf)$value
}

####TEST WHITTAKER FUNCTION
#i=1;y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i];
#whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))

#Full PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDF <- function(par,y1,y2) {
  
  mu1 = par[1]
  mu2 = par[2]
  a = par[3]
  b = par[4]
    
  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  #Returning lo
  C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))
  
}


#Test it works (at least values are between 0 and 1) 
#i=1;y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i];
#par=c(mu1,mu2,a,b)

###Working for a single value
optim(par=c(1,1,1,1), fn=bivGammaPDF,y1=gamma_c_mu1[i], y2=gamma_c_mu2[i]
      , control = c(maxit=10000,fnscale=-1),lower=c(0.1,0.1,0.1,0.1),upper=c(10,10,10,10),method = "L-BFGS-B")


maximiserForPDF <- function(par) {
  pdfEstimates <- vector()
  for (i in 1:length(gamma_c_mu1)) {
    tryCatch(
      {      print(i)
      y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i]; pdfEstimates[i] <- bivGammaPDF(par,y1,y2);
      },error=function(e){print("ERROR")})
  }
  sum(-log(pdfEstimates[!is.na(pdfEstimates)]))
}

optim(par=c(1,1,1,1), fn=maximiserForPDF
      , control = c(maxit=10000,trace=1))


################################## FULL FUNCTION

#Getting biased MLEs..?

optim_results <- list()
for (i in 1:100) {
  print(i)
  #1. Generate
  n=100; mu1=1;mu2=2; a=1; b=1;
  w<-rbeta(n,a,b)
  gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu1)
  gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu2)
  par=c(mu1,mu2,a,b)
  
  #2.Calculate full likelihood
  maximiserForPDF(par)

  #3. Get parameter estimates  
  optim_results[i]<-optim(par=c(3,3,3,3), fn=maximiserForPDF
        , control = c(maxit=10000,trace=1))
  
}

as.data.frame(optim_results)[1:100]



##################Numerical differentiation - try this next




