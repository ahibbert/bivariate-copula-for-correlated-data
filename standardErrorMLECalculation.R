
#Generate random variables from Nadarajah & Gupta's bivariate gamma
n=100; a=.25; b=1.75; mu1=1;mu2=2;
w<-rbeta(n,a,b)
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu1)
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu2)

#Implementating whittaker function
whittackerFunction <- function(l,m,p) {
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * integrate(subWhittackerFunction,lower=0,upper=Inf)$value
}

subWhittackerFunction <- function(t) {
  ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
}

#Full PDF with 4 parameters for the 2 covariates
bivGammaPDF <- function(y1,y2,mu1,mu2,a,b) {
    
  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))
  
}


#Full PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDF <- function(par,y1,y2) {
  
  mu1 = par[1]
  mu2 = par[2]
  a = par[3]
  b = par[4]
    
  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  #Returning lo
  -log(C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2)))
  
}

par=c(mu1,mu2,a,b)
#Test it works (at least values are between 0 and 1) 
bivGammaPDF(par,1,1)

#maximiserForPDF <- function(mu1,mu2,a,b) {
#  pdfEstimates <- vector()
#  for (i in 1:length(gamma_c_mu1)) {
#    y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i]; pdfEstimates[i] <- bivGammaPDF(y1,y2,mu1,mu2,a,b);
#  }
#  sum(-2*log(pdfEstimates))
#}


optim(par=c(1,1,1,1), fn=bivGammaPDF,y1=gamma_c_mu1[1], y2=gamma_c_mu2[1], control = c(maxit=10000,fnscale=-1))

