
#Generate random variables from Nadarajah & Gupta's bivariate gamma
n=100; a=.25; b=1.75; mu1=1;mu2=2;
w<-rbeta(n,a,b)
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu1)
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu2)

#Implementating whittaker function
whittackerFunction <- function(l,m,p) {
  subWhittackerFunction <- function(t) {
    ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
  }
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * integrate(subWhittackerFunction,0,Inf)$value
}

#Full PDF with 4 parameters for the 2 covariates
bivGammaPDF <- function(y1,y2,m1,m2,a,b) {
    
  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/m1)+(y2/m2))
  
}


#Test it works (at least values are between 0 and 1) 
bivGammaPDF(1,1,1,1,1,1)

#Alternative parameterisation
#m1=1;m2=2;a=.25;b=1.25;
#C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) );
#bivGammaPDFy1y2 <- function (y1,y2) {
#  C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/m1)+(y2/m2)) 
#}
#bivGammaPDFy1y2(1,2)

for i in 1:length(gamma_c_mu1)



f <- function (x, a) (x - a)^2
xmin <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)
xmin
