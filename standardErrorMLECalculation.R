############################################################FUNCTIONS############################################################

subWhittackerFunction <- function(t,l,m,p) {
  ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
}


whittackerFunction <- function(l,m,p) {
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * 
    integrate(subWhittackerFunction,l=l,m=m,p=p,lower=0,upper=Inf)$value
}

#Full PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDF <- function(par,y1,y2,output) {
  
  mu1 = par[1]
  mu2 = par[2]
  a = par[3]
  b = par[4]
  
  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  #Returning lo
  if(output=="pdf") {
    C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * tryCatch(
      whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2)), 
      error = function(e) return(NA))
  } else {
    log(C) + log(gamma(b)) + (a+b-1)*log(y1*y2) +(((a-1)/2)-(a+b))*log((y1/mu1)+(y2/mu2)) + (-0.5*((y1/mu1)+(y2/mu2))) + tryCatch(
      log(whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))), 
      error = function(e) return(NA))
  }
  
}

maximiserForPDF <- function(par, y1vector, y2vector) {
  pdfEstimates <- vector()
  for (i in 1:length(y1vector)) {
    tryCatch(
      {      
        y1=y1vector[i]; y2=y2vector[i]; pdfEstimates[i] <- bivGammaPDF(par,y1,y2,output="ll");
      },error=function(e){print("ERROR")})
  }
  sum(-(pdfEstimates[is.na(pdfEstimates)==FALSE]))
}

##########################################################TESTS########################################################################

#Generate random variables from Nadarajah & Gupta's bivariate gamma
n=1000; mu1=1;mu2=2; a=1; b=1;
w<-rbeta(n,a,b)
gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1)
gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2)
par=c(mu1,mu2,a,b)

#Implementating whittaker function


####TEST SUBWHITTAKER FUNCTION
i=2;y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i];
subWhittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),t=10)
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper = Inf)



####TEST WHITTAKER FUNCTION
#i=1;y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i];
#whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))


#Test it works (at least values are between 0 and 1) 
#i=1;y1=gamma_c_mu1[i]; y2=gamma_c_mu2[i];
#par=c(mu1,mu2,a,b)

#######maximiserForPDF output test + use in optim test
bivGammaPDF(c(1,1,1,1),1,1,output="ll")
maximiserForPDF(c(1,1,1,1))
optim(par=c(1.5,1.5,1.5,1.5), fn=maximiserForPDF
      , control = c(maxit=10000,trace=2)
      )


#######Optimisation

#gamma_c_mu1 <- 0
#gamma_c_mu2 <- 0

n=100; mu1=1;mu2=2; a=.75; b=1.25;
par=c(mu1,mu2,a,b)
optim_results_output<-mle_simulation(par,n=1000,sims=100)

mle_simulation <- function(par,n,sims) {
  
  #Make results blank
  optim_results <- data.frame()
  
  for (i in 1:sims) {
    print(i)
    set.seed(i)
    
    mu1 = par[1]
    mu2 = par[2]
    a = par[3]
    b = par[4]
    
    
    #1. Generate
    w<-rbeta(n,a,b)
    gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1)
    gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2)
    
    #2.Calculate full likelihood
    #maximiserForPDF(par)
  
    #3. Get parameter estimates  
    optim_result<-optim(par=c(1.5,1.5,1.5,1.5), fn=maximiserForPDF, y1vector=gamma_c_mu1, y2vector=gamma_c_mu2
          , control = c(maxit=10000,trace=1), lower=c(0.1,0.1,0.1,0.1))$par
    
    print(optim_result)
    
    optim_results[i,1] <- optim_result[1]
    optim_results[i,2] <- optim_result[2]
    optim_results[i,3] <- optim_result[3]
    optim_results[i,4] <- optim_result[4]
  }
  return(optim_results)
}

par(mfrow=c(2,2))
hist(optim_results[,1],main="mu_1")
hist(optim_results[,2],main=paste("mu_2")) ###Estimates /mu
hist(optim_results[,3],main=paste("a"))
hist(optim_results[,4],main=paste("b"))

##################Numerical differentiation - try this next

mu1=1;mu2=2; a=.5; b=1.25;
bivGammaPDF(c(mu1,mu2,a,b),y1,y2,output="pdf")

pdfEx(par,y1,y2)


###INtegration is going to infinity
y1=1; y2=1;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value
y1=.5; y2=.5;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value
y1=.1; y2=.1;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value
y1=.01; y2=.01;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value
y1=.001; y2=.001;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value
y1=.0001; y2=.0001;
integrate(subWhittackerFunction,l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2),lower=0,upper =Inf)$value


#####Do whole whittaker function


y=.5;y2=.5
whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))

y=.2;y2=.2
whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))

y=.1;y2=.1
whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))

y=.01;y2=.01
whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))

y=.0001;y2=.0001
whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))






