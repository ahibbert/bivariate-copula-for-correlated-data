###############################################BASE FUNCTIONS: PDF########################################################

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

#REPARAMETISEDFull PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDFRE <- function(par,output) {
  
  mu1 = par[1]
  mu2 = par[2]
  a = par[3]
  b = par[4]
  y1= par[5]
  y2= par[6]
  
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


##############################################################################MLE CALCULATIONS###########################################################

#########Calculates log-likelihood for full sample
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

#############Simulates bivariate gamma and estimates best parameters based on log likelihood
mle_simulation <- function(par,n,sims) {
  
  #Make results blank
  optim_results <- data.frame()
  
  for (i in 1:sims) {
    print(i)
    set.seed(i)
    
    mu1 = 1/par[1]
    mu2 = 1/par[2]
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
    optim_results[i,5] <- mu1
    optim_results[i,6] <- mu2
    optim_results[i,7] <- a
    optim_results[i,8] <- b
  }
  return(optim_results)
}

######Finding MLE for parameters across range of values of alpha, beta to estimate MLE SE
optim_results_output_combined<-data.frame()
n=100; mu1=1;mu2=2; a=NA; b=NA;

for (a in c(0.5,1,1.5)) {
  for (b in c(0.5,1,1.5)) {
    optim_results_output<-mle_simulation(c(mu1,mu2,a,b),n=n,sims=2)
    if (nrow(optim_results_output_combined)==0) {
      optim_results_output_combined = optim_results_output
    } else {
      optim_results_output_combined<-rbind(optim_results_output_combined,optim_results_output)  
    }
  }
}

load(optim_results_output_combined)

colnames(optim_results_output_combined)<-c("mu1_est","mu2_est","a_est","b_est","mu1_act","mu2_act","a_act","b_act")
optim_results_output_combined

#optim_results_output_075_125<-optim_results_output

par(mfrow=c(2,2))
hist(optim_results_output[c(1:7,9:100),1],main="mu_1")
hist(optim_results_output[c(1:7,9:100),2],main=paste("mu_2")) ###Estimates /mu
hist(optim_results_output[c(1:7,9:100),3],main=paste("a"))
hist(optim_results_output[c(1:7,9:100),4],main=paste("b"))

###Mean and variance
mean(optim_results_output[c(1:7,9:100),1]); sd(optim_results_output[c(1:7,9:100),1]) 
mean(optim_results_output[c(1:7,9:100),2]); sd(optim_results_output[c(1:7,9:100),2]) 
mean(optim_results_output[c(1:7,9:100),3]); sd(optim_results_output[c(1:7,9:100),3]) 
mean(optim_results_output[c(1:7,9:100),4]); sd(optim_results_output[c(1:7,9:100),4]) 

####
sd(optim_results_output[c(1:7,9:100),1]*optim_results_output[c(1:7,9:100),3])
sd(optim_results_output[c(1:7,9:100),2]*optim_results_output[c(1:7,9:100),3])

##################Numerical differentiation - try this next####################

mu1=1;mu2=2; a=.5; b=1.25; y1=1;y2=1;
par=c(mu1,mu2,a,b,y1,y2); h=c(0,0,0,0,0,0)

##############FOR TEST FUNCTION######

testFunction <- function(par) {
  mu1=par[1];mu2=par[2];a=par[3];b=par[4];y1=par[5];y2=par[6];
  1*mu1^1+2*mu2^2+3*a^3+4*b^4+5*y1^5+6*y2^6
}

hval=.00001; hpar=1
diffFunction<-function(par,hval,hpar) {
  hplus=c(0,0,0,0,0,0); hminus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;hminus[hpar]=-hval
  (testFunction(par+hplus) - testFunction(par+hminus))/(2*hval)
}

for (i in 1:6) {
  print(diffFunction(par,hval,hpar=i))
}

diff2Function<-function(par,hval,hpar) {
  hplus=c(0,0,0,0,0,0); hminus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;hminus[hpar]=-hval
  
  (testFunction(par+hplus)+testFunction(par-hplus)-2*testFunction(par))/(hval^2)
}

for (i in 1:6) {
  print(diff2Function(par,hval,hpar=i))
}


########FOR REAL FUNCTION###########

diffFunction<-function(par,hval,hpar) {
  hplus=c(0,0,0,0,0,0); hminus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;hminus[hpar]=-hval
  
  (bivGammaPDFRE(par+hplus,output="pdf")-bivGammaPDFRE(par-hplus,output="pdf"))/(2*hval)
}

diff2Function<-function(par,hval,hpar) {
  hplus=c(0,0,0,0,0,0); hminus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;hminus[hpar]=-hval
  
  (bivGammaPDFRE(par+hplus,output="pdf")+bivGammaPDFRE(par-hplus,output="pdf")-2*bivGammaPDFRE(par,output="pdf"))/(hval^2)
}

derivatives<-data.frame(matrix(0,nrow=1,ncol=18))

j=1
for (a in c(0.5,1,1.5)) {
  for (b in c(0.5,1,1.5)) {
    for (y1 in c(0.5,1,1.5)) {
      for (y2 in c(0.5,1,1.5)) {
        par = c(1,2,a=a,b=b,y1=y1,y2=y2)  
        derivatives[j,13]=par[1]
        derivatives[j,14]=par[2]
        derivatives[j,15]=par[3]
        derivatives[j,16]=par[4]
        derivatives[j,17]=par[5]
        derivatives[j,18]=par[6]
        for (i in 1:6) {
          derivatives[j,i]=(diffFunction(par,hval,hpar=i))
          derivatives[j,i+6]=(diff2Function(par,hval,hpar=i))
        }
        j=j+1
      }
    }
  }
}
derivatives

summary(derivatives)

par(mfrow=c(2,3))
hist(derivatives[,1])
hist(derivatives[,2])
hist(derivatives[,3])
hist(derivatives[,4])
hist(derivatives[,5])
hist(derivatives[,6])

par(mfrow=c(2,3))
hist(derivatives[,7])
hist(derivatives[,8])
hist(derivatives[,9])
hist(derivatives[,10])
hist(derivatives[,11])
hist(derivatives[,12])



###Need to calculate these derivatives at a given point

##########################################Reference: Whittaker function bounds#####################################
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






