####################################################### 1. BASE FUNCTIONS: PDF########################################################

options(scipen=999)

subWhittackerFunction <- function(t,l,m,p) {
  ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
}

whittackerFunction <- function(l,m,p) {
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * 
    integrate(subWhittackerFunction,l=l,m=m,p=p,lower=0,upper=Inf)$value
}

#Full PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDF <- function(par,output) {
  
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

#REPARAMETISEDFull PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDFRE <- function(par,output,parameterisation="B2") {
  
  a = par[3]
  b = par[4]
  y1= par[5]
  y2= par[6]
  mu1=exp(par[1])/a
  
  if (parameterisation=="B2") {
    #B_2 Parameterisation for GJRM
    mu2=exp(par[2])/a
  } else {
    #B_t Parameterisation for GJRM
    mu2=exp(par[2]+par[1])/a 
  }

  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
  #Returning lo
  if(output=="pdf") {
    C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * tryCatch(
      whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2)), 
      error = function(e) return(NA))
  } else {
    #log(C) + log(gamma(b)) + (a+b-1)*log(y1*y2) +(((a-1)/2)-(a+b))*log((y1/mu1)+(y2/mu2)) + (-0.5*((y1/mu1)+(y2/mu2))) + tryCatch(
    #  log(whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))), 
    #  error = function(e) return(NA))
    
    log(C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * tryCatch(
      whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2)), 
      error = function(e) return(NA)))
  }
  
}

bivGammaPDFRE_TEST <- function(par,output) {
  
  a = par[3]
  b = par[4]
  y1= par[5]
  y2= par[6]
  mu1=par[1]
  mu2=par[2]
  #mu1 = a/exp(par[1])
  #mu2 = a/exp(par[2])

  k=a
  theta=mu1
  
  #Returning lo
  if(output=="pdf") {
    (1/(gamma(k)*(theta^k)))*(y1^(k-1))*exp(-y1/theta)
  } else {
    log((1/(gamma(k)*(theta^k)))*(y1^(k-1))*exp(-y1/theta))
  }
}


######TEST
#par=c(10,12,1,1,1/10,1/12)
#bivGammaPDFRE(par,output="ll")
#bivGammaPDFRE_TEST(par,output="ll")

####################################################### 2. TESTS########################################################################

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


####################################################### 3. MLE CALCULATIONS###########################################################

########################################################## 3.1 MLE CALCULATION FUNCTIONS #############################################

#########Calculates log-likelihood for full sample
maximiserForPDF <- function(par, y1vector, y2vector) {
  pdfEstimates <- vector()
  for (i in 1:length(y1vector)) {
    tryCatch(
      {      
        y1=y1vector[i]; y2=y2vector[i]; pdfEstimates[i] <- bivGammaPDFRE(c(par[1:4],y1,y2),output="ll");
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
    
    a = par[3]
    b = par[4]
    #mu1 = par[1]
    #mu2 = par[2]
    mu1 = a/exp(par[1])
    mu2 = a/exp(par[2])
    #mu2 = mu1/exp(par[2])
    
    #1. Generate
    w<-rbeta(n,a,b)
    gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/mu1)
    gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/mu2)
    
    #2.Calculate full likelihood
    #maximiserForPDF(par)
  
    #3. Get parameter estimates  
    optim_result<-optim(par=c(log(1.5),log(1.5),1.5,1.5), fn=maximiserForPDF, y1vector=gamma_c_mu1, y2vector=gamma_c_mu2
          , control = c(maxit=10000,trace=1), lower=c(NA,NA,0.1,0.1))$par
    
    optim_results[i,1] <- optim_result[1]
    optim_results[i,2] <- optim_result[2]
    optim_results[i,3] <- optim_result[3]
    optim_results[i,4] <- optim_result[4]
    optim_results[i,5] <- log(a/exp(par[1]))
    optim_results[i,6] <- log(a/exp(par[2]))
    optim_results[i,7] <- a
    optim_results[i,8] <- b
  }
  return(optim_results)
}

########################################################## 3.2 MLE RUN #############################################

###Quick test of MLE simulation
set.seed(100)
a=1;b=1;
par=c(log(a/1),log(a/2),a,b,1,1)
optim_results_output<-mle_simulation(par[1:4],n=1000,sims=100)


######Finding MLE for parameters across range of values of alpha, beta to estimate MLE SE
optim_results_output_combined<-data.frame()
n=1000;

for (a in c(0.5,1,1.5)) {
  for (b in c(0.5,1,1.5)) {
    mu1=log(a/1);mu2=log(a/2);y1=1;y2=1;
    optim_results_output<-mle_simulation(c(mu1,mu2,a,b),n=n,sims=sims)
    if (nrow(optim_results_output_combined)==0) {
      optim_results_output_combined = optim_results_output
    } else {
      optim_results_output_combined<-rbind(optim_results_output_combined,optim_results_output)  
    }
  }
}

optim_results_output_combined_save<-optim_results_output_combined

#save(optim_results_output_combined_save,file="optim_results_output_combined_save231109.rds")

#########Load and analyse
#load(file="optim_results_output_combined_save231109.rds")
#optim_results_output<-optim_results_output_combined_save
colnames(optim_results_output)<-c("mu1_est","mu2_est","a_est","b_est","mu1_act","mu2_act","a_act","b_act")

aval=1;bval=1;
mu1_hat = optim_results_output[optim_results_output$a_act==aval & optim_results_output$b_act==bval & is.na(optim_results_output$mu1_est)==FALSE,1]
mu2_hat = optim_results_output[optim_results_output$a_act==aval & optim_results_output$b_act==bval & is.na(optim_results_output$mu1_est)==FALSE,2]
a_hat =   optim_results_output[optim_results_output$a_act==aval & optim_results_output$b_act==bval & is.na(optim_results_output$mu1_est)==FALSE,3]
b_hat =   optim_results_output[optim_results_output$a_act==aval & optim_results_output$b_act==bval & is.na(optim_results_output$mu1_est)==FALSE,4]

###Mean and variance

mean(mu1_hat); sd((mu1_hat)); 
mean(mu2_hat); sd((mu2_hat)); 
#mean(a_hat); sd(a_hat)
#mean(b_hat); sd(b_hat)


####Variance / sd for mean at time 1 and 2 are equivalent to alpha * mu1 and alpha * mu2

sd(mu1_hat*a_hat) ###equivalent to manual calculation of Var(XY)= E(X^2 Y^2)- (E(XY))^2
sd(mu2_hat*a_hat) ###SD for mu2 Var(XY)= E(X^2 Y^2)- (E(XY))^2???

sd(mu2_hat-mu1_hat)

################## 4. Numerical differentiation####################

##################### 4.1 FOR TEST FUNCTION######

testFunction <- function(par,output) {
  mu1=par[1];mu2=par[2];a=par[3];b=par[4];y1=par[5];y2=par[6];
  1*mu1^1+2*(mu2^2)*(3*a^3)+4*b^4+5*y1^5+6*y2^6
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

diff2Function<-function(par,hval,hpar) {
  
  d2 <- vector(length = 6)
  
  ###For given parameter
  
  for (i in 1:length(par)) {
    
    hplus1=c(0,0,0,0,0,0); hplus2=c(0,0,0,0,0,0)
    
    if (i == hpar) {
      hplus1[hpar]=hval
      d2[i]=
        (-2*testFunction(par,output="ll") + 
           testFunction(par+hplus1,output="ll") + 
           testFunction(par-hplus1,output="ll"))/(hval^2)    
    }    
    else {
      #### Add an amount to parameter i and add an amount to hpar
      
      #hpar=2;i=3; hplus1=c(0,0,0,0,0,0); hplus2=c(0,0,0,0,0,0)
      
      hplus1[hpar]=hval; hplus2[i]=hval;
      d2[i] <- (   testFunction(par+hplus1+hplus2,output="ll") -
                     testFunction(par+hplus1-hplus2,output="ll") -
                     testFunction(par-hplus1+hplus2,output="ll") +
                     testFunction(par-hplus1-hplus2,output="ll") 
      ) / (4*hval^2)
    }
  }
  return(d2)
}

d2matrix=matrix(nrow=6,ncol=6)
for (i in 1:6) {
  d2matrix[i,]=diff2Function(par,hval=.00001,hpar=i)
}
solve(d2matrix)



##################### 4.2 FOR REAL FUNCTION###########

diffFunction<-function(par,hval,hpar,test=FALSE) {
  
  if (test==TRUE) {pdfFunction<-bivGammaPDFRE_TEST} else {pdfFunction<-bivGammaPDFRE}

  hplus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;
  
  (pdfFunction(par+hplus,output="pdf")-pdfFunction(par-hplus,output="pdf"))/(2*hval)
}

diff2Function<-function(par,hval,hpar,test=FALSE,parameterisation="B2") {
  
  #par;test=TRUE;npar=c(3,1);hval=.00001;hpar=2
  #test=FALSE; npar=c(4,4); hval=.001; par=c(10,12,1,1,1/10,1/12)
  
  if (test==TRUE) {pdfFunction<-bivGammaPDFRE_TEST} else {pdfFunction<-bivGammaPDFRE}
  
  d2 <- vector(length = 6)
  
  ###For given parameter
  
  for (i in 1:length(par)) {
    
    hplus1=c(0,0,0,0,0,0); hplus2=c(0,0,0,0,0,0)
    
    if (i == hpar) {
      hplus1[hpar]=hval
      d2[i]=
          (2*pdfFunction(par,output="ll",parameterisation) -
           pdfFunction(par+hplus1,output="ll",parameterisation) - 
           pdfFunction(par-hplus1,output="ll",parameterisation))/(hval^2)    
    }    
    else {
      #### Add an amount to parameter i and add an amount to hpar
      
      hplus1[hpar]=hval; hplus2[i]=hval;
      d2[i] = (    pdfFunction(par+hplus1-hplus2,parameterisation,output="ll") +
                    pdfFunction(par-hplus1+hplus2,parameterisation,output="ll") -
                    pdfFunction(par-hplus1-hplus2,parameterisation,output="ll") -
                    pdfFunction(par+hplus1+hplus2,parameterisation,output="ll")
      ) / (4*hval*hval)
    }
  }
  return(d2)
}

numericalDerivativeSE_SAMPLE <- function(par,n,test=FALSE,npar=c(1:4),hval=.01,parameterisation="B2") {
  
  w<-rbeta(n=n,par[3],par[4])
  gamma_c_mu1<-w*rgamma(n,shape=par[3]+par[4],scale=exp(par[1])/par[3])
  gamma_c_mu2<-w*rgamma(n,shape=par[3]+par[4],scale=exp(par[2])/par[3])
  
  ses <- replicate(n, diag(6), simplify=F)
  
  for (j in 1:n) {
    par[5:6]<-c(gamma_c_mu1[j],gamma_c_mu1[j])
    d2matrix=matrix(nrow=6,ncol=6)
    for (i in 1:6) {
      d2matrix[i,]=diff2Function(par,hval=hval,hpar=i,test=test,parameterisation)
    }
    #ses[j,1:6]=sqrt(diag(solve(d2matrix)))
    ses[[j]]=d2matrix
  }
  
  avg_ses_n<-apply(simplify2array(ses), 1:2, mean,na.rm=TRUE)
  
  return(sqrt(diag(solve(avg_ses_n[npar,npar]))))
}

numericalDerivativeSE <- function(par,test=FALSE,npar=c(1:4),hval=.01,parameterisation="B2") {
  
  #test=FALSE; npar=c(1:4); hval=.001; par=c(10,12,1,1,1*10,1*12)
  #test=TRUE; npar=c(3,1); hval=.001; par=c(10,12,2,1,2*10,.1*12)
  
  d2matrix=matrix(nrow=6,ncol=6)
  for (i in 1:6) {
    d2matrix[i,]=diff2Function(par,hval=hval,hpar=i,test=test,parameterisation)
  }
  d2matrix[npar,npar]
  ###For test function
  #sqrt(diag(solve(d2matrix[npar,npar]))/1000)
  #eGamma(rgamma(n=1000,shape=par[3],scale=par[1]), method="numerical.MLE")
  return(diag(solve(d2matrix[npar,npar])))
}


#Prove sample is equivalent to choosing MLE value
n=1000
z=matrix(nrow=n/10,ncol=2)
for (i in 1:(n/10)) {
  print(i)
  z[i,1] = numericalDerivativeSE_SAMPLE(par,i*10)[1]
  z[i,2] = numericalDerivativeSE_SAMPLE(par,i*10)[2]
}

par(mfrow=c(1,2))
referenceNumDev=numericalDerivativeSE(par)
plot(1:(n/10),sqrt(z[,1]/n))
abline(h=sqrt(referenceNumDev[1]/n),col="red",xlab="n",ylab="mu1_se",main="Choose Y1 at MLE versus simulation")
plot(1:(n/10),sqrt(z[,2]/n))
abline(h=sqrt(referenceNumDev[2]/n),col="red",xlab="n",ylab="mu2_se",main="Choose Y2 at MLE versus simulation")

#For test function this should be psigamma(a,deriv=1), 1/(mu1^2), and 1/mu1 for other diagnals
# numDerivResults

###WORKING SECTION

numDerivResults <- matrix(nrow=400,ncol=7)
n=1000;
i = 1
origmu1=10; origmu2=12; a=1; b=1
for (a in .1+.1*1:20) {
  for (b in .1+.1*1:20) {
    print(i)
    
    #B_2 parameterisation
    mu1=log(a*1/origmu1); mu2=log(a*1/origmu2);  ###Changes a lot based on y1/y2
    par=c(mu1,mu2,a,b,exp(mu1)/a,exp(mu2)/a)
    numDerivResults[i,1:4] <- sqrt(numericalDerivativeSE(par,parameterisation="B2")/n) #sqrt(numericalDerivativeSE(par))/sqrt(n)
    
    #B_t parameterisation
    mu1=log(a*1/origmu1); mu2=log((1/origmu2)/(1/origmu1));  ###Changes a lot based on y1/y2
    par=c(mu1,mu2,a,b,exp(mu1)/a,exp(mu1+mu2)/a)
    numDerivResults[i,5] <- sqrt(numericalDerivativeSE(par,parameterisation="Bt")/n)[2] #sqrt(numericalDerivativeSE(par))/sqrt(n)
    
    numDerivResults[i,6] <- a
    numDerivResults[i,7] <- b
    i = i+1
  }
}
numDerivResults
#library(ExtDist)
#eGamma(rgamma(n=1000,shape=1,scale=10)) ##SE Benchmark using method of moments

colnames(numDerivResults) <- c("mu1_se","mu2_se_B2","a_se","b_se","mu2_se_Bt","a","b")
#save(numDerivResults,file="numDerivResults.rds")

par(mfrow=c(1,3))
plot(tau,numDerivResults[,1])
plot(tau,numDerivResults[,2])
plot(tau,numDerivResults[,5])


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




