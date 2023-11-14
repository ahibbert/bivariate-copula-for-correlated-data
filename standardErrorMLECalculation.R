####################################################### 1. BASE FUNCTIONS: PDF########################################################

subWhittackerFunction <- function(t,l,m,p) {
  ((t^(m-l-0.5))*((1+t)^(m+l-0.5)))*exp(-p*t)
}


whittackerFunction <- function(l,m,p) {
  (((p^(m+.5))*exp(-p/2))/(gamma(m-l+0.5))) * 
    integrate(subWhittackerFunction,l=l,m=m,p=p,lower=0,upper=Inf)$value
}

#Full PDF with 4 parameters for the 2 covariates - returns PDF
#bivGammaPDF <- function(par,output) {
#  
#  mu1 = par[1]
#  mu2 = par[2]
#  a = par[3]
#  b = par[4]
#  y1= par[5]
#  y2= par[6]
#  
#  C = 1 / ( ((mu1*mu2)^(a+b)) * gamma(a+b) * gamma(a) * gamma(b) )
#  #Returning lo
#  if(output=="pdf") {
#    C * gamma(b) * ((y1*y2)^(a+b-1)) * (((y1/mu1)+(y2/mu2))^(((a-1)/2)-(a+b))) * exp(-0.5*((y1/mu1)+(y2/mu2))) * tryCatch(
#      whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2)), 
#      error = function(e) return(NA))
#  } else {
#    log(C) + log(gamma(b)) + (a+b-1)*log(y1*y2) +(((a-1)/2)-(a+b))*log((y1/mu1)+(y2/mu2)) + (-0.5*((y1/mu1)+(y2/mu2))) + tryCatch(
#      log(whittackerFunction(l=a+((1-a)/2),m=a+b-(a/2),p=(y1/mu1)+(y2/mu2))), 
#      error = function(e) return(NA))
#  }
#  
#}

#REPARAMETISEDFull PDF with 4 parameters for the 2 covariates - returns PDF
bivGammaPDFRE <- function(par,output) {
  
  a = par[3]
  b = par[4]
  y1= par[5]
  y2= par[6]
  #mu1=par[1]
  #mu2=par[2]
  mu1 = a/exp(par[1])
  mu2 = a/exp(par[2])
  
  #mu2=a/(exp(par[1])*exp(par[2]))
  #mu1=1/par[1]
  #mu2=1/par[2]
  
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
    mu1 = par[1]
    mu2 = par[2]
    #mu1 = a/exp(par[1])
    #mu2 = a/exp(par[2])
    #mu2 = mu1/exp(par[2])
    
    #1. Generate
    w<-rbeta(n,a,b)
    gamma_c_mu1<-w*rgamma(n,shape=a+b,scale=1/(a/exp(par[1])))
    gamma_c_mu2<-w*rgamma(n,shape=a+b,scale=1/(a/exp(par[2])))
    
    #2.Calculate full likelihood
    #maximiserForPDF(par)
  
    #3. Get parameter estimates  
    optim_result<-tryCatch(
      {      optim(par=c(log(a/((a/exp(par[1])))),log(a/((a/exp(par[2])))),a+.1,b+.1), fn=maximiserForPDF, y1vector=gamma_c_mu1, y2vector=gamma_c_mu2
          , control = c(maxit=10000,trace=1), lower=c(-Inf,-Inf,0.05,0.05))$par },error=function(e){NA})
    
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
set.seed(99)
a=1;b=1;
par=c(log(a/10),log(a/12),a,b,a/10,a/12)
optim_results_output<-mle_simulation(par[1:4],n=1000,sims=100)

sd(optim_results_output[!is.na(optim_results_output[,1]),1])
sd(optim_results_output[!is.na(optim_results_output[,2]),2])

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

################## 4. Numerical differentiation - try this next####################

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

diffFunction<-function(par,hval,hpar) {
  hplus=c(0,0,0,0,0,0);
  hplus[hpar]=hval;
  
  (bivGammaPDFRE(par+hplus,output="pdf")-bivGammaPDFRE(par-hplus,output="pdf"))/(2*hval)
}

diff2Function<-function(par,hval,hpar) {
  
  d2 <- vector(length = 6)
  
  ###For given parameter
  
  for (i in 1:length(par)) {
    
    hplus1=c(0,0,0,0,0,0); hplus2=c(0,0,0,0,0,0)
    
    if (i == hpar) {
      hplus1[hpar]=hval
      d2[i]=
        (-2*bivGammaPDFRE(par,output="ll") +
           bivGammaPDFRE(par+hplus1,output="ll") + 
           bivGammaPDFRE(par-hplus1,output="ll"))/(hval^2)    
    }    
    else {
      #### Add an amount to parameter i and add an amount to hpar
      
      hplus1[hpar]=hval; hplus2[i]=hval;
      d2[i] <- (   bivGammaPDFRE(par+hplus1+hplus2,output="ll") -
                     bivGammaPDFRE(par+hplus1-hplus2,output="ll") -
                     bivGammaPDFRE(par-hplus1+hplus2,output="ll") +
                     bivGammaPDFRE(par-hplus1-hplus2,output="ll") 
      ) / (4*hval^2)
    }
  }
  d2[,2:3]
  return(d2)
}

numericalDerivativeSE_forMLE <- function(par) {
  
  ses=matrix(ncol=8,nrow=1)
  #par=c(mu1,mu2,a,b,y1,y2)
  #par=c(10,12,1,1,1/10,1/12)
  d2matrix=matrix(nrow=6,ncol=6)
  j=1
  for (i in 1:6) {
    d2matrix[i,]=diff2Function(par,hval=.0001,hpar=i)
  }
  
  ses[j,1:6]=diag(solve(d2matrix))
  ses[j,7]=d2matrix[3,1]
  ses[j,8]=d2matrix[1,3]
  return(ses)
}

###WORKING SECTION
numDerivResults <- matrix(nrow=400,ncol=10)
i = 1;
for (a in .1+.1*1:20) {
  for (b in .1+.1*1:20) {
    
    print(i)
    
    ###Change to estimate of amu??
    #mu1=log(a/10);mu2=log(10/12);  ###Changes a lot based on y1/y2
    #par=c(mu1,mu2,a,b,a/10,a/12)
    
    mu1=log(a/10);mu2=log(a/12);  ###Changes a lot based on y1/y2
    par=c(mu1,mu2,a,b,a/10,a/12)
    
    #mu1=10;mu2=12;
    #par=c(mu1,mu2,a,b,a/mu1,a/mu2)
    
    #mu1=1/10;mu2=1/12;
    #par=c(mu1,mu2,a,b,a*mu1,a*mu2)
    
    numDerivResults[i,1:8] <- numericalDerivativeSE_forMLE(par)
    
    numDerivResults[i,9] <- a
    numDerivResults[i,10] <- b
    i = i+1
  }
}
numDerivResults
sqrt(numDerivResults)
nrow(numDerivResults[is.na(numDerivResults[,1])|is.na(numDerivResults[,2]),])


#(numDerivResults[,1]/(10^2)) #Estimate for variance of log(1/X)
#(numDerivResults[,1]/(10^2))/(1000) #Variance of log(1/X) MLE for 1000 obs
#sqrt((numDerivResults[,1]/(10^2))/1000) #Variance of log(1/X)
#tauPlusNumDivResults<-as.data.frame(cbind(tau,sqrt((numDerivResults[,1]/(10^2))/1000),sqrt((numDerivResults[,2]/(12^2))/1000)))
tauPlusNumDivResults<-as.data.frame(cbind(tau,sqrt(numDerivResults[,1]/1000),sqrt(numDerivResults[,2]/1000)))
colnames(tauPlusNumDivResults) <- c("tau","mu1_se","mu2_se")

par(mfrow=c(1,2))
plot(tau,sqrt(numDerivResults[,1]/1000))
plot(tau,sqrt(numDerivResults[,2]/1000))

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






