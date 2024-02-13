generateBivDist <- function(n,a,b,c,mu1,mu2,dist) {
  
  set.seed(100)
  
  if(dist=="GA") {
    #Simulating bivariate random variable according to functional input
    w<-rbeta(n,a,b)
    margin_1<-w*rgamma(n,shape=a+b,scale=mu1)
    margin_2<-w*rgamma(n,shape=a+b,scale=mu2)
    
  }
  if(dist=="NO") {
    require(MASS)
    normData<-mvrnorm(n,mu=c(mu1,mu2),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    margin_1<-normData[,1]
    margin_2<-normData[,2]
  }
  
  #Transforming data to format required for random effect models
  patient<-as.factor(seq(1:n))
  dataset<-as.data.frame(rbind(cbind(patient,margin_1,0)
                               ,cbind(patient,margin_2,1)))
  colnames(dataset)<-c("patient","random_variable","time")
  
  dataset<-dataset[order(dataset$patient),]
  
  return(dataset)
}
