### to avoid underflow, use log scale for emission probabilities f1,f2,f3(case 1,2,4)
f1<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  nun<-nu+l1+l2
  kn<-k+l1+l2
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nuvarn<-xsq+ysq+l1*xbar^2+l2*ybar^2+nu*var+k*mu^2-(l1*xbar+l2*ybar+k*mu)^2/(l1+l2+k)
  nom<-gamma(nun/2)*sqrt(k)*(nu*var/2)^(nu/2)
  denom1<-log(2*pi)*(l1+l2)/2
  denom2<-log(gamma(nu/2)*sqrt(kn)*(nuvarn/2)^(nun/2))
  return(log(nom)-denom1-denom2)
}



f2<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nun<-nu+l1+l2
  w1<-(l1*xbar+k*mu)^2/(l1+k)
  w2<-(l2*ybar+k*mu)^2/(l2+k)
  nuvarn<-xsq+ysq+nu*var+l1*xbar^2+l2*ybar^2+2*k*mu^2-w1-w2
  nom<-k*gamma(nun/2)*(nu*var/2)^(nu/2)
  denom2<-log(sqrt(l1+k)*sqrt(l2+k)*gamma(nu/2)*(nuvarn/2)^(nun/2))
  denom1<-log(2*pi)*(l1+l2)/2
  return(log(nom)-denom1-denom2)
}




f3<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nun1<-nu+l1
  nun2<-nu+l2
  nuvarn1<-nu*var+xsq+l1*k*(xbar-mu)^2/(l1+k)
  nuvarn2<-nu*var+ysq+l2*k*(ybar-mu)^2/(l2+k)
  nom<-k*(nu*var/2)^nu*gamma(nun1/2)*gamma(nun2/2)
  denom1<-log(2*pi)*(l1+l2)/2
  denom2<-log(sqrt(l1+k)*sqrt(l2+k)*gamma(nu/2)^2*(nuvarn1/2)^(nun1/2)*(nuvarn2/2)^(nun2/2))
  return(log(nom)-denom1-denom2)
}




#'Emission probabilities of the data
#'
#'Generate the emission densities of 4 hidden states for all test sites in the proposed hierarchical mixture model setting.
#'@param nu The hyperparameter nu.
#'@param var The hyperparameter var.
#'@param k The hyperparameter k.
#'@param mu The hyperparameter mu.
#'@param A The input data.
#'@param n1 The sample size of group 1.
#'@param n2 The sample size of group 2.
#'@return The data frame that has emission densities of 4 hidden states as columns and test sites as rows.
#'@export

emission_probs <- function(nu, var, k, mu, A, n1, n2){

  data_df <- data.frame(A)

  emission_analytical<-apply(data_df,1,function(x){
    a<-f1(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    b<-f2(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    d<-f3(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    return(c(a,b,d))
  })
  emission_analytical <- t(emission_analytical)

  n <- dim(A)[1]
  emission_numeric <- rep(NA, n)
  for (s in 1:n){
    x<-as.numeric(data_df[s,])
    fint<-function(t){
      m1<-prod(dnorm(x[1:n1],t[2]/(1-t[2]^2),sqrt(t[3]/(1-t[3]))))
      m2<-prod(dnorm(x[(n1+1):(n1+n2)],t[2]/(1-t[2]^2),sqrt(t[1]/(1-t[1]))))
      m3<-dnorm(t[2]/(1-t[2]^2),mu,sqrt((t[3]/(1-t[3])+t[1]/(1-t[1]))/(4*k)))
      m4<-(nu*var/2)^nu*exp(-(nu*var)/(2*t[3]/(1-t[3])))*exp(-(nu*var)/(2*t[1]/(1-t[1])))
      m5<-gamma(nu/2)^2*(t[3]/(1-t[3]))^(1+nu/2)*(t[1]/(1-t[1]))^(1+nu/2)
      m6<-(1+t[2]^2)/(1-t[2]^2)^2
      m7<-1/(1-t[3])^2
      m8<-1/(1-t[1])^2
      m1*m2*m3*m4*m6*m7*m8/m5
    }
    emission_numeric[s] <- cubature::adaptIntegrate(fint,lowerLimit=c(0,-1,0),upperLimit=c(1,1,1),maxEval=5000)[1]$integral
    emission_numeric[s] <- as.numeric(emission_numeric[s])
  }

  emission <- cbind(exp(emission_analytical[,1]), exp(emission_analytical[,2]), emission_numeric, exp(emission_analytical[,3]))
  return(emission)
}

