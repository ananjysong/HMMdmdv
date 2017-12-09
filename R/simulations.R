#'Simulation parameter generation
#'
#'Generate true parameters for simulations
#'@param n_state Default: 4. Number of hidden states
#'@return The list of true parameters for simulations.
#'@export

truePara_generation <- function(n_state){
  initProb = runif(n_state)
  initProb = initProb/sum(initProb)
  transProb = matrix(NA, n_state, n_state)
  for (i in 1:n_state){
    transProb[i,] = runif(n_state)
    transProb[i,] = transProb[i,]/sum(transProb[i,])
  }
  nu0 = runif(1, 4, 6)
  var0 = runif(1, 18, 26)
  k0 = sample(2:4, 1)
  mu0 = runif(1, 0.5, 2)
  return(list(nu0, var0, k0, mu0, initProb, transProb))
}


#'Simulation data generation
#'
#'Generate data for simulations
#'@param nu0 Hyperparameter nu0.
#'@param var0 Hyperparameter var0.
#'@param k0 Hyperparameter k0.
#'@param mu0 Hyperparameter mu0.
#'@param n1 Hyperparameter n1.
#'@param n2 Hyperparameter n2.
#'@param initProb Initial state probabilities.
#'@param transProb Transition probability matrix.
#'@return The generated data and the true hidden states.
#'@export
#'
data_generation <- function(nu0, var0, k0, mu0, n, n1, n2, initProb, transProb){
  X<-matrix(NA, n, n1)
  Y<-matrix(NA, n, n2)
  Z<-rep(NA,n)
  Z[1]<-sample(1:4, 1, prob = initProb)

  for (i in 2:n){
    Z[i]<-sample(1:4, 1, prob = transProb[Z[i-1],])
  }


  for (i in 1:n){
    if(Z[i]==1)
    {
      sigma<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      mu<-rnorm(1,mu0,sqrt(sigma/k0))
      X[i,]<-rnorm(n1,mu,sqrt(sigma))
      Y[i,]<-rnorm(n2,mu,sqrt(sigma))
    }
    else if(Z[i]==2){
      sigma<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      mu1<-rnorm(1,mu0,sqrt(sigma/k0))
      mu2<-rnorm(1,mu0,sqrt(sigma/k0))
      X[i,]<-rnorm(n1,mu1,sqrt(sigma))
      Y[i,]<-rnorm(n2,mu2,sqrt(sigma))
    }
    else if(Z[i]==3){
      sigma1<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      sigma2<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      mu1<-rnorm(1,mu0,sqrt(sigma1/k0))
      mu2<-rnorm(1,mu0,sqrt(sigma2/k0))
      mu<-(mu1+mu2)/2
      X[i,]<-rnorm(n1,mu,sqrt(sigma1))
      Y[i,]<-rnorm(n2,mu,sqrt(sigma2))
    }
    else{
      sigma1<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      sigma2<-pscl::rigamma(1,nu0/2,nu0*var0/2)
      mu1<-rnorm(1,mu0,sqrt(sigma1/k0))
      mu2<-rnorm(1,mu0,sqrt(sigma2/k0))
      X[i,]<-rnorm(n1,mu1,sqrt(sigma1))
      Y[i,]<-rnorm(n2,mu2,sqrt(sigma2))
    }
  }

  A<-cbind(X,Y)
  print('data generated!')
  return(list(A, Z))
}
