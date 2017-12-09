#'EM on the data without case 3.
#'
#'Simply take the data without case 3 as independent Gaussian mixture and estimate the hyperparameters.
#'@param dat_df The input data frame.
#'@param EM_threshold The threshold of parameter diffrences among iterations. Default: 0.0001. Control the number of EM iterations.
#'@param n1 The sample size of group 1.
#'@param n2 The sample size of group 2.
#'@param init_para_est The initial parameter estimates.
#'@param niter The upper limit for number of EM iterations.
#'@return The number of iterations and the EM estimate of hyperparameters together with proportion parameters in the setting of no case 3.
#'@export

runEM_1 <- function(dat_df, EM_threshold = 0.0001, n1, n2, init_para_est, niter){

  init_para_est_1 <- rep(NA, 7)
  init_para_est_1[4:7] <- init_para_est[5:8]
  total <- init_para_est[1] + init_para_est[2] + init_para_est[4]
  prop1 <- init_para_est[1]/total
  prop2 <- init_para_est[2]/total
  prop4 <- init_para_est[4]/total
  init_para_est_1[1] <- prop1
  init_para_est_1[2] <- prop2
  init_para_est_1[3] <- prop4

  theta <- matrix(NA, nrow = niter, ncol = 7)
  theta[1,]<-c(0.5,0.25,0.25,6,4,2,0.2)
  theta[2,] <- init_para_est_1

  i <- 2
  while(max(abs(theta[i,]-theta[i-1,]))>0.0001){
    p1<-theta[i,1]
    p2<-theta[i,2]
    p3<-theta[i,3]
    nu<-theta[i,4]
    var<-theta[i,5]
    k<-theta[i,6]
    mu<-theta[i,7]

    epc<-apply(dat_df,1,function(x){
      d1<-f1(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      d2<-f2(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      d3<-f3(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      common_term<-log(p1)+d1+log(1+exp(log(p2)-log(p1)+d2-d1)+exp(log(p3)-log(p1)+d3-d1))
      fai1<-log(p1)+d1-common_term
      fai2<-log(p2)+d2-common_term
      fai3<-log(p3)+d3-common_term
      fai1<-exp(fai1)
      fai2<-exp(fai2)
      fai3<-exp(fai3)
      return(c(fai1,fai2,fai3))
    })
    epct<-t(epc)

    keep <- complete.cases(epct)
    epct <- epct[keep,]
    dat_df <- dat_df[keep,]
    theta[i+1,1]<-sum(epct[,1])/dim(dat_df)[1]
    theta[i+1,2]<-sum(epct[,2])/dim(dat_df)[1]
    theta[i+1,3]<-sum(epct[,3])/dim(dat_df)[1]


    loglik<-function(para){
      B<-apply(dat_df,1,function(x){
        a1<-f1(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        b1<-f2(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        c1<-f3(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        return(c(a1,b1,c1))
      })
      return(-sum(diag(epct%*%B)))}

    theta[i+1,4:7]<-optim(c(nu,var,k,mu),loglik)$par
    i <- i + 1
  }
  return(c(i,theta[i,]))
}



#'EM on the data without case 3.
#'
#'Simply take the data without case 3 as independent Gaussian mixture and estimate the hyperparameters.
#'@param dat_df The input data frame.
#'@param EM_threshold The threshold of parameter diffrences among iterations. Default: 0.0001. Control the number of EM iterations.
#'@param n1 The sample size of group 1.
#'@param n2 The sample size of group 2.
#'@param init_para_est The initial parameter estimates.
#'@param niter The upper limit for number of EM iterations.
#'@return The EM estimate of hyperparameters together with proportion parameters in the setting of no case 3.
#'@export

runEM <- function(dat_df, EM_threshold = 0.0001, n1, n2, init_para_est, niter){

  EM_para_est <- runEM_1(dat_df = dat_df, EM_threshold = EM_threshold, n1 = n1, n2 = n2, init_para_est = init_para_est, niter = niter)
  i <- EM_para_est[1]
  theta_est <- EM_para_est[-1]

  PI <- (1-init_para_est[3])*theta_est[1:3]
  PI<-c(PI[1],PI[2],init_para_est[3],PI[3])
  nu<-theta_est[4]
  var<-theta_est[5]
  k<-theta_est[6]
  mu<-theta_est[7]

  indep_para_est <- c(PI, nu, var, k, mu)
  print('independent EM parameter estimate:')
  print(indep_para_est)
  return(indep_para_est)
}
