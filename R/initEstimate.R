#'Subset out the case 3 data
#'
#'Subset out the case 3 data with same mean and different variance. The case 3 data
#'@param dat The input data (matrix) with rows as different test sites, columns as the samples from two groups with sample size as n1, n2.
#'@param mean_thresholdPV The p-value threshold for mean test to subset out case 3 data.
#'@param var_thresholdPV The p-value threshold for variance test to subset out case 3 data.
#'@param n1 The sample size of group 1.
#'@param n2 The sample size of group 2.
#'@return The data frame that has case 3 removed.
#'@export

remove_case3 <- function(dat, mean_thresholdPV, var_thresholdPV, n1, n2){
  # threshold to subset case3 data
  data_df <- data.frame(dat)
  pv3 <- apply(dat, 1, function(x) t.test(x[1:n1], x[(n1+1):(n1+n2)])$p.value)
  qv3 <- apply(dat, 1, function(x){
    l <- list(x[1:n1], x[(n1+1):(n1+n2)])
    bartlett.test(l)$p.value
  })

  data_df$pv3 <- pv3
  data_df$qv3 <- qv3
  case3 <- subset(data_df, pv3 >= mean_thresholdPV & qv3 <= var_thresholdPV)
  case3 <- case3[, -c(n1+n2+1, n1+n2+2)]

  rowname <- rownames(case3)
  rowno <- as.numeric(rowname)
  data_df <- data_df[-rowno, 1:(n1+n2)]
  return(data_df)
}



#'Initial parameter estimate
#'
#'Subset out the case 3 data with same mean and different variance. The case 3 data
#'@param dat_df The input data frame with case 3 removed.
#'@param mean_thresholdPV The p-value threshold for mean test to subset case 1 data.
#'@param var_thresholdPV The p-value threshold for variance test to subset case 1 data.
#'@param n1 The sample size of group 1.
#'@param n2 The sample size of group 2.
#'@param n The total number of tests (test sites, methylation probes).
#'@return Initial parameter estimate of state proportion parameters and hyper-parameters.
#'@export

init_est <- function(dat_df, mean_thresholdPV, var_thresholdPV, n1, n2, n){
  # input dat_df (case 3 removed)
  # threshold to subset case1 data
  pv <- apply(dat_df, 1, function(x) t.test(x[1:n1], x[(n1+1):(n1+n2)])$p.value)
  qv <- apply(dat_df, 1, function(x){
    l <- list(x[1:n1], x[(n1+1):(n1+n2)])
    bartlett.test(l)$p.value
  })
  dat_df$pv <- pv
  dat_df$qv <- qv
  case1 <- subset(dat_df, pv >= mean_thresholdPV & qv >= var_thresholdPV)
  case1 <- case1[, -c(n1+n2+1, n1+n2+2)]

  length1 <- dim(case1)[1]

  # use the average values in case 1 dataset for initial hyperparameter estimate
  mu_1 <- apply(case1, 1, mean)
  mu_1 <- as.numeric(mu_1)
  mu_0 <- mean(mu_1)

  mu_1_var <- var(mu_1)
  sigma_1 <- apply(case1, 1, var)
  sigma_1 <- as.numeric(sigma_1)
  sigma_1 <- sigma_1*(n1+n2-1)/(n1+n2-2)
  sigma_1_mean <- mean(sigma_1)
  sigma_1_var <- var(sigma_1)

  k_0 <- sigma_1_mean/mu_1_var
  nu_0 <- 4 + 2*sigma_1_mean^2/sigma_1_var
  var_0 <- (sigma_1_mean^3 + sigma_1_mean*sigma_1_var)/(sigma_1_mean^2 + 2*sigma_1_var)

  case2 <- subset(dat_df, pv <= mean_thresholdPV & qv >= var_thresholdPV)
  case2 <- case2[, -c(n1+n2+1, n1+n2+2)]
  length2 <- dim(case2)[1]

  p3_0 <- (n - dim(dat_df)[1])/n
  p1_0 <- length1/dim(dat_df)[1]
  p2_0 <- length2/dim(dat_df)[1]
  p4_0 <- 1 - p1_0 - p2_0
  p1_0 <- (1 - p3_0) * p1_0
  p2_0 <- (1 - p3_0) * p2_0
  p4_0 <- (1 - p3_0) * p4_0
  print('initial parameter estimate:')
  print(c(p1_0, p2_0, p3_0, p4_0, nu_0, var_0, k_0, mu_0))
  return(c(p1_0, p2_0, p3_0, p4_0, nu_0, var_0, k_0, mu_0))
}
