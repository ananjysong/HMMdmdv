# HMMdmdv
R package: Hidden Markov Model for Bayesian Inference of large-scale mean and variance testing

The HMMdmdv package implements Hidden Markov Model and makes Bayesian inference for large-scale testing of mean and variance in a proposed Bayesian hierarchical mixture model framework.




## Introduction

This is an example of using HMMdmdv package in R. HMMdmdv implements Hidden Markov Model for large-scale multiple testing of differential mean and variance in the framework of a proposed hierarchical mixture model. This vignette aims to demonstrate the usage of HMMdmdv through some example codes and example data. 


## Usage of HMMdmdv

### (1) installation
```
install.packages("HMMdmdv_0.1.1.tar.gz", repos = NULL, type="source")
```


### (2) loading example data
The example dataset has 1000 rows for 1000 test sites and 100 columns for 100 samples. In each individual test, we are comparing the first 50 samples with the second 50 samples. 
```{r}
library(HMMdmdv)
# Loading dataset
data(example_data)
# Check dataset
example_data[1:10, 1:8]
dim(example_data)
# 1000 tests, 50 samples in each group
n<-1000
n1 <- 50
n2 <- 50
```







### (3) Initial (heuristic) parameter estimate

```{r}
# remove case 3
data_df <- remove_case3(dat = example_data, mean_thresholdPV = 0.1, 
var_thresholdPV = 0.05, n1 = n1, n2 = n2)
# parameter estimate
init_para_est <- init_est(dat_df = data_df, mean_thresholdPV = 0.05, 
var_thresholdPV = 0.1, n1 = n1, n2 = n2, n = n)

```     


### (4) Independent Gaussian Mixture parameter estimate
```{r message=FALSE, warning=FALSE}
# upper limit of EM iterations
niter<-1000
# parameter estimate
indep_para_est <- runEM(dat_df = data_df, n1 = n1, n2 = n2, 
init_para_est= init_para_est, niter = niter)
```

### (5) Calculate emission densities 
```{r}
emissions <- emission_probs(indep_para_est[5], indep_para_est[6], 
indep_para_est[7], indep_para_est[8], example_data, n1 = n1, n2 = n2)
head(emissions, 10)
```

   


### (6) HMM parameter estimates 
```{r}
HMM_resList <- runHMM(emissions)

```

### (7) Posterior inference of hidden states 
```{r}
result <- posterior_inference(HMM_resList[[1]], HMM_resList[[2]], train = FALSE)

head(result, 10)
```



## More Information
Details for arguments and functions can be found by typing e.g. `help(package="HMMdmdv")`, `?runEM`.  

