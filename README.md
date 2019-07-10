# Gaussian Mirror for High-dimensional Controled Variable Selection


## Installation
```R
install_github("DeveloperName/PackageName")
```

## Examples
We use the following examples to illustrate the basic usage of the GM package. The GM method is designed for any design matrix. We give several examples based on design matrix with different correlation structures. For simplicity, we will use synthetic data constructed from a linear model such that the response only depends on a small fraction of the variables.

## Example 1
In this example, we generate the design matrix with i.i.d. rows, and eachrow is generated from AR(1) model with autocorrelation 0.5. 


```R
n <- 300  # Number of observations
p <- 1000  # Number of predictors included in model
q <- 60  # Number of true predictors
amplitude = 20 # signal amplitude (for noise level = 1) 


# Generate the variables from a multivariate normal distribution
times <- 1:p
sigma <- 1
rho = 0.5
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
x = mvrnorm(n, mu = rep(0,p), Sigma=V)

# Generate the response.
beta_q = rnorm(q, 0, amplitude)/sqrt(n)
true_index= sample(p,q)
beta = numeric(p)
beta[true_index]=beta_q
mu = x %*%  beta
y = mu + rnorm(n)

# Fit the Gaussian mirror model
fit = gm(y, x, ncores=4)
# fit = gm(y, x, ncores=4,var.est=T) # Run GM with estiamtion of the variance of FDR

```



## Example 2
In this example, we generate the design matrix with i.i.d. rows, and eachrow is generated from multivariate Gaussian distribution with mean zero and covariance matrix with constant correlation coefficient. 

```R
n <- 300  # Number of observations
p <- 1000  # Number of predictors included in model
q <- 60  # Number of true predictors
amplitude = 20 # signal amplitude (for noise level = 1) 


# Generate the variables from a multivariate normal distribution
times <- 1:p
sigma <- 1
rho = 0.5
V= matrix(rho, p,p)
diag(V)= 1
x = mvrnorm(n, mu = rep(0,p), Sigma=V)

# Generate the response.
beta_q = rnorm(q, 0, amplitude)/sqrt(n)
true_index= sample(p,q)
beta = numeric(p)
beta[true_index]=beta_q
mu = x %*%  beta
y = mu + rnorm(n)

# Fit the Gaussian mirror model
fit = gm(y, x, ncores=4)
# fit = gm(y, x, ncores=4,var.est=T) # Run GM with estiamtion of the variance of FDR

```


## Example 3
In this example, we consider non-Gaussian design matrix. We use SNP data set from a recent genome wide association study (GWAS) including a panel of 292 tomato accessions. We randomly select 1000 SNP as the design matrix.

```R
n = nrow(tomato)
p <- 1000  # Number of predictors included in model
q <- 60  # Number of true predictors
x = as.matrix(tomato[,sample(1:ncol(tomato),p)])
amplitude=40
beta_q = rnorm(q, 0, amplitude)/sqrt(n)
true_index= sample(p,q)
beta = numeric(p)
beta[true_index]=beta_q
mu = as.matrix(x) %*%  beta
y = mu + rnorm(n)
# Fit the Gaussian mirror model
fit = gm(y, x, ncores=4)
```







