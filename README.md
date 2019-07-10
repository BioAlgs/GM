# Gaussian Mirror for High-dimensional Controled Variable Selection


##Installation
```R
install_github("DeveloperName/PackageName")
```



##Example 1
This example illustrates the basic usage of the GM package. The GM method is designed for any design matrix. For simplicity, we will use synthetic data constructed from a linear model such that the response only depends on a small fraction of the variables.


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
```
