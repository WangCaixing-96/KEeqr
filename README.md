# KEeqr
This is an R package for "Estimation of Conditional Extremiles in Reproducing Kernel Hilbert Space"


## Installation

You can use `devtools` to directly install the latest version from Github

```R
# install.packages("devtools")
devtools::install_github("WangCaixing-96/KEeqr")
```

## Example

Here is a toy example:

```r
## R Packages loading
library(KEeqr)
library(kernlab)
library(VGAM)

## Extremile generating function
xpareto <- function(tau,scale,shape)
{
  f <- function(x,scale,shape,tau)
  {
    qpareto(x,scale,shape)*J_tau(tau,x)
  }

  res = integrate(f, lower = 0, upper = 1, shape=shape, scale=scale, tau=tau, rel.tol=5*10^8*.Machine$double.eps)$value
  return(res)
}

## Model Settings and Data Generalization
n = 200			# number of obsevations
d = 3			# number of covariates dimensions
k_1 = as.integer(n^(0.1))			# n-k_1 is the end point of extremile interval
k_n = floor(3*n^(1/2))			#  n-k_n is the start point of extremile interval
tau_extreme = 0.99			# Extreme extremile level
r = 1 # index of controling the homoscedastic or heteroscedastic cases
x = matrix(runif(n*d), ncol=d) # distribution of covarites
beta = c(runif(d))/10  # model parameters
f_x = apply(2 * pi * x, 2, sin)%*%beta # main function
epsilo = rpareto(n, scale = 1, shape = 3) # error terms
y = f_x + (1 + r * apply(x, 1, mean)) * epsilo # model setting
true_extremile = f_x + (1 + apply(x, 1, mean)) * xpareto(tau.extreme, scale=1, shape=3) # the true extremile in level tau_extremile
```

Then we first estimate Condition CDF by the quantile regression process

```R
## Parameters Setting
C = 2 # C=1/lambda*n controls the smoothness of the quantile model
s = 30 #number of the quantile levels
kernel = rbfdot(sigma=1) # The standard RBF kernel
##Dual Optimization via Quadratic solver  
CDF_hat = CDF_QR(x, y, kernel, C, s) # the estimated CDF of every sample

```

Next, you can calculate the estimated extremile by  
```R
## Extremile Estimation


lambda_2n = 0.01 # Parameter controls the smoothness of the extremile model
kernel = rbfdot(sigma=1) # The standard RBF kernel
extremile_hat = KEE_qr(x, y, lambda_2n, kernel, C, s, k_1, k_n, tau_extreme, CDF_hat)
```

you can also calculate the ordinary estimated extremile without extrapolation by 
```R
## Ordinary Extremile Estimation 


ord_extremile_hat = OKE_qr(x, y, lambda_2n, kernel, C, s, tau_extreme, CDF_hat)
```

Finally, we can evaluate the models by Abias and RMSE of the estimates
```R
## Abias and RMSE calculation  

Abias_extra = mean(abs(extremile_hat - true_extremile)) # Absolute bias of estimated extremile with extrapolation
RMSE_extra = sqrt(mean((extremile_hat - true_extremile)^2)) # Square root of MSE of estimated extremile with extrapolation

Abias_ord = mean(abs(extremile_hat - true_extremile)) # Absolute bias of estimated extremile without extrapolation
RMSE_ord = sqrt(mean((extremile_hat - true_extremile)^2)) # Square root of MSE of estimated extremile without extrapolation
```











