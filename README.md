# UBMM
Boosted EM algorithm: Uniform-Beta Mixture Model
-------------------------------------------------

**UBMM** applies a Boosted EM algorithm to fit a two-point mixture model of Uniform and Beta distributions. The results of evaluating **UBMM** returns a list of Weights, Beta parameters, and iterations to converge, respectively. The method **UBMM** is built in C++, which is quite fast and stable. The convergence depends on the initial values for weights and Beta shape parameters though. 

## Install UBMM 
* Install Rcpp, BH packages in R 
* Run `R CMD INSTALL UBMM` in cmd to install the package UBMM

## Arguments
-- | ---
x | A numeric vector which ranges between 0 and 1.
w | A vector of initial weights for the Uniform and Beta distributions in the mixture model.
a | Initial parameters for the Beta distribution.
Precision | Tolerance for convergence of the EM algorithm.
Iterations | Maximum number of iterations in the EM algorithm. Default is 10000L

## Run an example
```r
## generate a mixture of Uniform and Beta
## distribution with shape parameters 0.5
x=c(runif(9500),rbeta(500,0.5,0.5))
UBMM(x,c(0.5,0.5),c(1,2),1e-8)
```



