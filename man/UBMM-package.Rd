\name{UBMM}
\alias{Uniform-Beta Mixture Model}
\docType{package}
\title{
  Boosted EM algorithm: Uniform-Beta Mixture Model
}
\description{
\code{UBMM} applies a Boosted EM algorithm to fit a two-point mixture model of Uniform and Beta distributions. The results of evaluating \code{UBMM} returns a list of Weights, Beta parameters, and iterations to converge, respectively. The method \code{UBMM} is built in C++, which is quite fast and stable. The convergence depends on the initial values for weights and Beta shape parameters though. 
}
\usage{
UBMM(x, w, a, precision, Iterations=10000L)
}
\arguments{
  \item{x}{A numeric vector which ranges between 0 and 1.}
  \item{w}{A vector of initial weights for the Uniform and Beta distributions in the mixture model.}
  \item{a}{Initial parameters for the Beta distribution.}
  \item{Precision}{Tolerance for convergence of the EM algorithm.}
  \item{Iterations}{Maximum number of iterations in the EM algorithm. Default is 10000L.}
}
\details{
The \code{UBMM} provides a boosted EM algorithm to fit a two-point mixture of Uniform and Beta distributions. Instead of optimizing the shape parameters for the Beta distribution at each Maximization step, the boosted EM algorithm alternately updates those shape parameters by using the gradient method. This boosted EM algorithm is stable and much faster than the classical EM algorithm. 

The package depends R packages \code{Rcpp} and \code{BH}. Users need install the two packages before installing the \code{UBMM} package. 
}
\author{
Chong Ma <chongm@email.sc.edu>
}
\references{
Chakraborty, P. et.al. (2018) \emph{Asymptotic Conditional Update for Mixture Models Used in Large Scale Inference}. Statistics & Probability Letters.
}
\seealso{}
\examples{
\dontrun{
## generate a mixture of Uniform and Beta
## distribution with shape parameters 0.5
x=c(runif(9500),rbeta(500,0.5,0.5))
UBMM(x,c(0.5,0.5),c(1,2),1e-8)
}
}
\keyword{ package }
