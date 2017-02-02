\name{crossEstimation-package}
\alias{crossEstimation-package}
\docType{package}
\title{
\packageTitle{crossEstimation}
}
\description{
\packageDescription{crossEstimation}
}
\details{
\tabular{ll}{
Package: \tab crossEstimation \cr
Type: \tab Package \cr
Version: \tab 0.0 \cr
Date: \tab 2016-12-13 \cr
Imports: \tab glmnet \cr
}
}
\author{
\packageAuthor{crossEstimation}
}
\references{
S Wager, W Du, J Taylor, and R Tibshirani. "High-dimensional regression adjustments in randomized experiments". PNAS, November 8, 2016  vol. 113  no. 45  12673-12678

}
\keyword{ package }
\examples{
  # simulation with Gaussian covariates based on Figure 1 in reference paper
  set.seed(30)
  n <- 200
  p <- 500
  xmean <- 1
  xsigma <- 1
  sigma <- .1
  # set average treatment effect equal to one
  ymean0 <- 4
  ymean1 <- 3
  # set no heterogeneous treatment effects
  theta0 <- c(1, rep(0, p-1))
  theta1 <- c(1, rep(0, p-1))
  tau <- ymean1 - ymean0 + sum(xmean * theta1) - sum(xmean * theta0)
  # run 100 times and calculate coverage
  cover <- 0
  for (i in 1:100) {
    x <- matrix(rnorm(n * p, xmean, xsigma), n, p)
    T <- (runif(n) < 0.2)
    mu <- (ymean1 + x \%*\% theta1) * T + (ymean0 + x \%*\% theta0) * (1 - T)
    epsC <- rnorm(n, 0, sigma)
    epsT <- rnorm(n, 0, sigma)
    eps <- epsT * T + epsC * (1 - T)
    yobs <- mu + eps
    res <- ate.glmnet(x, yobs, T, alpha = 1, nfolds = 10, method = "joint", 
	  lambda.choice = "lambda.min")
    cover <- cover + (res$conf.int[1] < tau & tau < res$conf.int[2])
  }
  cover
}
