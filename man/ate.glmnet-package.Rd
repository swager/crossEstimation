\name{ate.glmnet-package}
\alias{ate.glmnet-package}
\docType{package}
\title{
\packageTitle{ate.glmnet}
}
\description{
\packageDescription{ate.glmnet}
}
\details{
\tabular{ll}{
Package: \tab ate.glmnet \cr
Type: \tab Package \cr
Version: \tab 0.0 \cr
Date: \tab 2016-12-13 \cr
Imports: \tab glmnet \cr
}
}
\author{
\packageAuthor{ate.glmnet}
}
\references{
S Wager, W Du, J Taylor, and R Tibshirani. "High-dimensional regression adjustments in randomized experiments". In: arXiv preprint arXiv:1607.06801 (2016)
}
\keyword{ package }
\examples{
 # simulation with Gaussian covariates
  n <- 200
  p <- 500
  p0 <- p
  theta0 <- rep(1/p0, p0)
  theta1 <- rep((1 + .1)/p0, p0)
  ymean0 <- 4
  ymean1 <- 3
  xmean <- 1
  xsigma <- 1
  sigma <- 1
  x <- matrix(rnorm(n * p, 1, xsigma), n, p)
  x0 <- x[,1:p0]
  T <- (runif(n) < 0.5)
  mu <- ymean1 * T + ymean0 * (1 - T) + x0 %*% theta1 * T + x0 %*% theta0 * (1 - T)
  epsC <- rnorm(n, 0, sigma)
  epsT <- rnorm(n, 0, sigma)
  eps <- epsT * T + epsC * (1 - T)
  yobs <- mu + eps
  ate.glmnet(x, yobs, T, alpha = 1, nfolds = 10, method = "joint", lambda.choice = "lambda.min")
}
