\name{ate.glmnet}
\alias{ate.glmnet}
\title{
Average treatment effect estimation for randomized experiments with high dimensional covariates
}
\description{
Use the cross estimation procedure with lasso in glmnet to estimate average treatment effects and confidence intervals
}
\usage{
ate.glmnet(X, Y, W, alpha = 1, nfolds = NULL, method = c("joint", "separate"), lambda.choice = c("lambda.1se", "lambda.min"), conf.level=0.9)
}
\arguments{
  \item{X}{
Data matrix with covariates, one observation per row
}
  \item{Y}{
Outcome vector for assigned treatment for each observation
}
  \item{W}{
Treatment vector for each observation
}
  \item{alpha}{
Elastic net parameter
}
  \item{nfolds}{
Number of folds
}
  \item{method}{
There are two methods, \code{method = "joint"} estimates coefficients from joint model and \code{method = "separate"} estimates coefficients for each treatment separately
}
  \item{lambda.choice}{
There are two cross validation choices, \code{lambda.choice = "lambda.1se"} uses the one standard deviation rule and \code{lambda.choice = "lambda.min"} uses the minimum
}
\item{conf.level}{Confidence level for confidence level. Default 0.9}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns list containing the following
\item{tau }{Average treatment effect estimate}
\item{var }{Variance estimate}
\item{conf.int }{Confidence interval for true tau}
\item{conf.level }{Confidence level}
}
\references{
S Wager, W Du, J Taylor, and R Tibshirani. "High-dimensional regression adjustments in randomized experiments". PNAS, November 8, 2016  vol. 113  no. 45  12673–12678
}
\author{
Stefan Wager, Wenfei Du, Jonathan Taylor, Rob Tibshirani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
  # run 1000 times and calculate coverage
  cover <- 0
  for (i in 1:1000) {
    x <- matrix(rnorm(n * p, xmean, xsigma), n, p)
    T <- (runif(n) < 0.2)
    mu <- (ymean1 + x %*% theta1) * T + (ymean0 + x %*% theta0) * (1 - T)
    epsC <- rnorm(n, 0, sigma)
    epsT <- rnorm(n, 0, sigma)
    eps <- epsT * T + epsC * (1 - T)
    yobs <- mu + eps
    res <- ate.glmnet(x, yobs, T, alpha = 1, nfolds = 10, method = "joint", lambda.choice = "lambda.min")
    cover <- cover + (res$conf.int[1] < tau & tau < res$conf.int[2])
  }
  cover
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
