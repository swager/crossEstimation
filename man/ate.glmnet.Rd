\name{ate.glmnet}
\alias{ate.glmnet}
\title{
Average treatment effect estimation for randomized experiments with high dimensional covariates
}
\description{
Use cross estimation procedure with lasso in glmnet to estimate average treatment effects and intervals
}
\usage{
ate.glmnet(X, Y, W, alpha = 1, nfolds = NULL, method = c("joint", "separate"), lambda.choice = c("lambda.1se", "lambda.min"))
}
\arguments{
  \item{X}{
Data matrix with covariates
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
There are two cross validation choices, \code{lambda.choice = "lambda.1se"} uses one standard deviation rule and \code{lambda.choice = "lambda.min"} uses minimum
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns list containing the following
\item{tau }{Average treatment effect estimate}
\item{var }{Variance estimate}
}
\references{
S Wager, W Du, J Taylor, and R Tibshirani. "High-dimensional regression adjustments in randomized experiments". In: arXiv preprint arXiv:1607.06801 (2016).
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
