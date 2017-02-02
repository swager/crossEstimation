\name{ate.randomForest}
\alias{ate.randomForest}
\title{
Average treatment effect estimation for randomized experiments with high dimensional covariates using random forest estimation
}
\description{
Use random forest to adjust covariates and estimate average treatment effects and confidence intervals
}
\usage{
ate.randomForest(X, Y, W, nodesize = 20, conf.level = 0.9)
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
  \item{nodesize}{
Node size for random forest. Default 20
}
\item{conf.level}{Confidence level for intervals. Default 0.9}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns list containing the following
\item{tau }{Average treatment effect estimate}
\item{var }{Variance estimate}
\item{conf.int }{Confidence interval for true tau}
\item{conf.level }{Confidence level}
}
\references{
S Wager, W Du, J Taylor, and R Tibshirani. "High-dimensional regression adjustments in randomized experiments". PNAS, November 8, 2016  vol. 113  no. 45  12673-12678
}
\author{
Stefan Wager, Wenfei Du, Jonathan Taylor, Rob Tibshirani
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
