rm(list = ls())

library(crossEstimation)

source("baselines.R")

n = 100
p = 500
tau = 10
eps = 0.2

X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, eps)
Y = rnorm(n) + X[,1] + tau * W

tau.hat = ate.glmnet(X, Y, W, lambda.choice = "lambda.min")
print(paste0("(", round(tau.hat$tau - 1.96 * sqrt(tau.hat$var), 2), ", ", round(tau.hat$tau + 1.96 * sqrt(tau.hat$var), 2), ")"))

tau.hat = ate.bloniarz(X, Y, W, lambda.choice = "lambda.min")
print(paste0("(", round(tau.hat$tau - 1.96 * sqrt(tau.hat$var), 2), ", ", round(tau.hat$tau + 1.96 * sqrt(tau.hat$var), 2), ")"))

tau.hat = ate.simple(X, Y, W)
print(paste0("(", round(tau.hat$tau - 1.96 * sqrt(tau.hat$var), 2), ", ", round(tau.hat$tau + 1.96 * sqrt(tau.hat$var), 2), ")"))

reps = replicate(100, {
	
X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, eps)
Y = rnorm(n) + X[,1] + tau * W

tau.hat.1 = ate.glmnet(X, Y, W, lambda.choice = "lambda.min")
print(paste0("cross est: (", round(tau.hat.1$tau - 1.96 * sqrt(tau.hat.1$var), 2), ", ", round(tau.hat.1$tau + 1.96 * sqrt(tau.hat.1$var), 3), ")"))
tau.hat.2 = ate.bloniarz(X, Y, W, lambda.choice = "lambda.min")
print(paste0("residual: (", round(tau.hat.2$tau - 1.96 * sqrt(tau.hat.2$var), 2), ", ", round(tau.hat.2$tau + 1.96 * sqrt(tau.hat.2$var), 3), ")"))
c(as.numeric(tau.hat.1), as.numeric(tau.hat.2))

})

print(mean(abs(reps[1,] - tau) / sqrt(reps[2,]) > 1.96))
print(mean(abs(reps[3,] - tau) / sqrt(reps[4,]) > 1.96))