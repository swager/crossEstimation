ate.bloniarz = function(X, Y, W, lambda.choice=c("lambda.min", "lambda.1se")) {
	
	fit.0 = cv.glmnet(X[W==0,], Y[W==0])
	fit.1 = cv.glmnet(X[W==1,], Y[W==1])
	
	resid.0 = Y[W==0] - predict(fit.0, newx=X[W==0,], s = lambda.choice)
	resid.1 = Y[W==1] - predict(fit.1, newx=X[W==1,], s = lambda.choice)
	
	tau.hat =  mean(resid.1) - mean(resid.0) +
	  mean(predict(fit.1, newx=X, s = lambda.choice)) -
	  mean(predict(fit.0, newx=X, s = lambda.choice))
	
	df.0 = sum(coef(fit.0, s = lambda.choice) != 0)
	df.1 = sum(coef(fit.1, s = lambda.choice) != 0)
	
	var.hat = mean(resid.0^2) / max(1, sum(W == 0) - df.0) +
	  mean(resid.1^2) / max(1, sum(W == 1) - df.1)
	  
	data.frame(tau=tau.hat, var=var.hat)
}

ate.simple = function(X, Y, W) {

	tau.hat = mean(Y[W==1]) - mean(Y[W==0])
	var.hat = var(Y[W==1]) / sum(W==1) + var(Y[W==0]) / sum(W==0)
	data.frame(tau=tau.hat, var=var.hat)
	
}
