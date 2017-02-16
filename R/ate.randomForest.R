ate.randomForest = function(X, Y, W, nodesize = 20, conf.level=.9) {
                      	
  if (prod(W %in% c(0, 1)) != 1) {                    	
    stop("Treatment assignment W must be encoded as 0-1 vector.")
  }
                        	
  nobs = nrow(X)
  pobs = ncol(X)
    
  yhat.0 = rep(NA, nobs)
  yhat.1 = rep(NA, nobs)
  
  if(length(unique(Y)) > 2) {  
  	
    rf.0 = randomForest::randomForest(X[W==0,], Y[W==0], nodesize = nodesize)
    rf.1 = randomForest::randomForest(X[W==1,], Y[W==1], nodesize = nodesize)

    yhat.0[W==0] = randomForest::predict(rf.0)
    yhat.0[W==1] = randomForest::predict(rf.0, newdata = X[W==1,])
    yhat.1[W==1] = randomForest::predict(rf.1)
    yhat.1[W==0] = randomForest::predict(rf.1, newdata = X[W==0,])
    
  } else {
  	
    rf.0 = randomForest::randomForest(X[W==0,], factor(Y)[W==0], nodesize = nodesize)
    rf.1 = randomForest::randomForest(X[W==1,], factor(Y)[W==1], nodesize = nodesize)

    yhat.0[W==0] = randomForest::predict(rf.0, type = "prob")[,2]
    yhat.0[W==1] = randomForest::predict(rf.0, newdata = X[W==1,], type = "prob")[,2]
    yhat.1[W==1] = randomForest::predict(rf.1, type = "prob")[,2]
    yhat.1[W==0] = randomForest::predict(rf.1, newdata = X[W==0,], type = "prob")[,2]
    
  }
  
  yhat.bar = (sum(W == 1) * yhat.0 + sum(W == 0) * yhat.1) / nobs
  
  tau.hat = mean((Y - yhat.bar)[W==1]) - mean((Y - yhat.bar)[W==0])
  var.hat = var((Y - yhat.bar)[W==1]) / sum(W == 1) +
    var((Y - yhat.bar)[W==0]) / sum(W == 0)
  ci=c(tau.hat-qnorm(1-(1-conf.level)/2)*sqrt(var.hat),tau.hat+qnorm(1-(1-conf.level)/2)*sqrt(var.hat))
  list(tau=tau.hat, var=var.hat, conf.int=ci, conf.level=conf.level)
}  
