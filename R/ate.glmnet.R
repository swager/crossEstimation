ate.glmnet = function(X, Y, W,
                      alpha=1,
                      nfolds=NULL,
                      method=c("joint", "separate"),
                      lambda.choice=c("lambda.min", "lambda.1se")) {
    
  method = match.arg(method)
  lambda.choice = match.arg(lambda.choice)
                      	
  if (prod(W %in% c(0, 1)) != 1) {                    	
    stop("Treatment assignment W must be encoded as 0-1 vector.")
  }
                        	
  nobs = nrow(X)
  pobs = ncol(X)
  
  if (is.null(nfolds)) {
  	nfolds = floor(max(3, min(10, sum(W==0)/5, sum(W==1)/5)))
  }
  
  # fold ID for cross-validation; balance treatment assignments
  foldid = rep(NA, nobs)
  foldid[W==0] = sample(rep(seq(nfolds), length = sum(W==0)))
  foldid[W==1] = sample(rep(seq(nfolds), length = sum(W==1)))
  
  # unique identifier for treatment status x foldid
  bucket = foldid + W * nfolds                    	
  
  # compute mean response for each bucket
  X.mean = matrix(NA, 2 * nfolds, pobs)
  Y.mean = rep(NA, 2 * nfolds)
  for (bucket.id in unique(bucket)) {
  	X.mean[bucket.id,] = colMeans(X[bucket == bucket.id, , drop = FALSE])
    Y.mean[bucket.id] = mean(Y[bucket == bucket.id])
  }
  
  # center data within each bucket
  X.centered = t(sapply(1:nobs, function(nn) {
  	X[nn,] - X.mean[bucket[nn],]
  }))
  Y.centered = Y - Y.mean[bucket]
  
  # run glmnet on centered data
  beta.0 = NA
  beta.1 = NA
  
  if (method == "joint") {
  
    regr.matrix = cbind(X.centered, (2 * W - 1) * X.centered)
    regr.fit = my.cv.glmnet(regr.matrix, Y.centered, foldid = foldid, alpha = alpha, intercept = FALSE, standardize = FALSE, lambda.choice = lambda.choice)
    betas.all = regr.fit$cv.betas[,-1]
    beta.0 = betas.all[,1:pobs] - betas.all[,pobs + (1:pobs)]
    beta.1 = betas.all[,1:pobs] + betas.all[,pobs + (1:pobs)]
   
  } else if (method == "separate") {
  	
    fit.0 = my.cv.glmnet(X.centered[W==0,], Y.centered[W==0], foldid = foldid[W==0], alpha = alpha, intercept = FALSE, standardize = FALSE, lambda.choice = lambda.choice)
    beta.0 = fit.0$cv.betas[,-1]
    
    fit.1 = my.cv.glmnet(X.centered[W==1,], Y.centered[W==1], foldid = foldid[W==1], alpha = alpha, intercept = FALSE, standardize = FALSE, lambda.choice = lambda.choice)
    beta.1 = fit.1$cv.betas[,-1]

  } else {
  	
  	stop("Invalid method")
  
  }
  
  # compute other statistics
  counts = outer(1:nfolds, 0:1, FUN = Vectorize(function(a, b) sum(bucket == a + nfolds * b)))
  
  beta.bar = t(sapply(1:nfolds, function(ID) {
  	(counts[ID, 2] * beta.0[ID,] + counts[ID, 1] * beta.1[ID,]) / (counts[ID, 2] + counts[ID, 1])
  }))
  
  X.fold.mean = matrix(NA, nfolds, pobs)
  for (ID in unique(foldid)) {
  	X.fold.mean[ID,] = colMeans(X[foldid == ID, , drop = FALSE])
  }
  
  # compute fold-wise ATE estimates
  tau = sapply(1:nfolds, function(ID) {
  	Y.mean[ID + nfolds] - Y.mean[ID] +
  	  sum((X.fold.mean[ID,] - X.mean[ID + nfolds,]) * beta.1[ID,]) -
  	  sum((X.fold.mean[ID,] - X.mean[ID,]) * beta.0[ID,])
  })
  tau.hat = mean(tau)
  
  # compute variance 	
  var.fold = sapply(1:nfolds, function(ID) {
  	var.0 = var(Y[W == 0 & foldid == ID] - X[W == 0 & foldid == ID,] %*% beta.bar[ID,])  / counts[ID, 1]
  	var.1 = var(Y[W == 1 & foldid == ID] - X[W == 1 & foldid == ID,] %*% beta.bar[ID,])  / counts[ID, 2] 
    var.0 + var.1
  })
  var.hat = mean(var.fold) / nfolds
    
  data.frame(tau=tau.hat, var=var.hat)
}  