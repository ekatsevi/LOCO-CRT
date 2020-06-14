# helper function to do cross-validated lasso, keeping track of active set
CV_lasso = function(X, y, nfolds = 10){
  n = nrow(X)
  m = ncol(X)
  foldid = rep_len(1:nfolds, length.out = n)
  fit_full = glmnet(X, y, standardize = FALSE, keep = TRUE)
  infold_predictions = predict.glmnet(fit_full, X)
  lambda = fit_full$lambda
  # lambda = get_lambda(X, y)
  infold_beta = array(FALSE, c(nfolds, m, 100))
  CV_errors = matrix(NA, 100, nfolds)
  all_outfold_predictions = matrix(NA, n, 100)
  for(fold in 1:nfolds){
    fit_infold = glmnet(X[foldid != fold,], y[foldid != fold], lambda = lambda, standardize = FALSE)
    outfold_predictions = predict.glmnet(fit_infold, X[foldid == fold,])
    CV_errors[1:length(fit_infold$lambda),fold] = 
      apply(outfold_predictions, 2, function(predictions)(mean((y[foldid == fold]-predictions)^2)))
    infold_beta[fold, ,1:length(fit_infold$lambda)] = array(fit_infold$beta != 0, c(1, m, length(fit_infold$lambda)))
    all_outfold_predictions[foldid == fold, 1:length(fit_infold$lambda)] = outfold_predictions
  }
  cvm = apply(CV_errors, 1, mean)
  cvsd = apply(CV_errors, 1, sd)/sqrt(nfolds)
  lambda_idx = which.min(cvm)
  active = apply(infold_beta, c(2,3), any)
  first_entries = apply(active, 1, 
                        function(active_variable)(ifelse(any(active_variable), min(which(active_variable)), Inf)))
  output = c()
  output$first_entries = first_entries
  output$cvm = cvm
  output$cvsd = cvsd
  output$fit_full = fit_full
  output$infold_predictions = infold_predictions
  output$outfold_predictions = all_outfold_predictions
  return(output)
}

# helper function for univariate Gaussian problem
get_univariate_CRT_pvalue = function(r, X_j, beta_j, X_cond_mean, X_cond_sd, inner_method){
  n = length(r)
  if(inner_method == "quadratic"){
    r_std = r/sqrt(sum(r*r))
    X_j_std = (X_j - X_cond_mean)/sqrt(sum((X_j - X_cond_mean)*(X_j - X_cond_mean)))
    sq_dot_product = sum(r_std*X_j_std)^2
    pvalue = pf(1/(n-1)*(1/sq_dot_product-1), n-1, 1)
  }
  if(inner_method == "quadratic_fixed"){
    r_original = r - beta_j*(X_j-X_cond_mean)
    pvalue = pchisq(sum(r_original*r_original)/(beta_j*beta_j*X_cond_sd*X_cond_sd),
                    df = n,
                    ncp = sum(r*r)/(beta_j*beta_j*X_cond_sd*X_cond_sd))
  }
  return(pvalue)
}

# full CRT test statistics
full_CRT_test_stats = function(X, y, j, full_crt_methods){
  n = nrow(X)
  p = ncol(X)
  
  num_methods = length(full_crt_methods$method_name) 
  
  test_stats = numeric(num_methods)
  names(test_stats) = full_crt_methods$method_name
  
  if(any(full_crt_methods$learning_method == "lasso")){
    fit_lasso = cv.glmnet(x = X, y = y, standardize = FALSE)
  }
  if(any(full_crt_methods$learning_method == "relaxed_lasso")){
    fit_relaxed_lasso = cv.glmnet(x = X, y = y, standardize = FALSE, penalty.factor = as.numeric(1:p != j))
  }
  
  # compute test statistics
  for(index in 1:num_methods){
    learning_method = full_crt_methods$learning_method[index]
    stopping_rule = full_crt_methods$stopping_rule[index]
    test_stat = full_crt_methods$test_stat[index]
    method_name = full_crt_methods$method_name[index]
    if(learning_method == "lasso"){
      fit = fit_lasso
    }
    if(learning_method == "relaxed_lasso"){
      fit = fit_relaxed_lasso
    }
    beta_hat = coef(fit, s = sprintf("lambda.%s", stopping_rule))[2:(p+1), 1, drop = FALSE]
    r = y - X%*%beta_hat
    if(test_stat == "insample_loss"){
      test_stats[method_name] = -mean(r*r) 
    }
    if(test_stat == "cv_loss"){
      lambda_idx = which(fit$lambda == fit[sprintf("lambda.%s", stopping_rule)])
      test_stats[method_name] = -fit$cvm[lambda_idx]
    }
    if(test_stat == "abs_coef"){
      test_stats[method_name] = abs(beta_hat[j])
    }
  }
  return(test_stats)
}

# sequential rule to select lambda
select_lambda = function(err_estimates, num_increasing){
  num_lambda = length(err_estimates)
  err_estimates_next_few = numeric(num_lambda)
  err_estimates_next_few[1:(num_lambda-1)] = sapply(1:(num_lambda-1), 
                                                    function(idx)(min(err_estimates[(idx+1):min(idx+num_increasing, num_lambda)])))
  err_estimates_next_few[num_lambda] = Inf
  return(min(which(err_estimates <= err_estimates_next_few)))
}

# helper function to find sequence of lambda for glmnet
# https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
get_lambda = function(X, y, nlambda = 100){
  n = nrow(X)
  m = ncol(X)
  lambda_max <- max(abs(colSums(X*as.vector(y-mean(y)))))/n
  lambda.min.ratio <- ifelse(n < m, 0.01, 1e-04)
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*lambda.min.ratio), 
                              length.out = nlambda)), digits = 10)
  return(lambdapath)
}