# LOCO-CRT, based on lasso with cross-validation
LOCO_CRT = function(X, y, vars, X_cond_means, X_cond_sds, loco_crt_methods, nfolds = 10){
  # extract problem dimensions
  n = nrow(X)
  m = ncol(X)
  num_vars = length(vars)
  num_methods = nrow(loco_crt_methods)
  if(num_methods == 0){
    return(NULL)
  }

  p_values = matrix(NA, num_methods, num_vars)
  
  if(num_vars == m){
    # run CV lasso on full set of variables, 
    # keeping track of active set
    cat(sprintf("Fitting lasso on all variables...\n"))
    output = CV_lasso(X,y,nfolds)
    first_entries = output$first_entries
    cvm_full = output$cvm
    cvsd_full = output$cvsd
    outsample_predictions_full = output$outfold_predictions
    insample_predictions_full = output$infold_predictions
    active = first_entries <= which.min(cvm_full)
    cat(sprintf("%d active variables in full lasso.\n", sum(active)))
  } else{
    active = !logical(num_vars)
  }
  

  # iterate over variables j to be tested
  for(var_idx in 1:num_vars){
    j = vars[var_idx]
    if(active[var_idx]){
      # fit the LOCO lasso
      cat(sprintf("Refitting lasso holding out variable %d...\n", j))
      foldid = rep_len(1:nfolds, length.out = n)
      fit_loco = cv.glmnet(x = X[,-j], y = y, foldid = foldid, standardize = FALSE, keep = TRUE)
      # extract information from LOCO fit
      cvm = fit_loco$cvm
      cvsd = fit_loco$cvsd
      outsample_predictions = fit_loco$fit.preval
      insample_predictions = predict(fit_loco, X[,-j,drop=FALSE], s = fit_loco$lambda)
    } else{
      # extract information from full fit
      cvm = cvm_full
      cvsd = cvsd_full
      outsample_predictions = outsample_predictions_full
      insample_predictions = insample_predictions_full
    }
    
    for(method_idx in 1:num_methods){
      inner_method = loco_crt_methods$inner_method[method_idx]
      stopping_rule = loco_crt_methods$stopping_rule[method_idx]
      fit_type = loco_crt_methods$fit_type[method_idx]
      if(stopping_rule == "1se"){
        lambda_idx_min = which.min(cvm)
        lambda_idx_loco = min(which(cvm <= cvm[lambda_idx_min] + cvsd[lambda_idx_min]))
      }
      if(stopping_rule == "min"){
        lambda_idx_loco = which.min(cvm)
      }
      if(fit_type == "in_sample"){
        r = y - insample_predictions[,lambda_idx_loco]
      }
      if(fit_type == "out_sample"){
        r = y - outsample_predictions[,lambda_idx_loco]      
      }
      p_values[method_idx, var_idx] = 
        get_univariate_CRT_pvalue(r, X[,j], NA, X_cond_means[,var_idx], X_cond_sds[var_idx], inner_method)
    }
  }
  return(p_values)
}

full_CRT = function(X, y, vars, X_cond_means, X_cond_sds, full_crt_methods){
  n = nrow(X)
  p = ncol(X)
  num_vars = length(vars)
  num_methods = nrow(full_crt_methods)
  if(num_methods == 0){
    return(NULL)
  }
  method_names = full_crt_methods$method_name

  stopifnot(length(unique(full_crt_methods$MC_reps)) == 1)
  MC_reps = full_crt_methods$MC_reps[1]
  
  p_values = matrix(NA, num_methods, num_vars, dimnames = list(method_names, vars))
  original_test_stats = array(0, dim = c(num_vars, num_methods), dimnames = list(vars, method_names))
  resampled_test_stats = array(0, dim = c(num_vars, num_methods, MC_reps), dimnames = list(vars, method_names, 1:MC_reps))

  # iterate over variables to test
  for(var_idx in 1:num_vars){
    j = vars[var_idx]
    cat(sprintf("Testing variable %d...\n", j))

    # compute original test statistics
    original_test_stats[var_idx,] = full_CRT_test_stats(X, y, j, full_crt_methods)

    # compute resampled test statistics
    for(MC_rep in 1:MC_reps){
      cat(sprintf("Working on Monte Carlo rep %d out of %d...\n", MC_rep, MC_reps))
      X_tilde_j = rnorm(n = n,  mean = X_cond_means[,var_idx], sd = X_cond_sds[var_idx])
      X_tilde = X
      X_tilde[,j] = X_tilde_j
      resampled_test_stats[var_idx,,MC_rep] = full_CRT_test_stats(X_tilde, y, j, full_crt_methods)
    }
    
    # compute p-values 
    for(var_idx in 1:num_vars){
      for(method in method_names){
        p_values[method, var_idx] = mean(c(Inf, resampled_test_stats[var_idx, method,]) >= original_test_stats[var_idx, method])
      }
    }
  }
  return(p_values)
}

marginal_CRT = function(X, y, vars, X_cond_means, X_cond_sds, marginal_crt_methods){
  num_vars = length(vars)
  num_methods = nrow(marginal_crt_methods)
  if(num_methods == 0){
    return(NULL)
  }
  method_names = marginal_crt_methods$method_name
  p_values = matrix(NA, num_methods, num_vars, dimnames = list(method_names, vars))
  for(method_idx in 1:num_methods){
    inner_method = marginal_crt_methods$inner_method[method_idx]
    for(var_idx in 1:num_vars){
      j = vars[var_idx]
      r = y
      p_values[method_idx, var_idx] = 
        get_univariate_CRT_pvalue(r, X[,j], NA, X_cond_means[,var_idx], X_cond_sds[var_idx], inner_method)
    }
  }
  return(p_values)
}

HRT = function(X, y, vars, X_cond_means, X_cond_sds, hrt_methods){
  n = nrow(X)
  p = ncol(X)
  num_vars = length(vars)
  num_methods = nrow(hrt_methods)
  if(num_methods == 0){
    return(NULL)
  }
  method_names = hrt_methods$method_name
  p_values = matrix(NA, num_methods, num_vars, dimnames = list(method_names, vars))
  
  for(method_idx in 1:num_methods){
    learning_method = hrt_methods$learning_method[method_idx]
    inner_method = hrt_methods$inner_method[method_idx]
    holdout_prop = hrt_methods$holdout_prop[method_idx]
    stopping_rule = hrt_methods$stopping_rule[method_idx]
    holdout_num = floor(holdout_prop*n)
    holdout = 1:holdout_num
    holdin = (holdout_num+1):n
    fit_holdin = cv.glmnet(x = X[holdin,], y = y[holdin], standardize = FALSE)
    beta_hat_holdin = coef(fit_holdin, s = sprintf("lambda.%s", stopping_rule))[2:(p+1), 1, drop = FALSE]
    for(var_idx in 1:num_vars){
      j = vars[var_idx]
      beta_hat = beta_hat_holdin
      if(learning_method == "lasso"){
        if(beta_hat[j] == 0){
          p_values[var_idx] = NA # screening
          next
        }
      }
      if(learning_method == "relaxed_lasso"){
        X_j_centered = X[holdin,j]-X_cond_means[holdin,var_idx]
        beta_hat[j] = sum(X_j_centered*(y[holdin,]-X[holdin,-j]%*%beta_hat[-j]))/sum(X_j_centered*X_j_centered)
      }
      r = y[holdout] - X[holdout,-j]%*%beta_hat[-j] - X_cond_means[holdout,var_idx]*beta_hat[j]
      p_values[method_idx, var_idx] = 
        get_univariate_CRT_pvalue(r, X[holdout,j], beta_hat[j], X_cond_means[holdout,var_idx], X_cond_sds[var_idx], inner_method)
    }
  }
  return(p_values)
}

oracle_HRT = function(X, y, beta, vars, X_cond_means, X_cond_sds, oracle_hrt_methods){
  n = nrow(X)
  p = ncol(X)
  num_vars = length(vars)
  num_methods = nrow(oracle_hrt_methods)
  if(num_methods == 0){
    return(NULL)
  }
  method_names = oracle_hrt_methods$method_name
  p_values = matrix(NA, num_methods, num_vars, dimnames = list(method_names, vars))
  
  for(method_idx in 1:num_methods){
    inner_method = oracle_hrt_methods$inner_method[method_idx]
    holdout_prop = oracle_hrt_methods$holdout_prop[method_idx]
    holdout_num = floor(holdout_prop*n)
    holdout = 1:holdout_num
    for(var_idx in 1:num_vars){
      j = vars[var_idx]
      if(beta[j] == 0){
        p_values[method_idx, var_idx] = 1
      } else{
        r = y[holdout] - X[holdout,-j]%*%beta[-j] - X_cond_means[holdout,var_idx]*beta[j]
        p_values[method_idx, var_idx] = 
          get_univariate_CRT_pvalue(r, X[holdout,j], beta[j], X_cond_means[holdout,var_idx], X_cond_sds[var_idx], inner_method)
      }
    }
  }
  return(p_values)
}