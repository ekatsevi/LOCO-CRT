run_one_experiment = function(experiment_name){
  ########## READ INPUT ##########
  # set up all the simulation parameters
  input_filename = sprintf("input_file_%s.R", experiment_name)
  input_mode = "experiment"
  source(input_filename, local = TRUE)
  
  ########## SET UP ARRAYS TO HOLD RESULTS ##########
  # various dimensions
  num_global_error_rates = 
    length(global_error_rates)
  num_methods = length(methods)
  num_parameters = nrow(parameters)
  
  # p-values and knockoff statistics
  statistics = array(NA, 
                     dim = c(num_methods, p, num_parameters), 
                     dimnames = list(methods, 1:p, NULL))
  # power and type-I error of global tests
  metrics = array(NA, 
                  dim = c(num_methods, 2, num_global_error_rates, num_parameters), 
                  dimnames = list(methods, c("power", "type_I_error"), global_error_rates, NULL))
  # rejections made by global tests
  rejections = array(FALSE, 
                     dim = c(num_methods, num_global_error_rates, p, num_parameters), 
                     dimnames = list(methods, global_error_rates, NULL, NULL))
  
  ########## LOOP OVER PARAMETER SETTINGS ##########
  for(index in 1:nrow(parameters)){
    A = parameters$A[index]
    k = parameters$k[index]
    rep = parameters$rep[index]
    rho = parameters$rho[index]
    cat(sprintf("Working on parameter index %d out of %d, corresponding to rep %d for A = %0.1f, rho = %0.1f, k = %d...\n", 
                index, num_parameters, index, A, rho, k))
    
    ########## PROBLEM SETUP ###########
    
    Sigma_chol = AR1_chol(rho, p)
    Sigma_inv = AR1_cov_inv(rho, p)
    
    # E[X_j|X_-j] = X[,-j]%*%t(cond_means[j,,drop=FALSE])
    cond_means = sparseMatrix(i = c(1,p, 2:(p-1), 2:(p-1)), 
                              j = c(1,p-1,1:(p-2),2:(p-1)), 
                              x = c(rep(rho,2),rep(-Sigma_inv[1,2]/Sigma_inv[2,2],2*(p-2))))
    
    # sd[X_j|X_-j] = cond_sds[j]
    cond_sds = sqrt(1/diag(Sigma_inv))
    
    # knockoffs equicorrelation constant
    s = min(1, 2*(1-rho)/(1+rho))
    
    # quantities needed for knockoff generation  
    knockoffs_cov = 2*s*sparseMatrix(i = 1:p, j = 1:p, x = rep(1, p)) - s*s*Sigma_inv
    knockoffs_chol = chol(knockoffs_cov)
    knockoffs_mean = sparseMatrix(i = 1:p, j = 1:p, x = rep(1, p)) - s*Sigma_inv
    
    # create design matrix with correlated columns
    X = matrix(rnorm(n*p), n, p)%*%Sigma_chol
    if(nrow(knockoff_methods) > 0){
      cat(sprintf("Creating knockoffs...\n"))
      knockoff_realizations = max(knockoff_methods$knockoff_realization)
      X_k_list = vector("list", knockoff_realizations)
      for(realization in 1:knockoff_realizations){
        X_k_list[[realization]] = X%*%knockoffs_mean + matrix(rnorm(n*p),n,p)%*%knockoffs_chol
      }
    }
    # calculate mean and standard deviation of each variable, conditional on the others
    X_cond_means = matrix(NA, n, p)
    for(j in 1:p){
      X_cond_means[,j] = as.numeric(X[,-j]%*%t(cond_means[j,,drop=FALSE]))
    }
    X_cond_sds = cond_sds
    
    ### create beta
    
    # equally spaced non-nulls, excluding edges
    buffer = ceiling(p/200)
    nonzeros = round(seq(buffer, p-buffer, length.out = k)) 
    beta = matrix(0, p, 1)
    # coefficients have magnitude A/sqrt(n) and random signs
    beta[nonzeros] = A/sqrt(n)*sample(x = c(-1,1),     
                                      size = length(nonzeros), 
                                      replace = TRUE)
    # create response y
    y = X%*%beta + rnorm(n = n, mean = 0, sd = sigma)
    
    ########## RUN METHODS ###########
    
    if(nrow(hrt_methods) > 0){
      cat(sprintf("Working on HRT...\n"))
      statistics[hrt_methods$method_name,,index] = 
        HRT(X, y, 1:p, X_cond_means, X_cond_sds, hrt_methods)
    }
    
    if(nrow(loco_crt_methods) > 0){
      cat(sprintf("Working on LOCO-CRT...\n"))
      statistics[loco_crt_methods$method_name,,index] = 
        LOCO_CRT(X, y, 1:p, X_cond_means, X_cond_sds, loco_crt_methods)
    }
    
    if(nrow(full_crt_methods) > 0){
      cat(sprintf("Working on full CRT...\n"))
      statistics[full_crt_methods$method_name,,index] = 
        full_CRT(X, y, 1:p, X_cond_means, X_cond_sds, full_crt_methods)
    }
    
    if(nrow(marginal_crt_methods) > 0){
      cat(sprintf("Working on marginal CRT...\n"))
      statistics[marginal_crt_methods$method_name,,index] = 
        marginal_CRT(X, y, 1:p, X_cond_means, X_cond_sds, marginal_crt_methods)
    }
    
    if(nrow(oracle_hrt_methods) > 0){
      cat(sprintf("Working on oracle HRT...\n"))
      statistics[oracle_hrt_methods$method_name,,index] = 
        oracle_HRT(X, y, beta, 1:p, X_cond_means, X_cond_sds, oracle_hrt_methods)
    }
    
    if(nrow(knockoff_methods) > 0){
      cat(sprintf("Working on knockoffs...\n"))
      for(method_idx in 1:nrow(knockoff_methods)){
        method_name = knockoff_methods$method_name[method_idx]
        test_stat = knockoff_methods$test_stat[method_idx]
        realization = knockoff_methods$knockoff_realization[method_idx]
        cat(sprintf("Computing %s for knockoff realization %d...\n", test_stat, realization))
        W = do.call(sprintf("stat.%s", test_stat), list(X, X_k_list[[realization]], as.numeric(y)))
        statistics[method_name, ,index] = W
      }
    }
    
    ########## EVALUATE TYPE-I AND TYPE-II ERRORS ###########
    
    # evaluate errors for global tests
    for(global_error_rate in global_error_rates){
      # FDR
      if(global_error_rate == "FDR"){
        # apply BH or knockoff filter to these p-values or knockoff statistics
        for(method in methods){
          if(grepl("knockoffs", method)){
            # run knockoff filter
            W = statistics[method, ,index]
            thres = knockoff.threshold(W, fdr=alpha_FDR, offset=1)
            selected = W >= thres
          } else{
            # run BH
            P = statistics[method, ,index]
            selected = logical(p)
            if(any(!is.na(P))){
              # to accommodate screening
              selected[!is.na(P)] = p.adjust(P[!is.na(P)], method = "BH") <= alpha_FDR
            }
          }
          rejections[method,"FDR",selected,index] = TRUE
          metrics[method, "power", "FDR", index] = sum(selected[nonzeros])/length(nonzeros)
          metrics[method, "type_I_error", "FDR", index] = sum(selected[-nonzeros])/max(1,sum(selected))
        }
      } 
      # FWER
      if(global_error_rate == "FWER"){
        # apply BH or knockoff filter to these p-values or knockoff statistics
        for(method in methods){
          if(grepl("knockoffs", method)){
            next
          } 
          # run BH
          P = statistics[method, ,index]
          selected = logical(p)
          if(any(!is.na(P))){
            # to accommodate screening
            selected[!is.na(P)]  = p.adjust(P[!is.na(P)], method = "bonferroni") <= alpha_FWER
          }
          rejections[method,"FWER",selected,index] = TRUE
          metrics[method, "power", "FWER", index] = sum(selected[nonzeros])/length(nonzeros)
          metrics[method, "type_I_error", "FWER", index] = as.numeric(sum(selected[-nonzeros]) > 0)
        }
      }
    }
  }

  ########## MASSAGE OUTPUT AND WRITE TO FILE ###########
  
  # create output directory, if necessary
  results_dir=sprintf("../results/%s", experiment_name)
  if(!dir.exists(results_dir)){
    dir.create(results_dir)
  }
  
  # statistics
  df = melt(statistics)
  names(df) = c("method", "variable", "index", "value")
  df = as_tibble(df)
  df = df %>% 
    mutate(A = parameters$A[index], 
           k = parameters$k[index],
           rho = parameters$rho[index],
           rep = parameters$rep[index]) %>%
    select(-index)
  statistics_filename = sprintf("../results/%s/statistics_%s.tsv", experiment_name, experiment_name)
  write_tsv(df, path = statistics_filename)
  
  # rejections
  df = melt(rejections)
  names(df) = c("method", "error_rate", "variable", "index", "rejection")
  df = as_tibble(df)
  df = df %>% mutate(A = parameters$A[index], 
                     k = parameters$k[index],
                     rho = parameters$rho[index],
                     rep = parameters$rep[index]) %>%
    select(-index)
  rejections_filename = sprintf("../results/%s/rejections_%s.tsv", experiment_name, experiment_name)
  write_tsv(df, path = rejections_filename)
  
  # metrics
  df = melt(metrics)
  names(df) = c("method", "metric", "error_rate", "index", "value")
  df = as_tibble(df)
  df = df %>% mutate(A = parameters$A[index], 
                     k = parameters$k[index],
                     rho = parameters$rho[index],
                     rep = parameters$rep[index]) %>%
    select(-index)
  metrics_filename = sprintf("../results/%s/metrics_%s.tsv", experiment_name, experiment_name)
  write_tsv(df, path = metrics_filename)
}