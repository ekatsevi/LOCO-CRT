rm(list = ls())
machine = "local"
source("setup.R")
n = 1000                          
p = 500                           
sigma = 1                         
A = 5
k = 50
rho = 0.5

precomp_filename = sprintf("%s/precomp/precomp_rho_%d_p_%d.Rda", base_dir, round(100*rho), p)
load(precomp_filename)

X = matrix(rnorm(n*p), n, p)%*%Sigma_chol
nonzeros = round(seq(20, p-20, 
                     length.out = k)) 
beta = matrix(0, p, 1)
beta[nonzeros] = A/sqrt(n)*sample(x = c(-1,1),     
                          size = length(nonzeros), 
                          replace = TRUE)
y = X%*%beta + rnorm(n = n, mean = 0, sd = sigma)

########### TEST MANUAL CV-LASSO FUNCTION #############

nfolds = 10
output = CV_lasso(X, y, nfolds)
fit_glmnet = cv.glmnet(X, y,  foldid = rep_len(1:nfolds, length.out = n), standardize = FALSE, keep = TRUE)

# check mean CV errors
plot(output$cvm[1:length(fit_glmnet$cvm)], fit_glmnet$cvm, pch = 19)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

# check sd CV errors
plot(output$cvsd[1:length(fit_glmnet$cvsd)], fit_glmnet$cvsd, pch = 19)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

# check out-of-sample predictions
plot(output$outfold_predictions[,10], fit_glmnet$fit.preval[,10], pch = 19)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

# check in-sample predictions
plot(output$infold_predictions[,10], predict(fit_glmnet, X, s = fit_glmnet$lambda)[,10], pch = 19)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

########### TEST AGREEMENT BETWEEN LASSO AND LOCO-LASSO ##########

nfolds = 10
fit_glmnet = cv.glmnet(X, y,  foldid = rep_len(1:nfolds, length.out = n), standardize = FALSE, keep = TRUE)
output_full = CV_lasso(X, y)

for(j in c(nonzeros[1], setdiff(which(output_full$first_entries < Inf), nonzeros)[1])){
  output_loco = CV_lasso(X[,-j], y)
  plot(log(fit_glmnet$lambda), fit_glmnet$cvm, type = "l", lwd = 2, col = 'black')
  points(log(fit_glmnet$lambda), output_loco$cvm[1:length(fit_glmnet$lambda)], type = "l", lwd = 2, col = 'red', lty = 2)
  if(output_full$first_entries[j] < Inf){
    abline(v = log(fit_glmnet$lambda[output_full$first_entries[j]]), lty = 2, lwd = 2)
  } else{
    print("This variable did not enter lasso path at all.")
  }  
}

########### TEST AR(1) CHOLESKY DECOMPOSITION AND INVERSE ##########

p = 100
rho = 0.5
Sigma = AR_cov(p, rho)
Sigma_chol = AR1_chol(rho, p)
max(abs(Sigma - t(Sigma_chol)%*%Sigma_chol))

Sigma_inv = AR1_cov_inv(rho, p)
max(abs(Sigma_inv%*%Sigma - diag(rep(1,p))))
max(abs(Sigma%*%Sigma_inv - diag(rep(1,p))))