# Construct AR(1) covariance matrix
AR_cov = function(p, rho){
  Sigma = outer(1:p, 1:p, function(i, j)(rho^(abs(i-j))))
  return(Sigma)
}

# Construct Cholesky factor of AR(1) covariance matrix
#  https://blogs.sas.com/content/iml/2018/10/03/ar1-cholesky-root-simulation.html
AR1_chol = function(rho, p){
  R = matrix(0, p,p)                # allocate p x p matrix
  R[1,] = rho^(0:(p-1))             # formula for 1st row
  c = sqrt(1 - rho*rho)             # scaling factor: c^2 + rho^2
  R2 = c * R[1,]                    # formula for 2nd row
  for(j in 2:p){                    # shift elements in 2nd row for remaining rows
      R[j, j:p] = R2[1:(p-j+1)]
  }
  R[abs(R) < .Machine$double.eps] = 0
  R_sparse = Matrix(R, sparse = TRUE)
  return(R_sparse)
}

# construct AR(1) precision matrix
AR1_cov_inv = function(rho, p){
  sparseMatrix(i = c(1:p, 1:(p-1), 2:p), 
               j = c(1:p, 2:p, 1:(p-1)), 
               x = c(1,rep(1+rho*rho, p-2), 1, rep(-rho, 2*(p-1)))/(1-rho*rho))
}