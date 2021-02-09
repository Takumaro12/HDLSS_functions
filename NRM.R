NRM <- function(X){
  d <- dim(X)[1]
  n <- dim(X)[2]
  r <- min(n-2, d)
  Mean <- apply(X, 1, mean)
  X <- sweep(X, 1, Mean, '-')
  Sd <- t(X) %*% X / (n - 1)
  eig <- eigen(Sd)
  dualval <- eig$values[1:r]
  dualvec <- eig$vectors
  
  nrmval <- numeric(r)
  nrmvec <- matrix(0, d, r)
  nrmscore <- matrix(0, n, r)
  
  for (i in 1:r){
    nrmval[i] <- dualval[i] - (sum(diag(Sd)) - sum(dualval[1:i])) / (n - i - 1)
    nrmvec[, i] <- X %*% dualvec[, i] / sqrt((n - 1) * nrmval[i])
    
    for (j in 1:n){
      nrmscore[j, i] <- dualvec[j, i] * sqrt(n * nrmval[i])
    }
  }
  return(list(values=nrmval, vectors=nrmvec, scores=nrmscore))
}