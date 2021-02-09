ASPCA <- function(X, lam, MeanZero=F){
  d <- dim(X)[1]
  n <- dim(X)[2]
  if (MeanZero){
    r <- min(n-1, d)
    Sd <- t(X) %*% X / n
    eig <- eigen(Sd)
    dualval <- eig$values[1:r]
    dualvec <- eig$vectors
    nrmval <- numeric(r)
    nrmvec <- matrix(0, d, r)
    aspca <- matrix(0, d, r)
    for (i in 1:r){
      nrmval[i] <- dualval[i] - (sum(diag(Sd)) - sum(dualval[1:i])) / (n - i)
      nrmvec[, i] <- X %*% dualvec[, i] / sqrt(n * nrmval[i])
      ord <- order(abs(nrmvec[, i]), decreasing=T)
      cri <- 0
      for (j in 1:d){
        cri <- cri + nrmvec[ord[j], i]^2
        aspca[ord[j], i] <- nrmvec[ord[j], i]
        if (cri >= 1){
          break
        }
      }
    }
  } else {
    r <- min(n-2, d)
    X <- sweep(X, 1, apply(X, 1, mean), '-')
    Sd <- t(X) %*% X / (n - 1)
    eig <- eigen(Sd)
    dualval <- eig$values[1:r]
    dualvec <- eig$vectors
    nrmval <- numeric(r)
    nrmvec <- matrix(0, d, r)
    aspca <- matrix(0, d, r)
    for (i in 1:r){
      nrmval[i] <- dualval[i] - (sum(diag(Sd)) - sum(dualval[1:i])) / (n - i - 1)
      nrmvec[, i] <- X %*% dualvec[, i] / sqrt((n - 1) * nrmval[i])
      ord <- order(abs(nrmvec[, i]), decreasing=T)
      cri <- 0
      for (j in 1:d){
        cri <- cri + nrmvec[ord[j], i]^2
        aspca[ord[j], i] <- nrmvec[ord[j], i]
        if (cri >= lam){
          break
        }
      }
    }
  }
  return(list(values=nrmval, vectors=nrmvec, asvectors=aspca))
}