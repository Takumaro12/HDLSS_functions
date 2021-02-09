DT <- function(X, sse_point=NULL, random='False'){
  #sse check function
  check_sse<- function(X){
    d <- dim(X)[1]
    n <- dim(X)[2]
    n1 <- as.integer(ceiling(n/2))
    n2 <- n - n1
    n12 <- list(n1, n2)
    r <- min(n2-1, d)
    
    
    if (random=='False'){
      index <- c(1:n)
      Xcdm <- list(X[, 1:n1], X[, (n1+1):n])
    } else if (random=='True'){
      pi <- matrix(1/n, 1, n)
      index <- sample(c(1:n), n, replace=FALSE, pi)
      Xcdm <- list(X[, index[1:n1]], X[, index[(n1+1):n]])
    }
    Mean <- list(apply(Xcdm[[1]], 1, mean), apply(Xcdm[[2]], 1, mean))
    
    for (i in 1:2){
      Xcdm[[i]] <- sweep(Xcdm[[i]], 1, Mean[[i]], '-')
    }
    
    Sd <- t(Xcdm[[1]]) %*% Xcdm[[2]] %*% t(Xcdm[[2]]) %*% Xcdm[[1]] / (n1 - 1) / (n2 - 1)
    eig <- eigen(Sd)
    cdmval <- eig$values[1:r]
    trSS <- sum(diag(Sd))
    psi <- numeric(n2)
    tau <- numeric(n2-1)
    khat <- n2-2
    
    psi[1] <- trSS 
    for (j in 2:(r+1)) {
      psi[j] <- psi[j-1]- cdmval[j-1]
      tau[j-1] <- psi[j]/psi[j-1]
    }
    
    dis <- 0
    for (i in 1:(n2-1)) {
      if(dis < 1){
        dis <- tau[i] * (1 + i*sqrt(log(n) / n))
        khat <- i - 1
      }
      else break
    }
    khat <- min(khat, n2-2)
    return(khat)
  } 
  #NRM
  NRM <- function(X, sse_point){
    d <- dim(X)[1]
    n <- dim(X)[2]
    r <- min(n-2, d)
    X <- sweep(X, 1, apply(X, 1, mean), '-')
    Sd <- t(X) %*% X / (n - 1)
    eig <- eigen(Sd)
    dualval <- eig$values[1:r]
    dualvec <- eig$vectors
    
    nrmval <- numeric(r)
    nrmvec <- matrix(0, d, r)
    nrmvec_self <- array(0, dim=c(d, n, sse_point))
    
    for (i in 1:r){
      nrmval[i] <- dualval[i] - (sum(diag(Sd)) - sum(dualval[1:i])) / (n - i - 1)
      nrmvec[, i] <- X %*% dualvec[, i] / sqrt((n - 1) * nrmval[i])
    }
    
    c <- sqrt(n-1) / (n-2)
    for (i in 1:sse_point){
      u_hat <- dualvec[, i]
      for (j in 1:n){
        u_hat[j] <- -u_hat[j]/(n-1)
        nrmvec_self[, j, i] <- c * X %*% u_hat / sqrt(nrmval[i])
        u_hat <- dualvec[, i]
      }
    }
    return ((list(values=nrmval, vectors=nrmvec, selfs=nrmvec_self)))
  }
  
  if (mode(sse_point)!="numeric"){
    M <- check_sse(X)
    #If sse_point had been entered, this value append M.
  } else {
    M <- sse_point
  }
  
  d <- dim(X)[1]
  n <- dim(X)[2]
  nrm <- NRM(X,M)
  nrmval <- nrm$values
  nrmvec <- nrm$vectors
  nrmvec_self <- nrm$selfs
  
  #xj project h(j,r)
  x_tilda <- matrix(0,M,n)
  for (r in 1:M){
    for (j in 1:n) {
      x_tilda[r,j] <- as.numeric(t(X[, j]) %*% nrmvec_self[, j, r])
    }
  }
  
  Xtilda <- X - nrmvec[,1:M]%*%x_tilda
  return (Xtilda)
}