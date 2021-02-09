cov_test <- function(X1, X2){
  NRM <- function(X){
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
    
    for (i in 1:r){
      nrmval[i] <- dualval[i] - (sum(diag(Sd)) - sum(dualval[1:i])) / (n - i - 1)
      nrmvec[, i] <- X %*% dualvec[, i] / sqrt((n - 1) * nrmval[i])
    }
    return ((list(values=nrmval, vectors=nrmvec)))
  }
  Delta <- function(X, split='False'){
    d <- dim(X)[1]
    n <- dim(X)[2]
    n1 <- as.integer(ceiling(n/2))
    n2 <- n - n1
    n12 <- list(n1, n2)
    r <- min(n2-1, d)
    Xcent <- sweep(X, 1, apply(X, 1, mean), '-')
    trS <- sum(diag(t(Xcent) %*% Xcent)) / (n-1)
    
    if (split=='False'){
      index <- c(1:n)
      Xcdm <- list(X[, 1:n1], X[, (n1+1):n])
    } else if (split=='True'){
      pi <- matrix(1/n, 1, n)
      index <- sample(c(1:n), n, replace=FALSE, pi)
      Xcdm <- list(X[, index[1:n1]], X[, index[(n1+1):n]])
    }
    Mean <- list(apply(Xcdm[[1]], 1, mean), apply(Xcdm[[2]], 1, mean))
    for (i in 1:2){
      Xcdm[[i]] <- sweep(Xcdm[[i]], 1, Mean[[i]], '-')
    }
    
    SDcross <- t(Xcdm[[1]]) %*% Xcdm[[2]]  / sqrt(n1 - 1) / sqrt(n2 - 1)
    svd <- svd(SDcross)
    cdmval <- svd$d[1]
    trSS <- sum(diag(SDcross %*% t(SDcross)))
    delta <- trSS - cdmval^2
    
    return(list(delta=delta))
  }
  
  X12 <- list(X1, X2)
  d <- dim(X1)[1]
  n12 <- list()
  nrmval <- list()
  nrmvec <- list()
  delta <- list()
  kappa <- list()
  for (i in 1:2){
    n12 <- append(n12, dim(X12[[i]])[2])
    nrm <- NRM(X12[[i]])
    delkappa <- Delta(X12[[i]])
    nrmval <- append(nrmval, list(nrm$values[1]))
    nrmvec <- append(nrmvec, list(nrm$vectors[, 1]))
    delta <- append(delta, list(delkappa$delta))
    
    Xcent <- sweep(X12[[i]], 1, apply(X12[[i]], 1, mean), '-')
    trS <- sum(diag(t(Xcent) %*% Xcent)) / (n12[[i]]-1)
    kappa <- append(kappa, trS-nrmval[[i]])
  }
  
  eta <- sqrt(delta[[1]]) / nrmval[[1]] + sqrt(delta[[2]]) / nrmval[[2]]
  T_NR <- ((nrmval[[1]] - nrmval[[2]])^2 + 
             2*nrmval[[1]]*nrmval[[2]]*(1-min(1, as.numeric(t(nrmvec[[1]]) %*% nrmvec[[2]])^2))^(1+eta)) /
    (2*nrmval[[1]]^2/(n12[[1]]-1) + 2*nrmval[[2]]^2/(n12[[2]]-1))
  
  kappa_ast <- (kappa[[1]]/kappa[[2]] + kappa[[2]]/kappa[[1]])/2
  delta_ast <- (delta[[1]]/delta[[2]] + delta[[2]]/delta[[1]])/2
  T_IYA <- kappa_ast * delta_ast * T_NR
  
  p_value <- pchisq(T_IYA, 1, lower.tail = FALSE)
  return(list(T_IYA=T_IYA, p_value=p_value))
}