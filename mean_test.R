mean_test <- function(X1, X2){
  check_sse<- function(X, split='False'){
    d <- dim(X)[1]
    n <- dim(X)[2]
    n1 <- as.integer(ceiling(n/2))
    n2 <- n - n1
    n12 <- list(n1, n2)
    r <- min(n2-1, d)
    
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
    
    Sd <- t(Xcdm[[1]]) %*% Xcdm[[2]] %*% t(Xcdm[[2]]) %*% Xcdm[[1]] / (n1 - 1) / (n2 - 1)
    eig <- eigen(Sd)
    cdmval <- eig$values[1:r]
    trSS <- sum(diag(Sd))
    
    psi <- numeric(n2)
    tau <- numeric(n2-1)
    khat <- n2-2
    
    for (i in 1:n2){
      if (i == 1){
        psi[i] <- trSS
      } else{
        term <- 0
        for (j in 1:(i-1)){
          term <- term + cdmval[j]
        }
        psi[i] <- trSS - term
      }
    }
    for (i in 1:(n2-1)){
      tau[i] <- psi[i+1] / psi[i] * (1 + i*sqrt(log(n) / n))
      if (tau[i] > 1){
        khat <- i - 1
        break
      }
    }
    khat <- min(khat, n2-2)
    
    return(list(khat=khat, psi=psi, tau=tau))
  }
  mean_test_nsse <- function(X1, X2){
    X <- list(X1, X2)
    n <- list()
    Mean <- list()
    trS <- list()
    
    for (i in 1:2){
      n <- append(n, dim(X[[i]])[2])
      Mean <- append(Mean, list(apply(X[[i]], 1, mean)))
      X_cent <- sweep(X[[i]], 1, Mean[[i]], '-')
      trS <- append(trS, list(sum(diag(t(X_cent) %*% X_cent)) / (n[[i]]-1)))
    }
    
    T_A <- norm(Mean[[1]]-Mean[[2]], type="2")^2 - trS[[1]]/n[[1]] - trS[[2]]/n[[2]]
    
    W <- function(X){
      d <- dim(X)[1]
      n <- dim(X)[2]
      X <- sweep(X, 1, apply(X, 1, mean), '-')
      M <- t(X) %*% X
      D <- numeric(n)
      for (i in 1:n){
        D[i] <- as.numeric(t(X[, i]) %*% X[, i])
      }
      term1 <- (n-1) * (n-2) * sum(diag(t(M) %*% M))
      term2 <- n * (n-1) * sum(D^2)
      term3 <- sum(D)^2
      w <- (term1 - term2 + term3) / n / (n-1) / (n-2) / (n-3)
      return(w)
    }
    term1 <- 2* (W(X[[1]]) / n[[1]] / (n[[1]]-1) + W(X[[2]]) / n[[2]] / (n[[2]]-1))
    term2 <- 4 * sum(diag(t(sweep(X[[1]], 1, apply(X[[1]], 1, mean), '-')) %*% sweep(X[[2]], 1, apply(X[[2]], 1, mean), '-') %*%
                            t(sweep(X[[2]], 1, apply(X[[2]], 1, mean), '-')) %*% sweep(X[[1]], 1, apply(X[[1]], 1, mean), '-'))) / 
      n[[1]] / n[[2]] / (n[[1]]-1) / (n[[2]]-1)
    K1_hat <- term1 + term2
    
    return (T_A / sqrt(K1_hat))
  }
  mean_test_sse <- function(X1, X2, sse_point){
    Dual <- function(X){
      d <- dim(X)[1]
      n <- dim(X)[2]
      r <- min(n-1, d)
      X <- sweep(X, 1, apply(X, 1, mean), '-')
      Sd <- t(X) %*% X / (n - 1)
      eig <- eigen(Sd)
      dualval <- eig$values[1:r]
      vec <- eig$vectors
      
      dualvec <- matrix(0, d, r)
      for (i in 1:r){
        dualvec[, i] <- X %*% vec[, i] / sqrt((n - 1) * dualval[i])
      }
      return ((list(values=dualval, vectors=dualvec)))
    }
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
    Psi <- function(X, split='False'){
      d <- dim(X)[1]
      n <- dim(X)[2]
      n1 <- as.integer(ceiling(n/2))
      n2 <- n - n1
      n12 <- list(n1, n2)
      r <- min(n2-1, d)
      
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
      
      Sd <- t(Xcdm[[1]]) %*% Xcdm[[2]] %*% t(Xcdm[[2]]) %*% Xcdm[[1]] / (n1 - 1) / (n2 - 1)
      eig <- eigen(Sd)
      cdmval <- eig$values[1:r]
      trSS <- sum(diag(Sd))
      
      psi <- numeric(n2)
      tau <- numeric(n2-1)
      khat <- n2-2
      
      for (i in 1:n2){
        if (i == 1){
          psi[i] <- trSS
        } else{
          term <- 0
          for (j in 1:(i-1)){
            term <- term + cdmval[j]
          }
          psi[i] <- trSS - term
        }
      }
      return(psi)
    }
    
    X <- list(X1, X2)
    K <- unlist(sse_point)
    d <- dim(X1)[1]
    n <- list()
    psi <- list()
    nrmval <- list()
    nrmvec <- list()
    nrmvec_self <- list()
    dualvec <- list()
    
    for (i in 1:2){
      n <- append(n, dim(X[[i]])[2])
      psi <- append(psi, list(Psi(X[[i]])))
      
      nrm <- NRM(X[[i]], K[i])
      nrmval <- append(nrmval, list(nrm$values))
      nrmvec <- append(nrmvec, list(nrm$vectors))
      nrmvec_self <- append(nrmvec_self, list(nrm$selfs))
      
      dual <- Dual(X[[i]])
      dualvec <- append(dualvec, list(dual$vectors))
    }
    
    
    term1 <- 0
    for (i in 1:2){
      term1_1 <- 0
      for (s in 1:n[[i]]){
        for (t in s:n[[i]]){
          if (s != t){
            term1_2 <- 0
            for (j in 1:K[i]){
              term1_2 <- term1_2 + as.numeric(t(nrmvec_self[[i]][, s, j]) %*% X[[i]][, s]) * as.numeric(t(nrmvec_self[[i]][, t, j]) %*% X[[i]][, t])
            }
            term1_1 <- term1_1 + as.numeric(t(X[[i]][, s]) %*% X[[i]][, t]) - term1_2
          }
        }
      }
      term1 <- term1 + 2*term1_1/n[[i]]/(n[[i]]-1)
    }
    
    term2 <- 0
    for (s in 1:n[[1]]){
      for (t in 1:n[[2]]){
        term2_1 <- numeric(d)
        term2_2 <- numeric(d)
        for (j in 1:K[1]){
          term2_1 <- term2_1 + as.numeric(t(nrmvec_self[[1]][, s, j]) %*% X[[1]][, s]) * nrmvec[[1]][, j]
        }
        for (j in 1:K[2]){
          term2_2 <- term2_2 + as.numeric(t(nrmvec_self[[2]][, t, j]) %*% X[[2]][, t]) * nrmvec[[2]][, j]
        }
        
        term2 <- term2 + as.numeric(t(X[[1]][, s] - term2_1) %*% (X[[2]][, t] - term2_2))
      }
    }
    term2 <- 2 * term2 / n[[1]] / n[[2]]
    
    T_ast <- term1 - term2
    
    
    term1 <- 0
    for (i in 1:2){
      term1 <- term1 + 2*psi[[i]][K[i] + 1] / n[[i]] / (n[[i]]-1)
    }
    A <- list()
    for (i in 1:2){
      term <- diag(0, d)
      for (j in 1:K[i]){
        term <- term + as.matrix(dualvec[[i]][, j]) %*% t(as.matrix(dualvec[[i]][, j]))
      }
      A <- append(A, list(diag(1, d) - term))
    }
    term2 <- 4 * sum(diag(t(sweep(X[[1]], 1, apply(X[[1]], 1, mean), '-')) %*% A[[1]] %*% sweep(X[[2]], 1, apply(X[[2]], 1, mean), '-') 
                          %*% t(sweep(X[[2]], 1, apply(X[[2]], 1, mean), '-')) %*% A[[2]] %*% sweep(X[[1]], 1, apply(X[[1]], 1, mean), '-'))) /
      n[[1]] / n[[2]] / (n[[1]]-1) / (n[[2]]-1)
    K_hat <- term1 + term2
    
    return(T_ast / sqrt(K_hat))
    
  }
  
  sse_point <- c(check_sse(X1)$khat, check_sse(X2)$khat)
  if (any(sse_point == 0)){
    print("NSSE model")
    ans <- mean_test_nsse(X1, X2)
  } else {
    print("SSE model")
    ans <- mean_test_sse(X1, X2, sse_point)
  }
  return(ans)
}
