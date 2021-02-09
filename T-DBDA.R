T_DBDA <- function(X, y, x0, sse_point=NULL, centering='True', random='False'){
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
  
  x0 <- as.matrix(x0)
  
  if (centering=='True'){
    X_mean <- apply(X, 1, mean)
    X <- sweep(X, 1, X_mean, "-")
    x0 <- sweep(x0, 1, X_mean, "-")
  }
  
  n0 <- dim(x0)[2]
  label <- levels(as.factor(y))
  nclass <- length(label)
  d <- dim(X)[1]
  
  cls <- list()
  X12 <- list()
  n12 <- list()
  Mean <- list()
  trS <- list()
  nrmval <- list()
  nrmvec <- list()
  nrmvec_self <- list()
  
  #Divide data matrix X by label y.
  for (i in label){
    X12 <- append(X12, list(X[, y==i]))
  }
  
  M <- list()
  #Change vector from list. 
  sse_point <- unlist(sse_point)
  #Result of check_see for each classes append M.
  
  if (mode(sse_point)!="numeric"){
    for (k in 1:nclass){
      M <- unlist(append(M, check_sse(X12[[k]])))
    }
    #If sse_point had been entered, this value append M.
  } else {
    M <- sse_point
  }
  
  for (k in 1:nclass){
    #sample size for each classes.
    n12 <- append(n12, dim(X12[[k]])[2])
    #mean vector size for each classes. 
    Mean <- append(Mean, list(apply(X12[[k]], 1, mean)))
    #centering
    Xcent <- sweep(X12[[k]], 1, Mean[[k]], '-')
    #value of tr(S) for each classes.
    trS <- append(trS, sum(diag(t(Xcent) %*% Xcent)) / (n12[[k]]-1))
    
    nrm <- NRM(X12[[k]], M[k])
    #eigenvalue by NRM
    nrmval <- append(nrmval, list(nrm$values))
    #eigenvector by NRM
    nrmvec <- append(nrmvec, list(nrm$vectors))
    #transformed eigenvector 
    nrmvec_self <- append(nrmvec_self, list(nrm$selfs))
  }
  dim(nrmvec_self[[2]])
  
  #xj project h(j,r)
  x_tilda <- list(matrix(0,M[1],n12[[1]]),matrix(0,M[2],n12[[2]])) 
  for (i in 1:2){
    for (r in 1:M[i]){
      for (j in 1:n12[[i]]) {
        x_tilda[[i]][r,j] <- as.numeric(t(X12[[i]][, j]) %*% nrmvec_self[[i]][, j, r])
      }
    }
  }
  
  #mean(x_tilda)*h_tilda
  x_varh <- list(numeric(d), numeric(d))
  for (i in 1:2){
    for (s in 1:M[i]){
      x_varh[[i]] <- x_varh[[i]] + mean(x_tilda[[i]][s,]) * nrmvec[[i]][, s]
    }
  }
  
  
  term2 <- list(0, 0)
  
  for (k in 1:2){#class
    for (r in 1:M[k]){#sse point
      for (i in 1:n12[[k]]){#sample
        for (j in 1:i){#sample
          if (i != j){
            term2[[k]] <- term2[[k]] + as.numeric(x_tilda[[k]][r,i] * x_tilda[[k]][r,j]) / n12[[k]] / (n12[[k]] - 1)
          }
        }
      }
    }
  }
  
  
  term1 <- list(0, 0)
  classifier <- numeric(n0)
  
  for (l in 1:n0){
    #calaculate W(x0)
    DBDA <- as.numeric(t(x0[, l] - (Mean[[1]]+Mean[[2]])/2) %*% (Mean[[2]]-Mean[[1]])) - trS[[1]]/(2*n12[[1]]) + trS[[2]]/(2*n12[[2]])
    
    for (k in 1:2){#class
      for (r in 1:M[k]){#sse point
        x_var <- apply(x_tilda[[k]], 1, mean)
        term1[[k]] <- term1[[k]] + as.numeric(t(x0[, l]) %*% nrmvec[[k]][, r]) * (x_var[r] - 1/2*as.numeric(t(nrmvec[[k]][, r]) %*% (Mean[[3-k]] - x_varh[[3-k]])))
      }
    }
    #evaluate l th sample
    classifier[l] <- DBDA + term1[[1]] - term1[[2]] - term2[[1]] + term2[[2]]
    #labeling test data
    if (classifier[l] < 0){
      cls <- unlist(append(cls, 1))
    } else if (classifier[l] >= 0){
      cls <- unlist(append(cls, 2))
    }
  }
  
  #reform answer to match original label y.
  ans <- list()
  for (l in 1:n0){
    for (i in 1:nclass){
      if (cls[l] == i){
        ans <- unlist(append(ans, list(label[[i]])))
      }
    }
  }
  
  if (mode(y) == "numeric"){
    ans <- as.numeric(ans)
  }
  
  return (list(khat = M, class=ans, classifier=classifier))
}