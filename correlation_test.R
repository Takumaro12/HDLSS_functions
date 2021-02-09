correlation_test <- function(X, x0){
  d <- dim(X)[1]
  n <- dim(X)[2]
  n1 <- as.integer(ceiling(n/2))
  n2 <- n - n1
  
  K <- c(3:(2*n-1))
  L <- length(K)
  Y <- array(0, dim=c(2, L, d))
  for (l in 1:L){
    V <- list()
    dv <- as.integer(floor(K[l]/2))
    
    if (dv >= n1){
      id <- c((dv-n1+1): dv)
      V <- append(V, list(id))
    } else{
      id <- append(c(1: dv), c((dv+n2+1): n))
      V <- append(V, list(id))
    }
    
    if (dv <= n1){
      id <- c((dv+1): (dv+n2))
      V <- append(V, list(id))
    } else{
      id <- append(c(1: (dv-n1)), c((dv+1): n))
      V <- append(V, list(id))
    }
    
    for (i in 1:2){
      Y[i, l, ] <- apply(X[, V[[i]]], 1, mean)
    }
  }
  
  u <- n1 * n2 / ((n1-1) * (n2-1))
  w <- 0
  t <- 0
  for (j in 1:n){
    for (i in 1:j){
      if (i != j){
        w <- w + (as.numeric((X[, i] - Y[1, (i+j-2), ]) %*% (X[, j] - Y[2, (i+j-2), ])))^2
        t <- t + as.numeric((X[, i] - Y[1, (i+j-2), ]) %*% (X[, j] - Y[2, (i+j-2), ])) * (x0[i] - y0[1, (i+j-2)]) * (x0[j] - y0[2, (i+j-2)])
      }
    }
  }
  W <-  2 * u / (n * (n - 1)) * w
  T_hat <- 2 * u / (n * (n - 1)) * t

  S0 <- var(x0)
  test <- T_hat / (S0 * sqrt(2*W) / n)
  p_value <-  pnorm(test, lower.tail = FALSE)
  
  return(list(test_stat=test, p_value=p_value))
}