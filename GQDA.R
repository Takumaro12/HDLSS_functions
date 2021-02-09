GQDA <- function(X, y, x0){
  x0 <- as.matrix(x0)
  n0 <- dim(x0)[2]
  label <- levels(as.factor(y))
  nclass <- length(label)
  d <- dim(X)[1]
  cls <- list()
  
  for (l in 1:n0){
    Y <- list()
    for (i in label){
      Xcls <- X[, y==i]
      n <- dim(Xcls)[2]
      Mean <- apply(Xcls, 1, mean)
      Xcls <- sweep(Xcls, 1, Mean, '-')
      trS <- sum(diag(t(Xcls) %*% Xcls / (n-1)))
      Y <- append(Y, list(d*norm(x0[, l] - Mean, type='2')^2/trS + d*log(trS) - d/n))
    }
    Y <- unlist(Y)
    cls <- unlist(append(cls, list(max(which(Y == min(Y))))))
  }
  
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
  return(list(class=ans))
}