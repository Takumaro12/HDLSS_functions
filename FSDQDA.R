
FS_DQDA <- function(X, y, x0, gamma=0.5){
  x0 <- as.matrix(x0)
  n0 <- dim(x0)[2]
  label <- levels(as.factor(y))
  nclass <- length(label)
  d <- dim(X)[1]
  
  cls <- list()
  X12 <- list()
  n12 <- list()
  Mean <- list()
  diagS <- list()
  for (i in label){
    X12 <- append(X12, list(X[, y==i]))
  }
  for (i in 1:nclass){
    n12 <- append(n12, list(dim(X12[[i]])[2]))
    Mean <- append(Mean, list(apply(X12[[i]], 1, mean)))
    Xcent <- sweep(X12[[i]], 1, Mean[[i]], '-')
    diagS <- append(diagS, list(diag(Xcent %*% t(Xcent) / (n12[[i]]-1))))
  }
  n_min <- min(unlist(n12))
  xi <- sqrt(log(d) / n_min)
  theta <- numeric(d)
  for (i in 1:nclass){
    for (j in 1:nclass){
      if (i != j){
        theta <- theta + ((Mean[[i]] - Mean[[j]])^2 + diagS[[i]])/(nclass*(nclass-1)*diagS[[j]])
      }
    }
  }
  theta <- theta - 1
  index <- which(theta > xi^gamma)
  
  for (l in 1:n0){
    Y <-list()
    for (i in 1:nclass){
      Y <- append(Y, sum((x0[index, l] - Mean[[i]][index])^2 / diagS[[i]][index] - 1 / n12[[i]] + log(diagS[[i]][index])))
    }
    
    argmin <- function(X){
      len <- length(X)
      m <- min(unlist(X))
      result <- list()
      for (l in 1:len){
        if (X[l] == m){
          result <- append(result, l)
        }
      }
      return(unlist(result))
    }
    cls <- unlist(append(cls, list(max(argmin(Y)))))
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
  return(list(class=ans,variables=index))
}
