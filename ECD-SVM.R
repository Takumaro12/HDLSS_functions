ECD-SVD <-function(X,p1,p2){
  n<-dim(X)[2]
  d<- dim(X)[1]
  n1<-ceiling(n/2)
  n2<-n-n1
  u<-2*n1*n2/((n1-1)*(n2-1)*n*(n-1))
  Y<-rbind(X,matrix(0,1,n))
  V1<-function(k,x){
    if (floor(k/2)>=n1){
      x[,(floor(k/2)-n1+1):floor(k/2)]
    } 
    else {
      cbind(x[,1:floor(k/2)],x[,(floor(k/2)+n2+1):n])
    }
  }
  V2<-function(k,x)
  {
    if (floor(k/2)<=n1)
    {
      x[,(floor(k/2)+1):(floor(k/2)+n2)]
    } 
    else
    {
      cbind(x[,1:(floor(k/2)-n1)],x[,(floor(k/2)+1):n])
    }
  }
  H1<-function(k)
  {
    apply(V1(k,Y),1, mean)
  }
  H2<-function(k)
  {
    apply(V2(k,Y),1, mean)
  }
  S<-c(3:(2*n-1))
  M1<-sapply(S,H1)
  M2<-sapply(S,H2)
  qq<-function(i,j)
  {
    (Y[1:p1,j]-M2[1:p1,i+j-2])%*%(t(Y[(p1+1):d,i]-M1[(p1+1):d,i+j-2])%*%(Y[(p1+1):d,j]-M2[(p1+1):d,i+j-2]))%*%(t(Y[1:p1,i]-M1[1:p1,i+j-2]))
  }
  
  pp<-function(i,j)
  {
    (Y[(p1+1):d,j]-M2[(p1+1):d,i+j-2])%*%(t(Y[1:p1,i]-M1[1:p1,i+j-2]))%*%(Y[1:p1,j]-M2[1:p1,i+j-2])%*%t(Y[(p1+1):d,i]-M1[(p1+1):d,i+j-2])
  }
  
  Q1 <- matrix(0,p1,p1)
  Q2 <- matrix(0,p2,p2)
  
  for (i in 1:(n-1)) {
    for(j in (i+1):n ){
      Q1 <- Q1 + (qq(i,j)+t(qq(i,j)))/2
      Q2 <- Q2 + (pp(i,j)+t(pp(i,j)))/2
    }
  }
  
  model <- eigen(u*Q1)
  val <- model$values
  U <- model$vectors
  
  model <- eigen(u*Q2)
  V <- model$vectors
  
  return(list(d=val,u=U,v=V))
}