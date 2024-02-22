#' @param data is a total  data set
#' @param m is the number of principal component
#' @return Apro,Dpro,Sigmahatpro
#' @export cov
#' @examples
#' PPC(data=data,m=m)
PPC=function(data,m){
  X=scale(data)
  n=nrow(X)
  p=ncol(X)
  P=as.matrix(diag(c(0,1),n,n))
  Xpro=scale(P%*%X)
  Sigmahatpro<-cov(Xpro)
  eig<-eigen(Sigmahatpro)
  lambdahat =eig$values[1:m]
  ind<-order(lambdahat,decreasing=T)
  lambdahat<-lambdahat[ind]
  Q <- eig$vectors
  Q<-Q[,ind]
  Qhat<-Q[,1:m]
  Apro <- matrix(0, nrow = p, ncol = m)
  for (j in 1:m) {Apro[, j] <- sqrt(lambdahat[j]) * Qhat[, j]}; Apro
  hpro <- diag(Apro %*% t(Apro))
  Dpro <- diag(Sigmahatpro - hpro)
  return(list(Apro=Apro,Dpro=Dpro,Sigmahatpro=Sigmahatpro))}

