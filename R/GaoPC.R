#' @param data is a total data set
#' @param m is the number of first layer principal component
#' @param pc is the number of second layer principal component
#' @return AG1,AG2,DG1,DG2,SigmahatG1,SigmahatG2
#' @export
#' @examples
#' GaoPC(data=ISE,m=3)
GaoPC=function(data,m){
  X=scale(data)
  n=nrow(X)
  SigmahatG1=cor(X)
  eig3<-eigen(SigmahatG1)
  lambda1hat =eig3$values[1:m]
  ind<-order(lambda1hat,decreasing=T)
  lambda1hat<-lambda1hat[ind]
  Q<- eig3$vectors
  Q<-Q[,ind]
  AG1<-Q[,1:m]
  hG1 <- diag(AG1 %*% t(AG1))
  DG1 <- diag(SigmahatG1 - hG1)
  p<-ncol(X)
  pc=2
  F1hat=X%*%AG1
  F1star<-F1hat/sqrt(n)
  SigmahatG2=cov(F1star)
  eig4<-eigen(SigmahatG2)
  lambda2hat = eig4$values[1:pc]
  ind<-order(lambda2hat,decreasing=T)
  lambda2hat<-lambda2hat[ind]
  Q<- eig4$vectors
  Q<-Q[,ind]
  AG2<-Q[,1:pc]
  hG2 <- diag(AG2 %*% t(AG2))
  DG2 <- diag(SigmahatG2 - hG2)
  Fhat=F1star%*%AG2
  XGhat=Fhat%*%t(AG2)%*%t(AG1)
  sGhat=cov(XGhat)
  hG3 <- diag(t(t(AG2)%*%t(AG1))%*%(t(AG2)%*%t(AG1)))
  DG3 <- diag(sGhat - hG3)                 #
  return(list(AG1=AG1,AG2=AG2,DG3=DG3,sGhat=sGhat))}
