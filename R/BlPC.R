#' @param data is a total data set
#' @param m is the number of principal component
#' @return ABr,ABc,DBr,DBc,SigmaB1hat,SigmaB2hat
#' @export
#' @examples
#' BlPC(data=ISE,m=3)
BlPC=function(data,m){
  Xrc=scale(data)
  n=nrow(Xrc)
  p=ncol(Xrc)
  Sigmahat=1/(n)*(Xrc-mean(Xrc))%*%t(Xrc-mean(Xrc))
  eig0<-eigen(Sigmahat)
  lambdar = eig0$values[1:p]
  ind<-order(lambdar,decreasing=T)
  lambdar<-lambdar[ind]
  Q1<-eig0$vectors
  Q1=Q1[,ind]
  Qr<- Q1[, 1:p]
  Ar0 <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {Ar0[, j] <- sqrt(lambdar[j]) * Qr[, j]}; Ar0
  hr <- diag(Ar0 %*% t(Ar0))
  Dr <- diag(Sigmahat - hr)
  Yr=scale(1/p*t(Xrc)%*%Ar0)
  Sigmahatw=1/(n)*t(Xrc-mean(Xrc))%*%(Xrc-mean(Xrc))
  eig1<-eigen(Sigmahatw)
  lambdac = eig1$values[1:m]
  ind<-order(lambdac,decreasing=T)
  lambdahatc<-lambdac[ind]
  Q2<-eig1$vectors
  Q2=Q2[,ind]
  Qc<- Q2[1:p, 1:m]
  Ac0 <- matrix(0, nrow = p, ncol =m)
  for (j in 1:m) {Ac0[, j] <- sqrt(lambdac[j]) * Qc[, j]}; Ac0
  hc <- diag(Ac0 %*% t(Ac0))
  Dc <- diag(Sigmahatw - hc)
  Yc=scale(1/n*(Xrc)%*%Ac0)
  SigmaB1hat=1/n*(Yc-mean(Yc))%*%t(Yc-mean(Yc))
  SigmaB2hat=1/p*(Yr-mean(Yr))%*%t(Yr-mean(Yr))
  eig2<-eigen(SigmaB1hat)
  lambdahatr = eig2$values[1:p]
  ind<-order(lambdahatr,decreasing=T)
  lambdahatr<-lambdahatr[ind]
  Qhat1<-eig2$vectors
  Qhat1=Qhat1[,ind]
  Qhatr<- Qhat1[, 1:p]
  ABr <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {ABr[, j] <- sqrt(lambdahatr[j]) * Qhatr[, j]}; ABr
  hBr <- diag(ABr %*% t(ABr))
  DBr <- diag(SigmaB1hat - hBr)
  eig3<-eigen(SigmaB2hat)
  lambdahatc = eig3$values[1:m]
  ind<-order(lambdahatc,decreasing=T)
  lambdahatc<-lambdahatc[ind]
  Qhat2<-eig3$vectors
  Qhat2=Qhat2[,ind]
  Qhatc<- Qhat2[1:p, 1:m]
  ABc <- matrix(0, nrow = p, ncol =m)
  for (j in 1:m) {ABc[, j] <- sqrt(lambdahatc[j]) * Qhatc[, j]}; ABc
  hBc<- diag(ABc%*% t(ABc))
  DBc<- diag(SigmaB2hat- hBc)
  return(
    list(ABr=ABr,ABc=ABc,DBr=DBr,DBc=DBc,
         SigmaB1hat=SigmaB1hat,SigmaB2hat=SigmaB2hat))}
