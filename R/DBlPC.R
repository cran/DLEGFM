#' @param data is total data set
#' @param m is the number of principal component
#' @param n1 is  the length of each data subset
#' @param K is the number of nodes
#' @return ABr,ABc,DBr,DBc,SigmaB1hat,SigmaB2hat
#' @export
#' @examples
#' DBlPC(data=ISE,m=3,n1=107,K=5)
DBlPC=function(data,m,n1,K){
  Sigmahat=list()
  ABr=list()
  ABc=list()
  DBr=list()
  DBc=list()
  Ar0=list()
  Ac0=list()
  Dr=list()
  Dc=list()
  Yr=list()
  Yc=list()
  Sigmahatw=list()
  SigmaB1hat=list()
  SigmaB2hat=list()
  for (i in 1:K) {
    n=nrow(data)
    p=ncol(data)
    L=matrix(rep(0,K*n1),ncol=n1)
    R=matrix(0,n1,n)
    L[i,]=sample(1:n,n1,replace=FALSE)
    r=matrix(c(1:n1,L[i,]),ncol=n1,byrow=T)
    R[t(r)]=1
    X1=R%*%as.matrix(data)
    X=scale(X1)
    Sigmahat[[i]]=1/(n)*(X-mean(X))%*%t(X-mean(X))#=Ar%*%t(Ar)Ac%*%t(Ac)+cov(E)
    eig0<-eigen(Sigmahat[[i]])
    lambdar = eig0$values[1:p]
    ind<-order(lambdar,decreasing=T)
    lambdar<-lambdar[ind]
    Q1<-eig0$vectors
    Q1=Q1[,ind]
    Qr<- Q1[, 1:p]
    Ar01 <- matrix(0, nrow = n1, ncol = p)
    for (j in 1:p) {Ar01[, j] <- sqrt(lambdar[j]) * Qr[, j]}
    Ar0[[i]]=Ar01
    hr <- diag(Ar0[[i]] %*% t(Ar0[[i]]))
    Dr[[i]] <- diag(Sigmahat[[i]] - hr)
    Yr[[i]]=scale(1/p*t(X)%*%Ar0[[i]])
    Sigmahatw[[i]]=1/(n)*t(X-mean(X))%*%(X-mean(X))
    eig1<-eigen(Sigmahatw[[i]])
    lambdac = eig1$values[1:m]
    ind<-order(lambdac,decreasing=T)
    lambdahatc<-lambdac[ind]
    Q2<-eig1$vectors
    Q2=Q2[,ind]
    Qc<- Q2[1:p, 1:m]
    Ac01 <- matrix(0, nrow = p, ncol =m)
    for (j in 1:m) {Ac01[, j] <- sqrt(lambdac[j]) * Qc[, j]}
    Ac0[[i]]=Ac01
    hc <- diag(Ac0[[i]] %*% t(Ac0[[i]]))
    Dc[[i]] <- diag(Sigmahatw[[i]] - hc)
    Yc[[i]]=scale(1/n*(X)%*%Ac0[[i]])
    SigmaB1hat[[i]]=1/n*(Yc[[i]]-mean(Yc[[i]]))%*%t(Yc[[i]]-mean(Yc[[i]]))
    SigmaB2hat[[i]]=1/p*(Yr[[i]]-mean(Yr[[i]]))%*%t(Yr[[i]]-mean(Yr[[i]]))
    eig2<-eigen(SigmaB1hat[[i]])
    lambdahatr = eig2$values[1:p]
    ind<-order(lambdahatr,decreasing=T)
    lambdahatr<-lambdahatr[ind]
    Qhat1<-eig2$vectors
    Qhat1=Qhat1[,ind]
    Qhatr<- Qhat1[, 1:p]
    ABr1 <- matrix(0, nrow = n/K, ncol = p)
    for (j in 1:p) {ABr1[, j] <- sqrt(lambdahatr[j]) * Qhatr[, j]}
    ABr[[i]]=ABr1
    hBr <- diag(ABr[[i]] %*% t(ABr[[i]]))
    DBr[[i]] <- diag(SigmaB1hat[[i]] - hBr)
    eig3<-eigen(SigmaB2hat[[i]])
    lambdahatc = eig3$values[1:m]
    ind<-order(lambdahatc,decreasing=T)
    lambdahatc<-lambdahatc[ind]
    Qhat2<-eig3$vectors
    Qhat2=Qhat2[,ind]
    Qhatc<- Qhat2[1:p, 1:m]
    ABc1 <- matrix(0, nrow = p, ncol =m)
    for (j in 1:m) {ABc1[, j] <- sqrt(lambdahatc[j]) * Qhatc[, j]}
    ABc[[i]]=ABc1
    hBc<- diag(ABc[[i]]%*% t(ABc[[i]]))
    DBc[[i]]<- diag(SigmaB2hat[[i]]- hBc)}
  return(
    list(ABr=ABr,ABc=ABc,DBr=DBr,DBc=DBc,
         SigmaB1hat=SigmaB1hat,SigmaB2hat=SigmaB2hat))}
