#' @param data is a total data set
#' @param m is the number of principal component
#' @param n1 is  the length of each data subset
#' @param K is the number of nodes
#' @return AG1,AG2,DG3,SGhat
#' @export
#' @examples
#' DGaoPC(data=ISE,m=3,n1=107,K=5)
DGaoPC=function(data,m,n1,K){
  SigmahatG1=list()
  SigmahatG2=list()
  AG1=list()
  AG2=list()
  DG1=list()
  DG2=list()
  DG3=list()
  F1hat=list()
  sGhat=list()
  for (i in 1:K) {
    n=nrow(data)
    pc=2
    L=matrix(rep(0,K*n1),ncol=n1)
    R=matrix(0,n1,n,ncol=n)
    L[i,]=sample(1:n,n1,replace=FALSE)
    r=matrix(c(1:n1,L[i,]),ncol=n1,byrow=T)
    R[t(r)]=1
    X1=R%*%as.matrix(data)
    X=scale(X1)
    SigmahatG1[[i]]<-cor(X)
    eig3<-eigen(SigmahatG1[[i]])
    lambda1hat =eig3$values[1:m]
    ind<-order(lambda1hat,decreasing=T)
    lambda1hat<-lambda1hat[ind]
    Q<- eig3$vectors
    Q<-Q[,ind]
    AG1[[i]]<-Q[,1:m]
    hG1 <- diag(AG1[[i]] %*% t(AG1[[i]]))
    DG1[[i]]<- diag(SigmahatG1[[i]] - hG1)
    p<-ncol(X)
    F1hat[[i]]=X%*%AG1[[i]]
    F1star<-F1hat[[i]]/sqrt(n)
    SigmahatG2[[i]]=cor(F1star)
    eig4<-eigen(SigmahatG2[[i]])
    lambda2hat = eig4$values[1:pc]
    ind<-order(lambda2hat,decreasing=T)
    lambda2hat<-lambda2hat[ind]
    Q<- eig4$vectors
    Q<-Q[,ind]
    AG2[[i]]<-Q[,1:pc]
    hG2 <- diag(AG2[[i]] %*% t(AG2[[i]]))
    DG2[[i]] <- diag(SigmahatG2[[i]] - hG2)
    Fhat=F1star%*%AG2[[i]]
    XGhat=Fhat%*%t(AG2[[i]])%*%t(AG1[[i]])
    sGhat[[i]]=cov(XGhat)
    hG3 <- diag(t(t(AG2[[i]])%*%t(AG1[[i]]))%*%(t(AG2[[i]])%*%t(AG1[[i]])))
    DG3[[i]] <- diag(sGhat[[i]] - hG3)}
  return(list(AG1=AG1,AG2=AG2,DG3=DG3,sGhat=sGhat))}
