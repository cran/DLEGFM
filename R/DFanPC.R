#' @param data is a total data set
#' @param m is the number of principal component
#' @param n1 is  the length of each data subset
#' @param K is the number of nodes
#' @return ABF,DF,SigmahatF
#' @export
#' @examples
#' DFanPC(data=ISE,m=3,n1=107,K=5)
DFanPC=function(data,m,n1,K){
  SigmahatF=list()
  AF=list()
  DF=list()
  for (i in 1:K) {
    n=nrow(data)
    L=matrix(rep(0,K*n1),ncol=n1)
    R=matrix(0,n1,n)
    L[i,]=sample(1:n,n1,replace=FALSE)
    r=matrix(c(1:n1,L[i,]),ncol=n1,byrow=T)
    R[t(r)]=1
    X1=R%*%as.matrix(data)
    X=scale(X1)
    SigmahatF[[i]]<-cor(X)
    eig2<-eigen(SigmahatF[[i]])
    lambdahat = eig2$values[1:m]
    ind<-order(lambdahat,decreasing=T)
    lambdahat<-lambdahat[ind]
    Q<- eig2$vectors
    Q<-Q[,ind]
    AF[[i]]<-Q[,1:m]
    hF <- diag(AF[[i]] %*% t(AF[[i]]))
    DF[[i]]<- diag(SigmahatF[[i]] - hF) }
  return(
    list(AF=AF,DF=DF,SigmahatF=SigmahatF))}
