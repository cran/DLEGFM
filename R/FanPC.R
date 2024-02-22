#' @param data is a total data set
#' @param m is the number of principal component
#' @return AF,DF,SigmahatF
#' @export
#' @examples
#' FanPC(data=ISE,m=3)
FanPC=function(data,m){
  X=scale(data)
  n=nrow(X)
  SigmahatF=cor(X)
  eig<-eigen(SigmahatF)
  lambdahat = eig$values[1:m]
  ind<-order(lambdahat,decreasing=T)
  lambdahat<-lambdahat[ind]
  Q<- eig$vectors
  Q<-Q[,ind]
  AF<-Q[,1:m]
  hF <- diag(AF %*% t(AF))
  DF <- diag(SigmahatF - hF)
  return(
    list(AF=AF,DF=DF,SigmahatF=SigmahatF))}
