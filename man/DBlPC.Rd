\name{DBlPC}
\alias{DBlPC}
\title{
Distributed Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.}
\usage{
DBlPC(data,m,n1,K)
}
\arguments{
\item{data}{ The data is total data set}
\item{m}{The  m is the number of first layer principal component}
\item{n1}{The n1 is  the length of each data subset}
\item{K}{The K is the number of nodes}
}
\value{
 \item{ABr }{estimation of load value}
 \item{ABc }{estimation of load value}
 \item{DBr }{estimation of error term}
 \item{DBc }{estimation of error term}
 \item{SigmaB1hat}{estimation of covariance}
  \item{SigmaB2hat}{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
DBlPC(data=ISE,m=3,n1=107,K=5)
}
