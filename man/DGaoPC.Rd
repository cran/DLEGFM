\name{DGaoPC}
\alias{DGaoPC}
\title{
Distributed Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.}
\usage{
DGaoPC(data,m,n1,K)
}
\arguments{
\item{data}{ The data is total data set}
\item{m}{The  m is the number of first layer principal component}
\item{n1}{The n1 is  the length of each data subset}
\item{K}{The K is the number of nodes}
}
\value{
 \item{AG1 }{estimation of load value}
 \item{AG2}{estimation of load value}
 \item{DG1 }{estimation of error term}
 \item{DG2 }{estimation of error term}
 \item{SigmahatG1}{estimation of covariance}
  \item{SigmahatG2}{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
DGaoPC(data=ISE,m=3,n1=107,K=5)
}

