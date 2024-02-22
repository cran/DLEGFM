\name{DFanPC}
\alias{DFanPC}
\title{
 Distributed Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.
}
\usage{
DFanPC(data,m,n1,K)
}
\arguments{
\item{data}{The data is total data set}
\item{m}{The m is the number of principal component}
\item{n1}{The n1 is  the length of each data subset}
\item{K}{The K is the number of nodes}
}
\value{
 \item{AF}{estimation of load value}
 \item{DF }{estimation of error term}
  \item{SigmahatF}{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
DFanPC(data=ISE,m=3,n1=107,K=5)
}
