\name{DGulPC}
\alias{DGulPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distributed Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.}
\usage{
DGulPC(data,m,n1,K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{The data is total data set}
\item{ m}{The  m is the number of first layer principal component}
\item{ n1}{ The n1 is  the length of each data subset}
\item{ K}{The K is the number of nodes}
}
\value{
 \item{AU1 }{estimation of load value}
 \item{AU2 }{estimation of load value}
 \item{DU3 }{estimation of error term}
 \item{S1hat }{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{DGulPC(data=ISE,m=3,n1=107,K=5)}

