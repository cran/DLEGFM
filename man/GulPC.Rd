\name{GulPC}
\alias{GulPC}
\title{Loading Estimation for General Factor Model}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.}
\usage{
GulPC(data,m)
}
\arguments{
\item{data}{The data is total data set}
\item{m}{The  m is the number of first layer principal component}
}
\value{
 \item{AU1 }{estimation of load value}
 \item{AU2 }{estimation of load value}
 \item{DU3 }{estimation of error term}
 \item{Shat }{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
GulPC(data=ISE,m=3) }

