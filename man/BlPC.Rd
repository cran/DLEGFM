\name{BlPC}
\alias{BlPC}
\title{
Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.}
\usage{
BlPC(data,m)
}
\arguments{
\item{data}{The data is total data set}
\item{m}{The  m is the number of first layer principal component}
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
BlPC(data=ISE,m=3)
}


