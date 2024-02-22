\name{PPC}
\alias{PPC}
\title{Loading Estimation for General Factor Model}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.
}
\usage{
PPC(data,m)
}
\arguments{
\item{data}{The data is total data set}
\item{m}{The m is the number of principal component}
}
\value{
 \item{Apro}{estimation of load value}
 \item{Dpro }{estimation of error term}
  \item{Sigmahatpro}{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
PPC(data=ISE,m=3)
}
