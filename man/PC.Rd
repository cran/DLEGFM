\name{PC}
\alias{PC}
\title{Loading Estimation for General Factor Model}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.
}
\usage{
PC(data,m)
}
\arguments{
\item{data}{The data is a highly correlated data set}

\item{m}{The m is the number of principal component}
}
\value{
 \item{Ahat}{estimation of load value}
 \item{Dhat }{estimation of error term}
  \item{Sigmahat}{estimation of covariance}
}
\author{
Guangbao Guo, Yaping Li
}
\examples{
PC(data=ISE,m=3)
}
