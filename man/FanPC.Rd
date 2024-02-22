\name{FanPC}
\alias{FanPC}
\title{
Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.
}
\usage{
FanPC(data,m)
}
\arguments{
\item{data}{The data is total data set}
\item{m}{The m is the number of principal component}
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
FanPC(data=ISE,m=3)
}
