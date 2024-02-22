\name{GaoPC}
\alias{GaoPC}
\title{
Loading Estimation for General Factor Model
}
\description{
This function estimates the load and residual terms based on the general factor model and calculates the estimated values.
}
\usage{
GaoPC(data,m)
}
\arguments{
\item{data}{The data is total data set}

\item{m}{The m is the number of principal component}
}
\value{
 \item{AG1}{estimation of load value}
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
GaoPC(data=ISE,m=3)
}

