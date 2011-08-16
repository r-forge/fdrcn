\name{fdrcn-package}
\alias{fdrcn-package}
\alias{fdrcn}
\docType{package}
\title{
False discovery rate estimation for copy number data.
}
\description{
    This package estimates false discovery rates for multi-sample copy number data.  
}
\details{
\tabular{ll}{
Package: \tab fdrcn \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-09\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

    
}

\examples{

    data(gbmsnippet)
    res=fdrcn(gbmsnippet)


}

\author{

Nancy R. Zhang

Maintainer: Nancy R. Zhang (nzhang@stanford.edu)
}
\references{
Efron and Zhang, 2011, False discovery rates and copy number variation, Biometrika V. 98 pp 251-271.
}

\keyword{false discovery rates}
\keyword{scan statistic}
\keyword{change-point}
\keyword{DNA copy number}
