\name{fdrcn}
\alias{fdrcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
    False discovery rate estimation for copy number data.
}
\description{
    This function computes false discovery rates for multi-sample copy number data.  It also finds positions that are enriched for gains (or losses) across profiles.
}
\usage{
fdrcn(x,win=11,scaleSamples=TRUE,nulltype=0,J=1,B=100,L=NA,alpha=.05,plots=TRUE,th=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \code{n} by \code{m} data matrix, for \code{n} locations and \code{m} samples.
  }
  \item{win}{
    window size for smoothing the data.
  }
  \item{scaleSamples}{
    Should data be scaled to have the same inter-quartile range across samples?
  }
  \item{nulltype}{
    0 for separate computation of FDRs for both gains and losses, 1 for only losses, 2 for only gains.
  }
  \item{J}{
    How many iterations should be used for estimating position wise signal proportion?
  }
  \item{B}{
    Number of bootstrap samples for computing the null quantiles of max(K).
  }
  \item{L}{
    Window-size used for block bootstrap.
  }
  \item{alpha}{
    Plotting parameter: the threshold for significance of max(K).
  }
  \item{plots}{
    Should plots be drawn?
  }
  \item{th}{
    Plotting parameter: the FDR threshold for the dotplot.
  }
}

\details{
    This is the main function of the \code{fdrcn} package.  It computes local false discovery rates for each position in each sample according to Efron and Zhang (2011).  Then, it finds
    the genome positions that are enriched for either gains or losses across samples.
}


\value{
 return(list(tdrGains=res2$Tdr,tdrLosses=res1$Tdr, KGains = K2, KLosses = K1, maxKGainsBoot =Kmax2, maxKLossesBoot = Kmax1))


   \item{tdrGains}{local true discovery rate for gains}
   \item{tdrLosses}{local true discovery rate for losses}
   \item{KGains}{K profile (sum of tdr across samples) for gains}
   \item{KLosses}{K profile (sum of tdr across samples) for losses}
   \item{maxKGainsBoot}{max(KGains) for each of the Bootstrap samples}
   \item{maxKLossesBoot}{max(KLosses) for each of the Bootstrap samples}
}
\references{
    Efron and Zhang, 2011, False discovery rates and copy number variation, Biometrika V. 98 pp 251-271.
}
\author{
Nancy R. Zhang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

    data(gbmsnippet)
    res=fdrcn(gbmsnippet)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
