\name{Kplot}
\alias{Kplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Plots the K profile for output from fdrcn.
}
\description{
    This function plots the K profile (sum of tdr across samples) that is computed by fdrcn.  It also 
    plots the significance line (specified by the user parameter alpha) for max(K) computed from block
    bootstrap.
}
\usage{
    Kplot(KGains, KLosses, maxKGainsBoot, maxKLossesBoot, alpha=0.95)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
   \item{KGains}{K profile (sum of tdr across samples) for 
       gains (see fdrcn)}
   \item{KLosses}{K profile (sum of tdr across samples) for 
       losses (see fdrcn)}
   \item{maxKGainsBoot}{max(KGains) for each of the Bootstrap 
       samples (see fdrcn)}
   \item{maxKLossesBoot}{max(KLosses) for each of the Bootstrap 
       samples (see fdrcn)}
   \item{alpha}{the threshold for significance of max(K)}

}


\references{
    Efron and Zhang, 2011, False discovery rates and copy 
    number variation, Biometrika V. 98 pp 251-271.
}
\author{
Nancy R. Zhang
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

    data(gbmsnippet)
    res=fdrcn(gbmsnippet)
    Kplot(res$KGains,res$KLosses, res$maxKGainsBoot,res$maxKLossesBoot,0.99)

}
% Add one or more standard keywords, see file 'KEYWORDS' in 
the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
