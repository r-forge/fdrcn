\name{Kplot}
\alias{Kplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Draws the dotplot using fdr matrix computed by fdrcn.
}
\description{
    This function draws a dot plot using fdr rates computed by fdrcn.  A dotplot places a dot (or any other simbol) at locations in the position by sample matrix where fdr is lower than alpha.
}
\usage{
    dotplot(fdr,th=.1,add=0,col=1,pch="-",...){
}
%- maybe also 'usage' for other objects documented here.

\arguments{
   \item{fdr}{matrix of local fdrs, 1- (tdr matrix returned by fdrcn)}
   \item{th}{threshold for placing a dot.}
   \item{add}{should we add to an existing dot plot?}
   \item{col}{color to be used.}
   \item{pch}{the character symbol to be used.}

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

%% ~Make other sections like Warning with \section{Warning 
}{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

    data(gbmsnippet)
    res=fdrcn(gbmsnippet)
    dotplot(1-gbmsnippet$tdrGains,th=0.05, add=FALSE, col="blue")
    dotplot(1-gbmsnippet$tdrLosses,th=0.05, add=TRUE, col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in 
the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
