\name{read.astr}
\alias{read.astr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Astrometry
}
\description{
Read the header of a provided FITS file, and return the astrometry associated with that file.
}
\usage{
read.astr(fitsname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitsname}{
string; relative path and filename of the FITS file to be analysed.
}
}
\details{
This procedure is designed to read all common astrometry keywords. Nonetheless, the source should be checked to ensure that your fits file keywords are included.
}
\value{
list; contains all keywords needed for determining locations in the fits image using astrometry.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
MeasureFluxes, fluxmeasurements
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
