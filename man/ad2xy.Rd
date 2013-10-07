\name{ad2xy}
\alias{ad2xy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RA/Dec to Pixel via Header Astrometry
}
\description{
Convert RA,Dec locations into X,Y (pixel numbers) by reading key astrometry values from an input list. 
}
\usage{
ad2xy(ra, dec, astr_struc, diagnostic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ra}{
numeric (n) vector; RA of objects
}
  \item{dec}{
numeric (n) vector; Dec values
}
  \item{astr_struc}{
list of Astrometry Keywords and values returned from read of a fits header
}
  \item{diagnostic}{
logical; if TRUE, diagnostic information will be output as type 'message'
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric (n,2) array; two numeric (n) vectors 'X' and 'Y', column-bound
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angus H Wright{ ICRAR} \\
 {\email{angus.wright@icrar.org}}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
read.astr, MeasureFluxes
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
