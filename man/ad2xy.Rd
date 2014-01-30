\name{ad2xy}
\alias{ad2xy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RA/Dec to Pixel via Header Astrometry
}
\description{
Converts RA/DEC (degrees) to X/Y (pixels) position using the TAN Gnomonic projection system.
}
\usage{
ad2xy(ra, dec, astr_struc, diagnostic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ra}{
numeric (n) vector; RA values
}
  \item{dec}{
numeric (n) vector; Dec values
}
  \item{astr_struc}{
list; Astrometry Keywords and values returned from read of a fits header. See "read.astr" function
}
  \item{diagnostic}{
logical; if TRUE, diagnostic information will be output as type 'message'
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric (n,2) array; two numeric (n) vectors 'X' and 'Y', column-bound
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR} \\
 {\email{angus.wright@icrar.org}}
}
\note{
 While this routine is designed for use with the TAN Gnomonic projection,
 an orthographic SIN projection *without rotation* is degenerate and therefore
 also valid for use
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{radec2xy}, \link{read.astr}, \link{MeasureFluxes}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
