\name{xy2radec}
\alias{xy2radec}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Pixel to RA/Dec
}
\description{
Convert pixel locations to RA and Dec by using a calibrating pixel and known pixel widths,
in the TAN Gnomonic coordinate system.
}
\usage{
xy2radec(x, y, ra0 = 0, dec0 = 0,
x0 = 0, y0 = 0, xscale = 1, yscale = 1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{
numeric (n) vector; x-coordinates of pixels
}
  \item{y}{
numeric (n) vector; y-coordinates of pixels
}
  \item{ra0}{
numeric; RA of the calibrating pixel (deg)
}
  \item{dec0}{
numeric; Dec of the calibrating pixel (deg)
}
  \item{x0}{
numeric; x-coordinate of the calibrating pixel
}
  \item{y0}{
numeric; y-coordinate of the calibrating pixel
}
  \item{xscale}{
numeric; width of pixels in RA (deg)
}
  \item{yscale}{
numeric; height of pixels in Dec (deg)
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric (n,2) array; RA and Dec vectors (column bound) of the input x,y coordinates
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
MeasureFluxes
}
\examples{

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Astrometry }
