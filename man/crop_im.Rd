\name{crop_im}
\alias{crop_im}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crop FITS Image
}
\description{
Crop a fits image to a sub-region determined by the input radius, RA, and Dec.
}
\usage{
crop_im(ra0 = -999, dec0 = -999, pathroot = "", inpim = NA, cutrad = 1, fitsoutname = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ra0}{
numeric; RA of focus of crop (deg). If -999 is provided, the crop will centre around the central RA pixel of the provided FITS image (ie. NAXIS[1]/2).
}
  \item{dec0}{
numeric; Dec of focus of crop (deg). If -999 is provided, the crop will centre around the central DEC pixel of the provided FITS image (ie. NAXIS[2]/2).
}
  \item{pathroot}{
string; Path to input FITS image
}
  \item{inpim}{
string; Filename of input FITS image
}
  \item{cutrad}{
numeric; Radius of cropped region (deg). Default value is 1 degree.
}
  \item{fitsoutname}{
string; Filename for output cropped image
}
}
\details{
The inputs ra0, dec0, inpim, cutrad, and fitsoutname can all be vectors of values if more than one crop is desired.
Note that each cropped image _must_ have its own individual fitsoutname - the function will not create extended FITS files.
For example, if you wish to crop one image in many places: you can provide a single inpim, but vectors of 'n' ra0s, dec0s, cutrads,
and fitsoutnames. If you wish to crop the same location out of many images, you can provide a single ra0, dec0, and cutrad, but
vectors of 'n' inpims and fitsoutnames.
}
\value{
NULL; The cropped image(s) are output to file.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
%%\note{
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{MeasureFluxes}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
