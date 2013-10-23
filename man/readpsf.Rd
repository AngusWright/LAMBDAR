\name{readpsf}
\alias{readpsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and Process PSF FITS Image
}
\description{
Reads an input fits image for a Point Spread Function (PSF) and performs calculations needed in MeasureFluxes.
}
\usage{
readpsf(filename, asperpix, stampsize_as, normalize = FALSE, gauss_fwhm_as = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
string; relative path to, and filename of, the PSF FITS file
}
  \item{asperpix}{
numeric; width of each image pixel in arcsec
}
  \item{stampsize_as}{
numeric; desired width of PSF stamp in arcsec
}
  \item{normalize}{
logical; if TRUE, PSF will be normalised to a height of 1
}
  \item{gauss_fwhm_as}{
numeric; FWHM of PSF gaussian. If supplied and greater than 0, PSF will be generated using a gaussian and not from input file
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric (n,n) array; image containing the read/generated PSF. The pixel length of each dimension of the array is floor(stampsize_as/asperpix). The function also creates
some global variables, the most useful of which is:
\item{sumpsf}{ numeric; the sum of pixels that make the psf function }
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
