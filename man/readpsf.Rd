\name{readpsf}
\alias{readpsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read and Process PSF FITS Image
}
\description{
__This function is not designed to be called directly, but rather is an internal function to MeasureFluxes__

Reads an input fits image for a Point Spread Function (PSF) and performs calculations needed in MeasureFluxes.
However, when a non-zero gaussian FWHM (in arcsec) is provided, a gaussian psf will be generated. This is useful
if no PSF map is present.
}
\usage{
readpsf(env=NULL,filename,asperpix,apsize,confidence,normalize=TRUE,gauss_fwhm_as=0,outenv=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{filename}{
string; relative path to, and filename of, the PSF FITS file
}
  \item{asperpix}{
numeric; width of each image pixel in arcsec
}
  \item{apsize}{
numeric; SemiMajor axis of the largest aperture in catalogue
}
  \item{confidence}{
numeric; desired maximum confidence percentage used in generating/saving the PSF and PSF convolved stamps.
         PSF limits beyond the radius corresponding to this confidence level will be truncated.
         If confidence is 1 when a non-zero gaussian FWHM is present, confidence is instead forced to be
         1-(1E-16); i.e. the highest non-unity value possible.
}
  \item{normalize}{
logical; if TRUE, PSF will be normalised to a height of 1
}
  \item{gauss_fwhm_as}{
numeric; FWHM of PSF gaussian. If supplied and greater than 0, PSF will be generated using a gaussian and not from input file
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
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
\note{
 If both a PSF map and a gaussian FWHM is provided, the procedure will generate a gaussian PSF in preference to reading a PSF from file.
}

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
