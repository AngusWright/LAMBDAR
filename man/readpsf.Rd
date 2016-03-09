\name{read.psf}
\alias{read.psf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Read and Process PSF FITS Image
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Reads an input fits image for a Point Spread Function (PSF) and performs calculations needed in measure.fluxes.
However, when a non-zero Gaussian FWHM (in arcsec) is provided, a Gaussian PSF will be generated. This is useful
if no PSF map is present.
}
\usage{
read.psf(outenv = parent.env(environment()), filename,
arcsec.per.pix, apsize, confidence, normalize = TRUE,
gauss.fwhm.arcsec = 0, env = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in `measure.fluxes'
}
  \item{filename}{
string; relative path to, and filename of, the PSF FITS file
}
  \item{arcsec.per.pix}{
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
  \item{gauss.fwhm.arcsec}{
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
numeric (n,n) array; image containing the read/generated PSF. The pixel length of each dimension of the array is floor(stampsize_as/arcsec.per.pix). The function also creates
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
measure.fluxes, flux.measurements
}
%\examples{
%
%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ PSF }
\keyword{ Imaging }% __ONLY ONE__ keyword per line
