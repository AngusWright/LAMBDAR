\name{fluxmeasurements}
\alias{fluxmeasurements}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Flux Measurements
}
\description{
Given image arrays and input object parameters, output a FITS table containing the fluxes of each object. Setup to this procedure is performed by the wrapping main function 
MeasureFluxes.
}
\usage{
fluxmeasurements(id_g, ra_g, dec_g, x_g, y_g, theta_g, a_g, b_g, im, hdr, ime, imm, pathroot, pathout, prfmap, lambda, field, conf, fluxcorr, fluxweight, beamarea_in, asperpix, Jybm, sourcemask, nopsf, filtcontam, contams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id_g}{
vector of (m) strings; The catalogue IDs of the objects to be measured
}
  \item{ra_g}{
numeric (m) vector; RAs of the objects to be measured in degrees
}
  \item{dec_g}{
numeric (m) vector; Decs of the objects to be measured in degrees
}
  \item{x_g}{
numeric (m) vector; x locations of the objects in the image-plane (in pixel number)
}
  \item{y_g}{
numeric (m) vector; y locations of the objects in the image-plane (in pixel number)
}
  \item{theta_g}{
numeric (m) vector; angles that the objects' elliptical apertures make with the y-axis (degrees; N0E90)
}
  \item{a_g}{
numeric (m) vector; lengths of the major axis' of the apertures in arcsec
}
  \item{b_g}{
numeric (m) vector; lengths of the minor axis' of the apertures in arcsec
}
  \item{im}{
numeric (n,n) array; input image
}
  \item{hdr}{
list; contains all keywords and values from the image FITS header 
}
  \item{ime}{
numeric (n,n) array; input image error map
}
  \item{imm}{
numeric (n,n) array; input image mask map
}
  \item{pathroot}{
string; relative path from current location to data/working directory
}
  \item{pathout}{
string; relative path to output directory
}
  \item{prfmap}{
string; filename of input prf/psf FITS image
}
  \item{conf}{
numeric; confusion noise in Janksy
}
  \item{fluxcorr}{
numeric; flux correction coefficient
}
  \item{fluxweight}{
numeric (m) vector; flux weighting applied to each object
}
  \item{beamarea_in}{
numeric; beam area of image
}
  \item{asperpix}{
numeric; width of each image pixel in arcsec
}
  \item{Jybm}{
logical; if TRUE, image pixels are measure in Jy per beam, as opposed to Jy per Pixel
}
  \item{sourcemask}{
logical; if TRUE, source mask FITS image will be output to file 
}
  \item{nopsf}{
logical; if TRUE, all objects will be convolved with a gaussian instead of an input PSF
}
  \item{filtcontam}{
logical; if TRUE, objects identified in the input catalogue as being 'contaminants' will be removed during object flux measurement
}
  \item{contams}{
numeric (m) vector; contaminant flags - if object 'i' is a contaminant, contams[i]=1; otherwise contams[i]=0
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
NULL; function creates global variables which may be used for further analysis after return, but otherwise all output is to file. 
The most relevant global variables are those defined last:
  \item{sfafluxout}{
numeric (m) vector; object fluxes
}
  \item{sfaerrout}{
numeric (m) vector; object flux errors
}
  \item{sfafluxwout}{
numeric (m) vector; weighted object fluxes
}
  \item{sfaerrwout}{
numeric (m) vector; weighted object flux errors
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
MeasureFluxes
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
