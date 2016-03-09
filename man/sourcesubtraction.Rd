\name{source.subtraction}
\alias{source.subtraction}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Source Subtraction from FITS image
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Perform a subtraction of sources and/or contaminants from a FITS image.
}
\usage{
source.subtraction(image, sfa_models, ap.lims.data.stamp,
fluxes, outputfile, outputheader, beamarea, inside.mask,
diagnostic = FALSE, verbose = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{image}{
numeric (j,k) array; the image to be used for subtraction
}
  \item{sfa_models}{
numeric length (m) list of (n,n) arrays; length (m) list of (n,n) stamps containing the apertures to be subtracted
}
  \item{ap.lims.data.stamp}{
numeric (m,4) array; limits of where the aperture stamps are to be placed, in image pixel space; [,1]=x_lo, [,2]=x_up, [,3]=y_lo, [,4]=y_up.
}
  \item{fluxes}{
numeric (m) vector; flux weights used in normalisation of apertures
}
  \item{outputfile}{
string; relative path to, and filename of, the output subtracted image
}
  \item{outputheader}{
list; the FITS header for the output image
}
  \item{beamarea}{
numeric; the beam area for the input image
}
  \item{inside.mask}{
logical (m) vector; if TRUE, then the object will be subtracted
}
  \item{diagnostic}{
logical; if TRUE, then perform diagnostic checks and outputs
}
  \item{verbose}{
logical; if TRUE, then perform verbose outputs
}
}
\details{
Procedure performs a subtraction of the provided sfa_models from `image'. The models are normalised so that the
integral of the ith model is fluxes[i]. This means that the source subtraction is only going to accurately remove all
object flux (leaving only noise) in cases where the object is perfectly modeled. As this is rarely the case, the
result of this source subtraction is generally the creation of `butterfly' patterns. This is an indication that objects
are not being modeled perfectly by the aperture-psf convolution.
}
\value{
NULL; all output is to file
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
\link{measure.fluxes}, \link{flux.measurements}
}
%\examples{
%
%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
