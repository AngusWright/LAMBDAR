\name{get.fwhm}
\alias{get.fwhm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Get PSF FWHM in pixels
}
\description{
Returns the pixel Full Width Half Max of an input image matrix
}
\usage{
get.fwhm(zdist)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{zdist}{
numeric (n,m) array; image for which FWHM is desired
}
}
\details{
The function determines the full width at half max of the input image, by searching outwards from the image maxima.
}
\value{
numeric; value of the FWHM in pixels
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
\examples{

  #Make Gaussian data: stdev 5 pix
  sd<-5
  x<-matrix(-10:10,21,21,byrow=TRUE)
  zdist<-matrix(exp(-1*(t(x)^2/(2*sd^2)+x^2/(2*sd^2))),21,21)
  print(paste(2*sqrt(2*log(2))*sd,"~",get.fwhm(zdist)))

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
\keyword{ PSF }% __ONLY ONE__ keyword per line
