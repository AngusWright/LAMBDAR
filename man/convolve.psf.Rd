\name{convolve.psf}
\alias{convolve.psf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
2D Convolution
}
\description{
Perform the 2D convolution of two 2D arrays of equal dimension, using Fast Fourier Transforms (fft).
This function is designed to be used for convolutions of Point Spread Functions (PSFs) with object
apertures in photometric measurements. As such, the positional information of the first array (or
more accurately the phase information of the first array's fft) is discarded.
}
\usage{
convolve.psf(arr1, arr2, normalise = TRUE, renorm.arr2 = FALSE,
zapdig = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{arr1}{
(n,n) dimensional array; Usually will contain the PSF function. Phase information of this array is discarded in convolution
}
  \item{arr2}{
(n,n) dimensional array; Phase information of this array is retained in convolution. Magnitude information is optionally included/excluded
}
  \item{normalise}{
logical; if TRUE, the normalised Fourier Transform will be calculated
}
  \item{renorm.arr2}{
logical; if TRUE, then the magnitude information of arr2 is removed prior to FFT via normalisation
}
  \item{zapdig}{
integer; the precision to be used in a call to `zapsmall'. Used for removing (unwanted) small variations in near-zero areas
of the array return by the FFT. If zapdig is NULL, getOption("digits") is used.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
(n,n) dimensional array containing the convolved & `zapped' image.
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
\link{measure.fluxes}, \link{zapsmall}
}
\examples{

#Load LAMBDAR
library(LAMBDAR)

#Generate a 2D Gaussian PSF
#Stamp size (in pix)
n<-100
#Sigma
sigma.pix<-4
#Centre
x0<-y0<-n/2
y<-matrix(1:n,n,n)
x<-t(y)
psf<-exp(-(((x - x0)^2/(2 * sigma.pix^2)) + ((y -
            y0)^2/(2 * sigma.pix^2))))

#Make a Delta Function
ap<-matrix(0,n,n)
ap[n/2,n/2]<-1

#Show the convolution
layout(cbind(1,2,3))
image(x=1:n,y=1:n,z=psf,col=grey.colors(1E3),useRaster=TRUE,xlab='X (pix)',
ylab="Y (pix)",main="PSF Function",asp=1)
image(x=1:n,y=1:n,z=ap,col=grey.colors(1E3),useRaster=TRUE,xlab='X (pix)',
ylab="Y (pix)",main="Aperture Function",asp=1)
image(x=1:n,y=1:n,z=convolve.psf(psf,ap),col=grey.colors(1E3),useRaster=TRUE,
xlab='X (pix)',ylab="Y (pix)",main="Convolved Aperture Function",asp=1)

#Make an Aperture Function
#Pixel Grid
ap.pix<-expand.grid(seq(1,n),seq(1,n))
#Generate Aperture Pixel values
ap.pix.vals <- generate.aperture(x=ap.pix[,1],y=ap.pix[,2],xstep=1,ystep=1,
xcen=n/2,ycen=n/2,majax=10,axrat=0.5,axang=45)
#Construct Aperture Matrix
ap<-matrix(ap.pix.vals[, 3], ncol = n)

#Show the convolution
layout(cbind(1,2,3))
image(x=1:n,y=1:n,z=psf,col=grey.colors(1E3),useRaster=TRUE,xlab='X (pix)',
ylab="Y (pix)",main="PSF Function",asp=1)
image(x=1:n,y=1:n,z=ap,col=grey.colors(1E3),useRaster=TRUE,xlab='X (pix)',
ylab="Y (pix)",main="Aperture Function",asp=1)
image(x=1:n,y=1:n,z=convolve.psf(psf,ap),col=grey.colors(1E3),useRaster=TRUE,
xlab='X (pix)',ylab="Y (pix)",main="Convolved Aperture Function",asp=1)

#Phase information in the PSF-matrix is discarded during convolution,
#meaning that the PSF can be arbitrarily placed on its stamp:
#Big PSF, located in the stamp corner
big.psf<-matrix(0,n*3,n*3)
big.psf[1:n,1:n]<-psf
#A Big Aperture
big.ap.pix<-expand.grid(seq(1,n*3),seq(1,n*3))
big.ap.pix.vals <- generate.aperture(x=big.ap.pix[,1],y=big.ap.pix[,2],xstep=1,
ystep=1,xcen=n*3/2,ycen=n*3/2,majax=30,axrat=0.5,axang=45)
big.ap<-matrix(big.ap.pix.vals[, 3], ncol = n*3)

#Show the convolution
layout(cbind(1,2,3))
image(x=1:(n*3),y=1:(n*3),z=big.psf,col=grey.colors(1E3),useRaster=TRUE,
xlab='X (pix)',ylab="Y (pix)",main="PSF Function",asp=1)
image(x=1:(n*3),y=1:(n*3),z=big.ap,col=grey.colors(1E3),useRaster=TRUE,xlab='X (pix)',
ylab="Y (pix)",main="Aperture Function",asp=1)
image(x=1:(n*3),y=1:(n*3),z=convolve.psf(big.psf,big.ap),col=grey.colors(1E3),
useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)",main="Convolved Aperture Function",asp=1)


}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
\keyword{ Aperture }
\keyword{ Convolution }% __ONLY ONE__ keyword per line
\keyword{ PSF }
