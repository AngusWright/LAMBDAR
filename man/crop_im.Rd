\name{crop_im}
\alias{crop_im}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Crop FITS Image
}
\description{
Crop a fits image to a sub-region determined by the input radius, RA, and Dec.
}
\usage{
crop_im(ra0 = -999, dec0 = -999, pathroot = "./",
inpim = NA, cutrad = 1,fitsoutname = NA, extn=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ra0}{
numeric; RA of focus of crop (deg). If -999 is provided, the crop will centre around the central RA pixel of the provided FITS image (ie. NAXIS[1]/2).
}
  \item{dec0}{
numeric; Dec of focus of crop (deg). If -999 is provided, the crop will centre around the central DEC pixel of the provided FITS image (ie. NAXIS[2]/2).
}
  \item{pathroot}{
string; Path to input FITS image. Default is "./", i.e. the current R working directory.
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
  \item{extn}{
numeric; FITS extension to crop from
}
}
\details{
The inputs ra0, dec0, inpim, cutrad, and fitsoutname can all be vectors of values if more than one crop is desired.
Note that each cropped image _must_ have its own individual fitsoutname - the function will not create extended FITS files.
For example, if you wish to crop one image in many places: you can provide a single inpim, but vectors of `n' ra0s, dec0s, cutrads,
and fitsoutnames. If you wish to crop the same location out of many images, you can provide a single ra0, dec0, and cutrad, but
vectors of `n' inpims and fitsoutnames.
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
\examples{

#Load a Sample Image and Catalogue
data("SDSS.sample",envir=environment())
data("ApCat.sample",envir=environment())
#Write Sample Image to File
write.fits(file="SampleImage.fits",SDSS.sample)

#Crop a 40" (radius) cutout from the centre of the sample image:
crop_im(inpim="SampleImage.fits",cutrad=40/3600,fitsoutname="SampleImage_cropped.fits")

#Read the cropped sample image
SDSS.sample.cropped<-read.fits("SampleImage_cropped.fits")

#Image the Uncropped and Cropped Images
layout(cbind(1,2))
#Full Image
image(x=1:length(SDSS.sample$dat[[1]][,1]),y=1:length(SDSS.sample$dat[[1]][1,]),z=SDSS.sample$dat[[1]],
col=grey.colors(1E3),zlim=quantile(SDSS.sample$dat[[1]],c(0,0.999)),asp=1,useRaster=TRUE,
xlab='X (pix)',ylab="Y (pix)")
#Cropped Image
image(x=1:length(SDSS.sample.cropped$dat[[1]][,1]),y=1:length(SDSS.sample.cropped$dat[[1]][1,]),
z=SDSS.sample.cropped$dat[[1]],col=grey.colors(1E3),zlim=quantile(SDSS.sample.cropped$dat[[1]],c(0,0.99)),asp=1,
useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)")

#The cropped image has been modified to reflect the new astrometry:
#Get Object Positions in the cropped and uncropped images
xy<-ad2xy(ApCat.sample$RAdeg,ApCat.sample$DECdeg,read.astr("SampleImage.fits"))
xy.cropped<-ad2xy(ApCat.sample$RAdeg,ApCat.sample$DECdeg,read.astr("SampleImage_cropped.fits"))
aspp<-abs(read.astr("SampleImage.fits")$CD[1,1])*3600
#Plot Full Image w/ Object positions
image(x=1:length(SDSS.sample$dat[[1]][,1]),y=1:length(SDSS.sample$dat[[1]][1,]),z=SDSS.sample$dat[[1]],
col=grey.colors(1E3),zlim=quantile(SDSS.sample$dat[[1]],c(0,0.999)),asp=1,useRaster=TRUE,
xlab='X (pix)',ylab="Y (pix)")
points(xy,pch=3,col='red')
#Plot Cropped Image w/ Object positions
image(x=1:length(SDSS.sample.cropped$dat[[1]][,1]),y=1:length(SDSS.sample.cropped$dat[[1]][1,]),
z=SDSS.sample.cropped$dat[[1]],col=grey.colors(1E3),zlim=quantile(SDSS.sample.cropped$dat[[1]],c(0,0.99)),asp=1,
useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)")
points(xy.cropped,pch=3,col='red')


}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
\keyword{ Output }
