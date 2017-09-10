\name{radec.to.xy}
\alias{radec.to.xy}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
RA/Dec to Pixel
}
\description{
Maps a provided array of RAs and Decs to pixel locations X and Y using the TAN Gnomonic projection,
when provided with calibrating central RA/Dec, X/Y, and scales.
}
\usage{
radec.to.xy(ra, dec, ra0 = 0, dec0 = 0,
x0 = 0, y0 = 0, xscale = 1, yscale = 1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ra}{
numeric (n) vector; the (n) RAs to be converted (deg)
}
  \item{dec}{
numeric (n) vector; the (n) Decs to be converted (deg)
}
  \item{ra0}{
numeric; the RA at the calibrating pixel (deg)
}
  \item{dec0}{
numeric; the Dec at the calibrating pixel (deg)
}
  \item{x0}{
numeric; x-coordinate of the calibrating pixel
}
  \item{y0}{
numeric; x-coordinate of the calibrating pixel
}
  \item{xscale}{
numeric; RA-interval between successive x pixels (deg)
}
  \item{yscale}{
numeric; Dec-interval between successive y pixels (deg)
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric (n,2) array; two column bound vectors containing the desired x pixel locations [,'X'] and y pixel locations [,'Y']
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Aaron Robotham{ ICRAR}
 {\email{aaron.robotham@icrar.org}}
}
\note{
  This function is called by the ad.to.xy function, which extracts relevant information from an astrometry data.frame
  internally. It is therefore often more convenient to call ad.to.xy rather than to call radec.to.xy directly.
}
\note{
  This function is designed to work with the TAN Gnomonic projection. However, due to degeneracy the procedure
  can also be used with a orthographic SIN projection *without rotation*.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

#Load LAMBDAR
library(LAMBDAR)

#Load sample image
data(SDSS.sample)
#Load sample catalogue
data(ApCat.sample)
#Save sample image to File
write.fits(file="SampleImage.fits",SDSS.sample)
#Read Astrometry
astr<-read.astrometry("SampleImage.fits")
#Get the relevent parameters from the Astrometry Header
### NB: This snippet is within the ad.to.xy function, which is a
###     wrapper to the radec.to.xy function
ra0<-astr$CRVAL[1]
dec0<-astr$CRVAL[2]
x0<-astr$CRPIX[1]
y0<-astr$CRPIX[2]
xscale<-astr$CD[1,1]
yscale<-astr$CD[2,2]
#Determine the pixel-locations of the aperture catalogue objects
xy<-radec.to.xy(ra=ApCat.sample$RAdeg,dec=ApCat.sample$DECdeg, ra0=ra0,
dec0=dec0, x0=x0, y0=y0, xscale=xscale, yscale=yscale)
#Print sample
print(head(xy))

#Show the positions on the image
#Plot the Image
image(x=1:length(SDSS.sample$dat[[1]][,1]),
y=1:length(SDSS.sample$dat[[1]][1,]),
z=SDSS.sample$dat[[1]],col=grey.colors(1E3),
zlim=quantile(SDSS.sample$dat[[1]],c(0,0.995)),useRaster=TRUE,
xlab="X (pix)",ylab="Y (pix)",asp=1)
#Overlay the points
points(xy,pch=3,col='red')

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Astrometry }
