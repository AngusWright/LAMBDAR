\name{ran.cor}
\alias{ran.cor}
\alias{plot.ran.cor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Perform LAMBDAR randoms/blanks estimation
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

This function performs a parallelised random/blank aperture analysis on the image stamp for each aperture under consideration.
}
\usage{
ran.cor(data.stamp,ap.stamp,mask.stamp=NULL,ap.stamp.lims=NULL,
data.stamp.lims=NULL,rem.mask=FALSE,numIters=1E2,mpi.opts="",
sigclip=3,nclip=0,cat.x=NULL,cat.y=NULL,rand.x=NULL,rand.y=NULL,
ran.main.mask.lim=0.99)
plot.ran.cor(data.stamp,ap.stamp,mask.stamp=NULL,ap.stamp.lims=NULL,
data.stamp.lims=NULL,mask.stamp.lims=NULL,toFile=FALSE,rem.mask=FALSE,
numIters=1E2,path="./",plot.sci=FALSE,contams=NULL,plot.all=FALSE,
sigclip=3,nclip=3,res=120,cat.id=NULL,cat.x=NULL,cat.y=NULL,
rand.x=NULL,rand.y=NULL,ran.main.mask.lim=0.99)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.stamp}{
numeric [[m]][n,n] list; stamps of the input image
}
  \item{ap.stamp}{
numeric [[m]][k,k] list; aperture stamps
}
  \item{mask.stamp}{
numeric [[m]][n,n] list; stamps of the input mask image
}
  \item{ap.stamp.lims}{
numeric [m,4] array; indicies of the aperture stamps in the data.stamp space
}
  \item{data.stamp.lims}{
numeric [m,4] array; indicies of the data stamps in data image space (used
for calibrating cat.x and cat.y positions).
}
  \item{rem.mask}{
logical; if TRUE, Blanks will be measured. Otherwise, randoms will be measured
}
  \item{numIters}{
numeric; number of blanks/randoms
}
  \item{mpi.opts}{
list; the MPI options to be passed to the MPI backend, if it is in use.
}
  \item{sigclip}{
numeric; sigma to use for clipping of randoms/blanks, to improve randoms/blanks stdev estimate. [Not currently implemented in .par file]
}
  \item{nclip}{
numeric; number of iterative sigma-clips to do when calculating the stdev of randoms/blanks, to improve randoms/blanks stdev estimate. [Not currently implemented in .par file]
}
  \item{cat.x}{
numeric [m] vector; x (pixel) coordinates of sources in the data image space
}
  \item{cat.y}{
numeric [m] vector; y (pixel) coordinates of sources in the data image space
}
  \item{toFile}{
logical; output images as PNG?
}
  \item{path}{
character; path to use for output, if images are being written to file.
}
  \item{plot.sci}{
logical; should all of the science objects be plotted?
}
  \item{contams}{
logical [m] vector; the flag for sources that are contaminants (TRUE) or science targets (FALSE); used when plot.sci == TRUE?
}
  \item{plot.all}{
logical; should all objects be plotted?
}
  \item{res}{
numeric; resolution (ppi) of output PNG
}
  \item{cat.id}{
character [m] vector; id's of sources, for naming output files
}
  \item{rand.x}{
numeric [m] vector; x (pixel) coordinates of randoms in the data image space. If NULL, x values are generated at random.
}
  \item{rand.y}{
numeric [m] vector; y (pixel) coordinates of randoms in the data image space. If NULL, y values are generated at random.
}
  \item{ran.main.mask.lim}{
numeric; integration limit in [0,1] for masking the main source in the image
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%return=data.frame(randMean.mean=mean(flux,na.rm=T),randMean.SD=sd(flux,na.rm=T),randMean.MAD=mad(flux,na.rm=T),nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
Returns a data frame with 7 columns for randoms/blanks, containing the mean of the mean pixel values per aperture (___Mean.mean),
the standard deviation of the mean pixel values (___Mean.SD), the median absolute deviation from median of the mean pixel values (___Mean.MAD),
the mean area (after masking) of the apertures used for randoms/blanks measurements (___Ap.mean), the standard deviation of the aperture areas used
for randoms/blanks measurements (___Ap.SD), the median absolute deviation from median of the aperture areas used for randoms/blanks measurements
(___Ap.MAD), and the number of randoms/blanks that had finite measurements within the aperture (n____).
}
\author{
Angus H. Wright { ICRAR}
 {\email{angus.wright@icrar.org}}
}

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
#Determine the pixel-locations of the aperture catalogue objects
xy<-ad.to.xy(ra=ApCat.sample$RAdeg,dec=ApCat.sample$DECdeg,astr.struc=astr)
#Get the Main Object's ID
main<-which.min(sqrt((xy[,1]-astr$NAXIS[1]/2)^2+(xy[,2]-astr$NAXIS[2]/2)^2))
#Make the Main Object's aperture
#Make aperture grid
grid<-expand.grid(1:length(SDSS.sample$dat[[1]][,1]),1:length(SDSS.sample$dat[[1]][1,]))
SDSS.pixel.res<-0.339
#Generate Values
aper<-generate.aperture(x=grid[,1],y=grid[,2],xstep=1,ystep=1,
xcen=xy[main,1],ycen=xy[main,2],
axrat=ApCat.sample$radminasec[main]/ApCat.sample$radmajasec[main],
axang=-1*ApCat.sample$rotN2E[main],majax=ApCat.sample$radmajasec[main]/SDSS.pixel.res,
resample.iterations=0)
#Convert return back into matrix
aper<-matrix(aper[,3],nrow=(length(SDSS.sample$dat[[1]][,1])))

#Perform a randoms estimate the Primary target
randoms<-ran.cor(data.stamp=SDSS.sample$dat[[1]],ap.stamp=aper,
cat.x=xy[main,1],cat.y=xy[main,2])
#Show randoms estimate
print(randoms[,c('randMean.mean','randMean.SD','randMean.MAD')])

#Perform a randoms estimate for purely gaussian noise
gauss<-array(rnorm(length(SDSS.sample$dat[[1]])),dim=dim(SDSS.sample$dat[[1]]))
randoms.gauss<-ran.cor(data.stamp=gauss,ap.stamp=aper,
cat.x=xy[main,1],cat.y=xy[main,2])
#Show randoms estimate
print(randoms.gauss[,c('randMean.mean','randMean.SD','randMean.MAD')])

#Plot the randoms estimate
randoms<-plot.ran.cor(data.stamp=gauss,ap.stamp=aper,
cat.x=xy[main,1],cat.y=xy[main,2])
}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{Imaging}
\keyword{Background}
