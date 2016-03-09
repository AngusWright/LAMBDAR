\name{skyback.par}
\alias{skyback.par}
\alias{PlotSkyback}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Function to find local sky values
}
\description{
This function calculates local sky values using a combination source masking, sigma clipping of extreme pixels and iterating the median sky
values of weighted annuli that are located around the ra and dec of the object of interest. The function is designed to cope with quite complex
sky problems and is robust to outliers. It returns a single sky value and error for each requested ra and dec location.
}
\usage{
skyback.par(x_p, y_p, cutlo = 0, cuthi = 100, im_mask,
imm_mask=NULL, remmask = TRUE, radweight = 1, clipiters = 5,
PSFFWHMinPIX = 2/0.339, hardlo = 3, hardhi = 10,
probcut = 3, mpiopts = "")
PlotSkyback(id_g,x_p,y_p,stamplims=NULL,cutlo=0,
cuthi=100,im_mask,imm_mask=NULL,remmask=TRUE,radweight=1,
clipiters=5,PSFFWHMinPIX=2/0.339,hardlo=3,hardhi=10,
probcut=3,plotall=FALSE,path=NULL,
toFile=TRUE,res=120)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x_p}{
numeric (i) array; Vector of x-pixel values for the sky locations of interest (typically these will be the positions of objects).
}
  \item{y_p}{
numeric (i) array; Vector of y-pixel values for the sky locations of interest (typically these will be the positions of objects).
}
  \item{cutlo}{
numeric (i) array; Vector of low radius values (in pixels) where the code will start to calculate the sky annuli around the object.
Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius.
}
  \item{cuthi}{
numeric (i) array; Vector of low radius values (in pixels) where the code will stop calculating sky annuli around the object.
Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius.
}
  \item{im_mask}{
numeric (n,m) array; The data image stamps to use in sky estimation.
}
  \item{imm_mask}{
numeric (n,m) array; The sourcemask image stamps. This must be pixel matched to im_mask, with pixel values of 1 where the pixels are likely to
be sky, and 0 where they are likely to be sources. The SourceMask output of LAMBDAR can be used.
}
  \item{remmask}{
logical; should the maskim be used to remove known sources.
}
  \item{radweight}{
numeric; What radius power-law weighting should be used to bias the sky towards sky annuli nearer to the source. Larger values
weight the sky value more towards central values (default 1, so linear weighting).
}
  \item{clipiters}{
numeric; How many 3-sigma clips of the sky will be made.
}
  \item{PSFFWHMinPIX}{
numeric; The FWHM of the PSF for origim in pixels.
}
  \item{hardlo}{
numeric; The lower limit for the inner sky annuli, where the limit used is hardlo * PSFFWHMinPIX.
}
  \item{hardhi}{
numeric; The lower limit for the outer sky annuli, where the limit used is hardhi * PSFFWHMinPIX.
}
  \item{probcut}{
numeric; probability-level used for iterative clipping of sky-pixels in estimation.
}
  \item{mpiopts}{
list; options to be passed to the MPI foreach backend, if it is being used.
}
  \item{id_g}{
character (n) vector; id's of the sources.
}
  \item{stamplims}{
numeric (n,4) matrix; indicies of the image stamps in the full image space. If stamps are not cutup (im_mask is a matrix),
then the stamplims can be NULL.
}
  \item{plotall}{
logical; should all objects be plotted?
}
  \item{path}{
character; path to use for output, if images are being written to file.
}
  \item{toFile}{
logical; should figures be sent to file?
}
  \item{res}{
numeric; resolution in pixels per inch of the output PNGs.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a data frame with 5 columns containing the sky value (sky), sky error (skyerr) and the number of sky annuli that have
error bars encompassing the final sky (Nnearsky), the RMS of the sky pixels (skyRMS), and the Pearson's Normality Test p-value for the sky
pixels (skyRMSpval). Nnearsky can be at most 10, and will be 6-7 for a well behaved sky. skyRMSpval is a good diagnostic tool for determining
when the sky background is highly non-gaussian, usually indicating the presence of an unrecognised contaminant (e.g. stellar diffraction spike),
and can be used to screen out poorly derived/constrained fluxes.
}
\author{
Aaron Robotham{ ICRAR}
 {\email{aaron.robotham@icrar.org}};
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}

\examples{

#Load sample image
data(SDSS.sample)
#Load sample catalogue
data(ApCat.sample)
#Save sample image to File
write.fits(file="SampleImage.fits",SDSS.sample)
#Read Astrometry
astr<-read.astr("SampleImage.fits")
#Determine the pixel-locations of the aperture catalogue objects
xy<-ad2xy(ra=ApCat.sample$RAdeg,dec=ApCat.sample$DECdeg,astr_struc=astr)
#Get the Main Object's ID
main<-which.min(sqrt((xy[,1]-astr$NAXIS[1]/2)^2+(xy[,2]-astr$NAXIS[2]/2)^2))

#Perform a sky estimate for the Primary target
skyest<-skyback.par(xy[main,1],xy[main,2],im_mask=SDSS.sample$dat[[1]],remmask=FALSE)
#Show Sky Estimate
print(skyest[,c('sky','skyerr','skyRMS','Nnearsky')])
#Note that the sky is likely non-gaussian:
print(skyest[,'skyRMSpval'])

#Perform a sky estimates for purely gaussian noise
gauss<-array(rnorm(length(SDSS.sample$dat[[1]])),dim=dim(SDSS.sample$dat[[1]]))
skyest.gauss<-skyback.par(xy[,1],xy[,2], cutlo=rep(0,length(xy[,1])), cuthi=rep(100,length(xy[,1])),
im_mask=gauss,remmask=FALSE)
#But the p-value is still systematically low:
hist(skyest.gauss[,'skyRMSpval'],breaks=20,xlab='Sky RMS Normality p-value',col='red',density=20,angle=-45)
#This is becuase of the sigma-clipping:
skyest.gauss.noclip<-skyback.par(xy[,1],xy[,2], cutlo=rep(0,length(xy[,1])), cuthi=rep(100,length(xy[,1])),
,im_mask=gauss,remmask=FALSE,clipiters=0)
hist(skyest.gauss.noclip[,'skyRMSpval'],breaks=20,xlab='Sky RMS Normality p-value',add=TRUE,col='blue',density=10)
legend('topright',legend=c("With Clipping","Without Clipping"),col=c('red','blue'),lty=1,pch=0)

#Plot the Sky Estimate
skyest<-PlotSkyback(id_g="",x_p=xy[main,1],y_p=xy[main,2],im_mask=gauss,remmask=FALSE,toFile=FALSE)

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{Imaging}
\keyword{Background}
