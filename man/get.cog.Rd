\name{get.cog}
\alias{get.cog}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Generate Curve of Growth
}
\description{
When provided with an image matrix, generate the 2D curve of growth.
}
\usage{
get.cog(zdist, centre=NULL,sample=NULL,proj=NULL,SNR=FALSE,poly.degree=4,flexible=TRUE,na.rm=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{zdist}{
numeric (n,m) matrix; image
}
  \item{centre}{
numeric (length 2) vector; pixel location to centre COG on. If NULL, CoG is centred on the
brightest pixel. Can be approximate if flexible=TRUE; the function will find the max closest to this point. 
}
  \item{sample}{
numeric; Number of radial points to sample COG at. Dictates output length.
}
  \item{proj}{
numeric (length 2) vector; Parameters for calculating de-projected COG. Must be length 2
with first value Axel Ratio (b/a) and second value position angle (PA).
}
  \item{SNR}{
logical; Do you want the CoG returned as a SNR rather than straight flux? 
}
  \item{poly.degree}{
numeric; degree of the polynomial used for fits and derivation of slope/concavity. 
}
  \item{flexible}{
logical; Is the centre only approxmately correct? This allows the centre to be moved slightly. 
}
  \item{na.rm}{
logical; Remove bad data prior to cog measurement?  
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
List of length 2;
    $x: numeric vector; Radius values
    $y: numeric vector; COG values
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

#Load LAMBDAR
library(LAMBDAR)

# nxn Dimension of the image
n<-101
#Gaussian (PSF) signal with Gaussian noise
#Pixel Indicies
y<-matrix(1:n,n,n)
x<-t(y)
#PSF Parameters
x.sigma=1.5
y.sigma=2.5
flux<-300
#Make Point Source
z<-exp(-(((x-n/2)^2/(2*x.sigma^2))+((y-n/2)^2/(2*y.sigma^2))))
z<-z/sum(z)*flux

##Random Noise
sky<-matrix(rnorm(n^2,sd=0.1), ncol=n)
#Get Curve of Growth
cog<-get.cog(z+sky, c(n/2,n/2))
#Plot
dev.new(height=5,width=8)
layout(cbind(1,2,3))
magplot(cog$avg$x,cog$avg$y, main="Curve of growth", xlab="Radius (pix)",
ylab="Flux within Radius", type='l')
abline(h=flux,lty=2,col='darkgreen')

#Gauusian (PSF) signal with Gaussian noise and negative sky
##Random Noise & negative sky
sky<-matrix(rnorm(n^2,mean=-0.02,sd=0.1), ncol=n)
#Get Curve of Growth
cog<-get.cog(z+sky, c(n/2,n/2))
#Plot
magplot(cog$avg$x,cog$avg$y, main="Curve of growth; -ve sky", xlab="Radius (pix)",
ylab="Flux within Radius", type='l')
abline(h=flux,lty=2,col='darkgreen')

#Gauusian (PSF) signal with Gaussian noise and positive sky
##Random Noise & negative sky
sky<-matrix(rnorm(n^2,mean=+0.02,sd=0.1), ncol=n)
#Get Curve of Growth
cog<-get.cog(z+sky, c(n/2,n/2))
#Plot
magplot(cog$avg$x,cog$avg$y, main="Curve of growth; +ve sky", xlab="Radius (pix)",
ylab="Flux within Radius", type='l')
abline(h=flux,lty=2,col='darkgreen')

#Exponential signal with Gaussian noise
#Exponential parameters:
bn<-1.678
R.Eff.pix<-10
pa<-45
axel.ratio<-0.5
flux<-300
#Pixel Indicies
xy<-expand.grid(1:(2*n),1:(2*n))
xy[,1]<-xy[,1]-n
xy[,2]<-xy[,2]-n
#Generate Exponential Profile
xy<-rotate.data.2d(xy[,1],xy[,2],pa)
xy[,2]<-xy[,2]*axel.ratio
exp.val<-exp(-bn*(sqrt((xy[,1])^2+(xy[,2])^2)/R.Eff.pix-1))
z<-matrix(exp.val,n*2,n*2)
z<-z[seq(ceiling(n*0.5),floor(n*1.5)),seq(ceiling(n*0.5),floor(n*1.5))]
z<-z/sum(z)*flux

##Random Noise
sky<-matrix(rnorm(n^2,sd=0.1), ncol=n)
#Get Curve of Growth
cog<-get.cog(z+sky, c(n/2,n/2))
#Plot
layout(cbind(1,2))
magplot(cog$avg$x,cog$avg$y, main="Curve of growth; Exponential", xlab="Radius (pix)",
ylab="Flux within Radius", type='l')
abline(h=flux,lty=2,col='darkgreen')

#Curve of Growth can be measured using 'Deprojected' radii,
#if apertures are robust:
#Get Deprojected Curve of Growth
deproj.cog<-get.cog(z+sky, c(n/2,n/2),proj=c(axel.ratio,pa))
magplot(cog$avg$x,cog$avg$y, main="Deprojected Curve of growth; Exponential",
xlab="Deprojected Radius (pix)", ylab="Flux within Dep. Radius",
type='l')
abline(h=flux,lty=2,col='darkgreen')


#CoGs can be used to measure half-light radii,
#which are useful for aperture diagnostics
#Plot Source Image
image(x=seq(ceiling(-n/2),floor(n/2)),y=seq(ceiling(-n/2),floor(n/2)),z=z,
xlab="X (pix)",ylab="Y (pix)",main="Source",col=grey.colors(1E3),asp=1)
#Get Half Light Radius
hlr<-min(cog$avg$x[which(cog$avg$y>=max(cog$avg$y,na.rm=TRUE)/2)])
deproj.hlr<-min(deproj.cog$avg$x[which(deproj.cog$avg$y>=max(deproj.cog$avg$y,na.rm=TRUE)/2)])
#Overplot Input Source Half Light Radius:
lines(ellipse(a=R.Eff.pix,e=1-axel.ratio,pa=pa),col='forestgreen',lty=1)
#Overplot Image half light radius:
lines(ellipse(a=hlr,e=0,pa=0),col='red',lty=3)
#Overplot Deprojected half light radius:
lines(ellipse(a=deproj.hlr,e=1-axel.ratio,pa=pa),col='gold',lty=2)

#Plot Source Image
image(x=seq(ceiling(-n/2),floor(n/2)),y=seq(ceiling(-n/2),floor(n/2)),
z=z+sky,xlab="X (pix)",ylab="Y (pix)",main="Source+Noise",
col=grey.colors(1E3),asp=1)
#Get Half Light Radius
hlr<-min(cog$avg$x[which(cog$avg$y>=max(cog$avg$y,na.rm=TRUE)/2)])
deproj.hlr<-min(deproj.cog$avg$x[which(deproj.cog$avg$y>=max(deproj.cog$avg$y,na.rm=TRUE)/2)])
#Overplot Input Source Half Light Radius:
lines(ellipse(a=R.Eff.pix,e=1-axel.ratio,pa=pa),col='forestgreen',lty=1)
#Overplot Image half light radius:
lines(ellipse(a=hlr,e=0,pa=0),col='red',lty=3)
#Overplot Deprojected half light radius:
lines(ellipse(a=deproj.hlr,e=1-axel.ratio,pa=pa),col='gold',lty=2)

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Photometry }
\keyword{ Imaging }% __ONLY ONE__ keyword per line
\keyword{ Background }
\keyword{ Aperture }
