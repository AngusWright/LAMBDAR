\name{interp.2d}
\alias{interp.2d}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Interpolate 2D function
}
\description{
An interpolation function for 2D distributions
}
\usage{
interp.2d(x,y,obj)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{
numeric (n) array; x-indices of destination grid
}
  \item{y}{
numeric (n) array; y-indices of destination grid
}
  \item{obj}{
list; elements {x: the x indices of the original grid; y: the y indices of the original grid; zobj: the original 2D image matrix}
}
}
\details{
Function is used for interpolating images from one grid onto a new grid of arbitrary resolution. The function is not able to
perform rotations etc, but is particularly useful for modifying astronomical images to match each other in pixel space.
}
\value{
length(n,3) data.frame; x indicies, y indicies, and the interpolated image values at (x,y).

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}};
Aaron Robotham{ ICRAR}
 {\email{aaron.robotham@icrar.org}}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Generate the original grid/image
x.min<-y.min<--50
x.max<-y.max<-50
x<-seq(x.min,x.max)
y<-seq(y.min,y.max)
old.grid<-expand.grid(x,y)
#Make image matrix of r values
z<-matrix(sqrt(old.grid[,1]^2+old.grid[,2]^2),length(x),length(y))

#Generate the new grid to interpolate onto
x.new.min<-0
y.new.min<--15
x.new.max<-25
y.new.max<-15
x.new<-seq(x.new.min,x.new.max,by=0.1)
y.new<-seq(y.new.min,y.new.max,by=0.1)
new.grid<-expand.grid(x.new,y.new)

#Interpolate
z.new<-interp.2d(x=new.grid[,1],y=new.grid[,2],obj=list(x=x,y=y,z=z))

#Reconstruct image
z.new<-matrix(z.new[,3],length(x.new),length(y.new))

#Plot the two images
layout(cbind(1,2))
image(x=x,y=y,z=z,col=rainbow(1E1),xlab='X (pix)',
ylab='Y (pix)',xlim=c(x.min,x.max)*2/3,ylim=c(x.min,x.max)*2/3,asp=1,zlim=c(0,max(z)/2))
image(x=x.new,y=y.new,z=z.new,col=rainbow(1E1),zlim=c(0,max(z)/2),add=TRUE)
rect(x.new.min,y.new.min,x.new.max,y.new.max,col='gold',density=0,lwd=3)

#Using this implementation, the function can also be used
#to interpolate and sample from a function simultaneously.

#Generate sampling points (in the middle 2/3 of the image)
x.samp<-runif(350,min=x.min,max=x.max)*2/3
y.samp<-runif(350,min=y.min,max=y.max)*2/3
samp.grid<-cbind(x.samp,y.samp)

#Interpolate
z.samp<-interp.2d(x=samp.grid[,1],y=samp.grid[,2],obj=list(x=x,y=y,z=z))

#Reconstruct image
z.samp.im<-matrix(NA,length(x),length(y))
z.samp.im[cbind(x.samp-x.min+1,y.samp-y.min+1)]<-z.samp[,3]

#Plot the image
image(x=x,y=y,z=z,col=grey.colors(1E1),xlab='X (pix)',
ylab='Y (pix)',xlim=c(-50,50),ylim=c(-50,50),zlim=c(0,max(z)/2),asp=1)
image(x=x,y=y,z=z.samp.im,col=rainbow(1E1),zlim=c(0,max(z)/2),add=TRUE)


}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
\keyword{ Aperture }
