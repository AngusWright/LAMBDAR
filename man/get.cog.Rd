\name{get.cog}
\alias{get.cog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Curve of Growth
}
\description{
When provided with an image matrix, generate the 2D curve of growth.
}
\usage{
get.cog(zdist, centre = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zdist}{
numeric (n,m) matrix; image
}
  \item{centre}{
numeric (length 2) vector; pixel location to centre COG on
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
Angus H Wright{ ICRAR} \\
 {\email{angus.wright@icrar.org}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#Load Library
library(LAMBDAR)

#Gauusian signal with Gaussian noise
#Create data
z<-matrix(0, ncol=100, nrow=100)
##Add Signal
Amplitude=10
x.sigma=1.5       #pix
y.sigma=2.5       #pix
centre=c(50,50)   #pix
for (i in  30:70) {
  for (j in 30:70) {
    z[i,j]=Amplitude*exp(-(((i-centre[1])^2/(2*x.sigma^2))+((j-centre[2])^2/(2*y.sigma^2))))
  }
}
##Random Noise
sky<-matrix(rnorm(1E4,sd=0.1), ncol=100)
#Get Curve of Growth
cog<-get.cog(z+sky, centre)
#Plot
plot(cog, main="Curve of growth", xlab="Radius (pix)", ylab="COG", type='l')

#Gauusian signal with Gaussian noise and negative sky
##Random Noise & negative sky
sky<-matrix(rnorm(1E4,mean=-0.02,sd=0.1), ncol=100)
#Get Curve of Growth
cog<-get.cog(z+sky, centre)
#Plot
plot(cog, main="Curve of growth; -ve sky", xlab="Radius (pix)", ylab="COG", type='l')

#Gauusian signal with Gaussian noise and positive sky
##Random Noise & negative sky
sky<-matrix(rnorm(1E4,mean=+0.02,sd=0.1), ncol=100)
#Get Curve of Growth
cog<-get.cog(z+sky, centre)
#Plot
plot(cog, main="Curve of growth; +ve sky", xlab="Radius (pix)", ylab="COG", type='l')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
