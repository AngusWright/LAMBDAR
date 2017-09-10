\name{fill.grid}
\alias{fill.grid}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Fill a box with overlapping tiles 
}
\description{
Tile a box of radius 'rad' within the limits of 'ra' and 'dec', with overlaps of size 'buffer'. 
}
\usage{
fill.grid(ra, dec, buffer=0.1, rad=sqrt(3)/2, buffer.edge=TRUE, plot=TRUE, main="Block Coverage Map")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ra}{
numeric vector; The RA limits of the box, in degrees. 
}
  \item{dec}{
numeric vector; The DEC limits of the box, in degrees. 
}
  \item{buffer}{
numeric; The overlap between adjacent tiles, in degrees. 
}
  \item{rad}{ 
numeric; The radius of each tile, in degrees. 
}
  \item{buffer.edge}{ 
logical; do we need to buffer the limits of the box as well?  
}
  \item{plot}{ 
logical; do you want a plot of the box and tiles?  
}
  \item{main}{ 
string; the title for the plot.  
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
array of dimension 2; the central RA and DEC positions of each tile. 
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

ra<-c(120,140)
dec<-c(-3,2)
tile.cen<-fill.grid(ra,dec,plot=TRUE)

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }% __ONLY ONE__ keyword per line
