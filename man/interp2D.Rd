\name{interp2D}
\alias{interp2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate 2D function
}
\description{
An interpolation function for 2D distributions
}
\usage{
interp2D(x,y,obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric (n) array; x-indicies of destination grid
}
  \item{y}{
numeric (n) array; y-indicies of destination grid
}
  \item{obj}{
list; elements {x: the x indicies of the original grid; y: the y indicies of the original grid; zobj: the original 2D image matrix}
}
}
\details{
Function is used for interpolating images from one grid onto a new grid of arbitrary resolution. The function is not able to
perform rotations etc, but is particularly useful for modifying astronomical images to match each other in pixel space.
}
\value{
numeric (n,m) matrix; the interpolated image with axes n=length(x), m=length(y)

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angus H Wright, Aaron S.G. Robotham
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
