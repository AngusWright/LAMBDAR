\name{sph2car}
\alias{sph2car}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spherical to Cartesian Conversion
}
\description{
Does what is says on the tin
}
\usage{
sph2car(long, lat, radius = 1, deg = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long}{
numeric; logitude values
}
  \item{lat}{
numeric; latitude values
}
  \item{radius}{
numeric; radius values
}
  \item{deg}{
logical; are the values in degree [TRUE] or radians [FALSE]
}
}
%\details{
%
%}
\value{
data.frame with cartesian coordinates x, y, z.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aaron S.G. Robotham
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
