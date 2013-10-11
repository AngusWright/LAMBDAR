\name{radec2xy}
\alias{radec2xy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RA/Dec to Pixel
}
\description{
Maps a provided array of RAs and Decs to pixel locations X and Y, when provided with calibrating central RA/Dec, X/Y, and scales.
}
\usage{
radec2xy(ra, dec, ra0 = 0, dec0 = 0, x0 = 0, y0 = 0, xscale = 1, yscale = 1)
}
%- maybe also 'usage' for other objects documented here.
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
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
