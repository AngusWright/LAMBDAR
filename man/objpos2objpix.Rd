\name{objpos2objpix}
\alias{objpos2objpix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pixel Rounding
}
\description{
Converts floating point (precise) pixel values into integer pixel values.
}
\usage{
objpos2objpix(x_g, y_g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_g}{
numeric (n) vector; precise x pixel values for (n) objects
}
  \item{y_g}{
numeric (n) vector; precise y pixel value for (n) objectss
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric (2,n) array; row-bound vectors containing the integer pixel values for x and y
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
MeasureFluxes
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
