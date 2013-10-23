\name{opencatalogue}
\alias{opencatalogue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Open and Read Catalogue
}
\description{
Open and read an input catalogue for the MeasureFluxes procedure.
}
\usage{
opencatalogue(pathroot, catalogue, forcepointsources)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathroot}{
string; relative path from current directory to data/working directory
}
  \item{catalogue}{
string; filename of the catalogue
}
  \item{forcepointsources}{
logical; if TRUE, aperture parameters are not read from catalogue and standard gaussian values are assigned
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
NULL; The procedure itself returns nothing. However it sets global variables that are of import:
  \item{id_g}{
string (n) vector; catalogue IDs for the (n) objects
}
  \item{ra_g}{
numeric (n) vector; RAs for the (n) objects
}
  \item{dec_g}{
numeric (n) vector; Decs for the (n) objects
}
  \item{a_g}{
numeric (n) vector; Aperture major axis for the (n) objects (arcsec)
}
  \item{b_g}{
numeric (n) vector; Aperture minor for the (n) objects (arcsec)
}
  \item{theta_g}{
numeric (n) vector; Aperture major axis orientation for the (n) objects. Angle from y-axis (degrees; N0E90)
}
  \item{fluxweight}{
numeric (n) vector; weights for the (n) objects. Used in scaling and in deblending
}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
