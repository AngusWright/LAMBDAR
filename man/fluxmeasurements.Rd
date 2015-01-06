\name{fluxmeasurements}
\alias{fluxmeasurements}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Flux Measurements
}
\description{
__This Function is not designed to be called directly; it is a subroutine designed to be called directly from MeasureFluxes__

Given an environment containing all the requisite parameters, the function output a FITS table containing measured
photometric fluxes measured from an input FITS file, and using input aperture parameters.
Setup to this procedure is performed by a host of functions, the main of which is 'MeasureFluxes'.

As with MeasureFluxes, this function will either return NULL or a list containing the measured fluxes, determined
by the user in the parameter file.
}
\usage{
fluxmeasurements(env=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
environment; The environment where all parameters are stored. See 'readparfile' for details.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  NULL; or a list containing the following:
  \item{sfafluxout}{
numeric (m) vector; object fluxes
}
  \item{sfaerrout}{
numeric (m) vector; object flux errors
}
  \item{dfafluxout}{
numeric (m) vector; deblended object fluxes
}
  \item{dfaerrout}{
numeric (m) vector; deblended object flux errors
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
\link{MeasureFluxes}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
