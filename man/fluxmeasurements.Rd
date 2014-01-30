\name{fluxmeasurements}
\alias{fluxmeasurements}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Flux Measurements
}
\description{
Given an environment containing all the requisite parameters, the function output a FITS table containing measured
fluxes. Setup to this procedure is performed by a host of functions, the main of which is 'MeasureFluxes'. Generally speaking,
this function is not designed to be called directly.
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
  \item{sfafluxout}{
numeric (m) vector; object fluxes
}
  \item{sfaerrout}{
numeric (m) vector; object flux errors
}
  \item{sfafluxwout}{
numeric (m) vector; weighted object fluxes
}
  \item{sfaerrwout}{
numeric (m) vector; weighted object flux errors
}
  \item{dfafluxout}{
numeric (m) vector; deblended object fluxes
}
  \item{dfaerrout}{
numeric (m) vector; deblended object flux errors
}
  \item{dfafluxwout}{
numeric (m) vector; deblended weighted object fluxes
}
  \item{dfaerrwout}{
numeric (m) vector; deblended weighted object flux errors
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
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
