\name{make.catalogue.apertures}
\alias{make.catalogue.apertures}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Make Aperture Stamps
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Create a length (m) list of (n,n) arrays, containing object apertures given input aperture parameters
}
\usage{
make.catalogue.apertures(outenv = parent.env(environment()), env = NULL)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in `measure.fluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
%- maybe also `usage' for other objects documented here.
\details{
Aperture parameters, which are used in the creation of the aperture stamps, are stored in the parameter environment `env'.
The returned list of aperture stamps is sent to `outenv'. If not supplied, `outenv' is set to be the same environment as `env'.
}
\value{
numeric [[m]][n,n] list; aperture stamps with pixel range [0,1]
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
\link{measure.fluxes}, \link{flux.measurements}
}
%%\examples{

%%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Aperture }
