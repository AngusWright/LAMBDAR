\name{make_esa_mask}
\alias{make_esa_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Exponential Profile Stamps
}
\description{
__This function is not designed to be called directly, but rather is an internal function to CreateSim__

Create a length (m) list of (n,n) arrays, containing exponential profiles determined by the input aperture parameters.
The exponential profiles are constructed assuming that the input apertures are 2.5*Kron apertures, such as those
output by SExtractor. Using these parameters, the function then creates realistic exponential profiles using a
Monte Carlo integration of simulated photon counts on a per-object basis. In this way, the simulation correctly
models the central region of the exponential profile, while simultaneously generating realistic shot-noise variations.
Objects that are supplied with zero width (ie. point sources) are simulated as perfect gaussians.
}
\usage{
make_esa_mask(env=NULL,outenv=NULL)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
%- maybe also 'usage' for other objects documented here.
\value{
numeric [[m]][n,n] list; stamps containing exponential profiles of simulated galaxies
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
\link{MeasureFluxes}, \link{fluxmeasurements}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
