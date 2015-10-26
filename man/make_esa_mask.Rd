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
make_esa_mask(outenv = parent.env(environment()),
env = NULL, ObsParm, padGals, col.corr = 0, confuse = FALSE)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
  \item{ObsParm}{
List containing Observation Parameters needed in creating simulated galaxies; Exposure Time (sec; exp=),
Collecting Area (meters squared; area=), Effective wavelength of filter (Å; lamEff=), and Effective width of filter (Å; WEff=).
For example, SDSS r-band has parameters ObsParm=list(exp=53.9,area=4.91,lamEff=1111.2,Weff=6165).
}
  \item{padGals}{
logical; if TRUE then the simulated image will have galaxies added beyond those in the input catalogue, such that
the galactic number counts follow the power law log(N)=0.38*m_r-4.78.
}
  \item{col.corr}{
numeric; The colour correction need to convert SDSS r-band magnitude to magnitude in the band being analysed.
This is needed when padding galaxies because our power law is defined in the r-band. Naturally, this is a 0th
order correction, and should be noted when determining the accuracy of simulations that are well beyond the
SDSS r-band.
}
  \item{confuse}{
logical; do you want to confuse the simulation such that the source density is greater than 1 object per FWHM^2
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
