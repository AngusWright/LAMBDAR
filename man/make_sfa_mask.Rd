\name{make.convolved.apertures}
\alias{make.convolved.apertures}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Make Weighted and Convolved Apertures
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Convolve the input aperture stamps with the PSF, and weight the resultant apertures by the input flux.weightin.
}
\usage{
make.convolved.apertures(outenv = parent.env(environment()),
sa_mask, flux.weightin = NULL, immask = NULL, env = NULL, subs = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in `measure.fluxes'
}
  \item{sa_mask}{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) aperture stamps
}
  \item{flux.weightin}{
numeric vector; can be either a length (m) vector of aperture weights, or a single numeric applied to all apertures as a scale-factor.
}
  \item{immask}{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) mask stamps
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
  \item{subs}{
numeric length(k) vector; the subset of apertures to convolve and weight. Default NULL is treated as 1:length(m).
}
}
\value{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) convolved & weighted aperture stamps
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
\link{measure.fluxes}, \link{flux.measurements}, \link{convolve.psf}
}
%%\examples{

%%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Aperture }
\keyword{ Deblend }
\keyword{ Convolution }% __ONLY ONE__ keyword per line

