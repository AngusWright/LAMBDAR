\name{make_sfa_mask}
\alias{make_sfa_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Weighted and Convolved Apertures
}
\description{
Convolve the input aperture stamps with the PSF, and weight the resultant apertures by the input fluxweightin.
}
\usage{
make_sfa_mask(sa_mask, fluxweightin = NULL, nopsf = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sa_mask}{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) aperture stamps
}
  \item{fluxweightin}{
numeric vector; can be either a length (m) vector of aperture weights, or a single numeric applied to all apertures as a scale-factor.
}
  \item{nopsf}{
logical; if TRUE, the apertures are convolved with a gaussian instead of an input PSF
}
}
\details{
Many of the parameters required by this function are setup in the global environment by the function MeasureFluxes and its subsidiaries.
}
\value{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) convolved, weighted, aperture stamps
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
