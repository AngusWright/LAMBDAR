\name{make_deblended_weightmap}
\alias{make_deblended_weightmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deblend Photometeric Objects
}
\description{
Deblends objects which are overlapping in a photometric image. Procedure is called by fluxmeasurements, part of the MeasureFluxes routine.
The procedure takes individual weighted object apertures, and determines what proportion of the blended flux should belong to the object in question.
It performs this analysis on all apertures from the input catalogue.
}
\usage{
make_deblended_weightmap(wsfa, wfa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wsfa}{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) weighted object aperture stamps
}
  \item{wfa}{
numeric (p,q) array; full mask of all weighted object apertures, placed in their approriate locations, coadded
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) weighted, deblended, object aperture stamps
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
MeasureFluxes, fluxmeasurements
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
