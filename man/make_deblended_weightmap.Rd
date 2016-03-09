\name{make_deblended_weightmap}
\alias{make_deblended_weightmap}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Deblend Photometric Objects
}
\description{
__This function is not designed to be called directly, but rather is an internal function to MeasureFluxes__

Deblends objects which are overlapping in a photometric image. Procedure is called by fluxmeasurements, part of the MeasureFluxes routine.
The procedure takes individual weighted object apertures, and determines what proportion of the blended flux should belong to the object in question.
It performs this analysis on all apertures from the input catalogue.
}
\usage{
make_deblended_weightmap(outenv=parent.env(environment()),env=NULL,subs=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in `MeasureFluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
  \item{subs}{
numeric length(k) vector; the subset of apertures to calculate deblending for. Default NULL is treated as all (i.e. 1:ngal).
}
}
\details{
The most relevant parameters to the execution of this function are the weighted & convolved aperture stamps (wsfa), and the
weighted & convolved full aperture mask (wfa):

- wsfa: numeric length (m) list of (n,n) arrays; a length (m) list of (n,n) weighted & convolved object aperture stamps
-  wfa: numeric (p,q) array; full mask of all weighted & convolved object apertures, placed in their appropriate locations, coadded

The deblended weightmaps are determined by what fraction of each aperture contributes to the full mask model at each pixel
in the objects stamp. I.e. Deblended weightmap @ (i,j)th pixel of an aperture stamp = wsfa[i,j]/wfa[n,m], where (n,m) is the
corresponding location of (i,j) in the full mask.
}
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
\link{MeasureFluxes}, \link{fluxmeasurements}
}
%%\examples{

%%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Deblend }
\keyword{ Aperture }% __ONLY ONE__ keyword per line
