\name{rancor.par}
\alias{rancor.par}
\alias{plotRanCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform LAMBDAR randoms/blanks correction
}
\description{
__This function is not designed to be called directly, but rather is an internal function to MeasureFluxes__

This function performs a parallelised random/blank aperture analysis on the image stamp for each aperture under consideration.
}
\usage{
rancor.par(im_mask,imm_mask,ap_mask,stamplims,
masklims,remask=TRUE,numIters=1E3,mpiopts="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im_mask}{
numeric [[m]][n,n] list; stamps of the input image
}
  \item{imm_mask}{
numeric [[m]][n,n] list; stamps of the input mask image
}
  \item{ap_mask}{
numeric [[m]][k,k] list; aperture stamps
}
  \item{stamplims}{
numeric [m,4] array; indicies of the aperture stamps in the im_mask space
}
  \item{masklims}{
numeric [m,4] array; indicies of the aperture stamps in the imm_mask space
}
  \item{remask}{
logical; if TRUE, Blanks will be measured. Otherwise, randoms will be measured
}
  \item{numIters}{
numeric; number of blanks/randoms
}
  \item{mpiopts}{
list; the MPI options to be passed to the MPI backend, if it is in use.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%randMean.mean=mean(flux,na.rm=T),randMed.mean=mean(fluxmed,na.rm=T),randMean.SD=sd(flux,na.rm=T),randMed.SD=sd(fluxmed,na.rm=T),nRand=length(which(is.finite(flux))))
Returns a data frame with 5 columns for randoms/blanks, containing the mean of the mean pixel values per aperture (___Mean.mean),
the mean of the median pixel values per aperture (___Med.mean), the standard deviation of the mean pixel values (___Mean.SD), the
standard deviation of the median pixel values (___Med.SD), and the number of randoms/blanks that had finite measurements within the aperture (n____).
}
\author{
Angus H Wright
}

\examples{
##No examples to run, perhaps in final release version of LAMBDAR
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sky}
\keyword{background}
