\name{ran.cor}
\alias{ran.cor}
\alias{plot.ran.cor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Perform LAMBDAR randoms/blanks estimation
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

This function performs a parallelised random/blank aperture analysis on the image stamp for each aperture under consideration.
}
\usage{
ran.cor(data.stamp,mask.stamp,ap.stamp,stamplims,
masklims,remask=TRUE,numIters=1E3,mpi.opts="",sigclip=3,nclip=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data.stamp}{
numeric [[m]][n,n] list; stamps of the input image
}
  \item{mask.stamp}{
numeric [[m]][n,n] list; stamps of the input mask image
}
  \item{ap.stamp}{
numeric [[m]][k,k] list; aperture stamps
}
  \item{stamplims}{
numeric [m,4] array; indicies of the aperture stamps in the data.stamp space
}
  \item{masklims}{
numeric [m,4] array; indicies of the aperture stamps in the mask.stamp space
}
  \item{remask}{
logical; if TRUE, Blanks will be measured. Otherwise, randoms will be measured
}
  \item{numIters}{
numeric; number of blanks/randoms
}
  \item{mpi.opts}{
list; the MPI options to be passed to the MPI backend, if it is in use.
}
  \item{sigclip}{
numeric; sigma to use for clipping of randoms/blanks, to improve randoms/blanks stdev estimate. [Not currently implemented in .par file]
}
  \item{nclip}{
numeric; number of iterative sigma-clips to do when calculating the stdev of randoms/blanks, to improve randoms/blanks stdev estimate. [Not currently implemented in .par file]
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%return=data.frame(randMean.mean=mean(flux,na.rm=T),randMean.SD=sd(flux,na.rm=T),randMean.MAD=mad(flux,na.rm=T),nRand=length(which(is.finite(flux))),randAp.mean=mean(sumap,na.rm=T),randAp.SD=sd(sumap,na.rm=T),randAp.MAD=mad(sumap,na.rm=T))
Returns a data frame with 7 columns for randoms/blanks, containing the mean of the mean pixel values per aperture (___Mean.mean),
the standard deviation of the mean pixel values (___Mean.SD), the median absolute deviation from median of the mean pixel values (___Mean.MAD),
the mean area (after masking) of the apertures used for randoms/blanks measurements (___Ap.mean), the standard deviation of the aperture areas used
for randoms/blanks measurements (___Ap.SD), the median absolute deviation from median of the aperture areas used for randoms/blanks measurements
(___Ap.MAD), and the number of randoms/blanks that had finite measurements within the aperture (n____).
}
\author{
Angus H. Wright { ICRAR}
 {\email{angus.wright@icrar.org}}
}

\examples{
##No examples to run, perhaps in final release version of LAMBDAR
}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{Imaging}
\keyword{Background}
