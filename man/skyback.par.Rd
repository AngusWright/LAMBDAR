\name{skyback.par}
\alias{skyback.par}
\alias{PlotSkyback}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find local sky values
}
\description{
This function calculates local sky values using a combination source masking, sigma clipping of extreme pixels and iterating the median sky
values of weighted annuli that are located around the ra and dec of the object of interest. The function is designed to cope with quite complex
sky problems and is robust to outliers. It returns a single sky value and error for each requested ra and dec location.
}
\usage{
skyback.par(x_p, y_p, cutlo = 0, cuthi = 100, im_mask,
imm_mask, remmask = TRUE, radweight = 1, clipiters = 5,
PSFFWHMinPIX = 2/0.339, hardlo = 3, hardhi = 10,
probcut = 3, mpiopts = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x_p}{
numeric (i) array; Vector of x-pixel values for the sky locations of interest (typically these will be the positions of objects).
}
  \item{y_p}{
numeric (i) array; Vector of y-pixel values for the sky locations of interest (typically these will be the positions of objects).
}
  \item{cutlo}{
numeric (i) array; Vector of low radius values (in pixels) where the code will start to calculate the sky annuli around the object.
Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius.
}
  \item{cuthi}{
numeric (i) array; Vector of low radius values (in pixels) where the code will stop calculating sky annuli around the object.
Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius.
}
  \item{im_mask}{
numeric (n,m) array; The data image stamps to use in sky estimation.
}
  \item{imm_mask}{
numeric (n,m) array; The sourcemask image stamps. This must be pixel matched to im_mask, with pixel values of 1 where the pixels are likely to
be sky, and 0 where they are likely to be sources. The SourceMask output of LAMBDAR can be used.
}
  \item{remmask}{
logical; should the maskim be used to remove known sources.
}
  \item{radweight}{
numeric; What radius power-law weighting should be used to bias the sky towards sky annuli nearer to the source. Larger values
weight the sky value more towards central values (default 1, so linear weighting).
}
  \item{clipiters}{
numeric; How many 3-sigma clips of the sky will be made.
}
  \item{PSFFWHMinPIX}{
numeric; The FWHM of the PSF for origim in pixels.
}
  \item{hardlo}{
numeric; The lower limit for the inner sky annuli, where the limit used is hardlo * PSFFWHMinPIX.
}
  \item{hardhi}{
numeric; The lower limit for the outer sky annuli, where the limit used is hardhi * PSFFWHMinPIX.
}
  \item{probcut}{
numeric; pribability-level used for iterative clipping of sky-pixels in estimation.
}
  \item{mpiopts}{
list; options to be passed to the MPI foreach backend, if it is being used.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a data frame with 5 columns containing the sky value (sky), sky error (skyerr) and the number of sky annuli that have
error bars encompassing the final sky (Nnearsky), the RMS of the sky pixels (skyRMS), and the Pearson's Normality Test p-value for the sky
pixels (skyRMSpval). Nnearsky can be at most 10, and will be 6-7 for a well behaved sky. skyRMSpval is a good diagnostic tool for determining
when the sky background is highly non-gaussian, usually indicating the presence of an unrecognised contaminant (e.g. stellar diffraction spike),
and can be used to screen out poorly derived/constrained fluxes.
}
\author{
Aaron Robotham, Angus H Wright
}

\examples{
##No examples to run, perhaps in final release version of LAMBDAR
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sky}
\keyword{background}
