\name{skyback}
\alias{skyback}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find local sky values
}
\description{
This function calculates local sky values using a combination source masking, sigma clipping of extreme pixels and iterating the median sky values of weighted annuli that are located around the ra and dec of the object of interest. The function is designed to cope with quite complex sky problems and is robust to outliers. It returns a single sky value and error for each requested ra and dec location.
}
\usage{
skyback(ra, dec, cutlo = 0, cuthi = 100, origim, astrom, maskim, remmask = TRUE, radweight = 1, clipiters = 5, PSFFWHMinPIX = 2/0.339, hardlo = 2, hardhi = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ra}{
Vector of Right Ascensions for the sky locations of interest (typically these will be the positions of objects).
}
  \item{dec}{
Vector of Declinations for the sky locations of interest (typically these will be the positions of objects). Should be same length as ra.
}
  \item{cutlo}{
Vector of low radius values (in pixels) where the code will start to calculate the sky annuli around the object. Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius. Should be same length as ra.
}
  \item{cuthi}{
Vector of low radius values (in pixels) where the code will stop calculating sky annuli around the object. Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius. Should be same length as ra.
}
  \item{origim}{
The data image to use.
}
  \item{astrom}{
Astrometry structure for the image being used, as provided by the read.astr function in LAMBDAR.
}
  \item{maskim}{
The source maskim. This must be pixel matched to origim, with pixel values of 1 where the pixels are likely to be sky, and 0 where they are likely to be sources. The SourceMask output of LAMBDAR can be used.
}
  \item{remmask}{
Logical- should the maskim be used to remove known sources.
}
  \item{radweight}{
What radius power-law weighting should be used to bias the sky towards sky annuli nearer to the source. Larger values weight the sky value more towards central values (default 1, so linear weighting).
}
  \item{clipiters}{
How many 3-sigma clips of the sky will be made.
}
  \item{PSFFWHMinPIX}{
The FWHM of the PSF for origim in pixels.
}
  \item{hardlo}{
The lower limit for the inner sky annuli, where the limit used is hardlo * PSFFWHMinPIX.
}
  \item{hardhi}{
The lower limit for the outer sky annuli, where the limit used is hardhi * PSFFWHMinPIX.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a data frame with 3 columns containing the sky value (sky), sky error (skyerr) and the number of sky annuli that have error bars encompassing the final sky (Nnearsky), where this can be at most 10, and will be 6-7 for a well behaved sky.
}
\author{
Aaron Robotham
}

\examples{
##No examples to run, perhaps in final release version of LAMBDAR
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sky}
\keyword{background}
