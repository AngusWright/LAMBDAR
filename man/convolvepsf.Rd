\name{convolvepsf}
\alias{convolvepsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
2D Convolution
}
\description{
Perform the 2D convolution of two 2D arrays of equal dimension, using Fast Fourier Transforms (fft).
This function is designed to be used for convolutions of Point Spread Functions (PSFs) with object
apertures in photometric measurements. As such, the positional information of the first array (or
more accurately the phase information of the first array's fft) is discarded.
}
\usage{
convolvepsf(arr1, arr2, normalise = TRUE, mag2 = TRUE, zapdig = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{arr1}{
(n,n) dimensional array; Usually will contain the PSF function. Phase information of this array is discarded in convolution
}
  \item{arr2}{
(n,n) dimensional array; Phase information of this array is retained in convolution. Magnitude information is optionally included/excluded
}
  \item{normalise}{
logical; if TRUE, the normalised Fourier Transform will be calculated
}
  \item{nomag2}{
logical; if TRUE, then the magnitude information of arr2 is removed prior to FFT via normalisation
}
  \item{zapdig}{
integer; the precision to be used in a call to 'zapsmall'. Used for removing (unwanted) small variations in near-zero areas
of the array return by the FFT. If zapdig is NULL, getOption("digits") is used.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
(n,n) dimensional array containing the convolved & 'zapped' image.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR} \\
 {\email{angus.wright@icrar.org}}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{MeasureFluxes}, \link{zapsmall}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
