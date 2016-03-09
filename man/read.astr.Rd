\name{read.astr}
\alias{read.astr}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Read Astrometry
}
\description{
Read the header of a provided FITS file, and return the astrometry associated with that file.
}
\usage{
read.astr(fitsname,hdu=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fitsname}{
string; relative path and filename of the FITS file to be analysed.
}
  \item{hdu}{
integer; fits extension number
}
}
\details{
This procedure is designed to read all common astrometry keywords. Nonetheless, the source should be checked to ensure that your fits file keywords are included.
}
\value{
list; contains all keywords needed for determining locations in the fits image using astrometry.
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
\examples{

#Load sample image
data(SDSS.sample)
#Save sample image to File
write.fits(file="SampleImage.fits",SDSS.sample)
#Show the Full FITS Header
hdr<-read.fitshdr("SampleImage.fits")
print(hdr)
#Read Astrometry
astr<-read.astr("SampleImage.fits")
#Show the Astrometry Parameters
print(astr)

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Astrometry }
