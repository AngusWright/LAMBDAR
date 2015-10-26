\name{readimage}
\alias{readimage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Input Images
}
\description{
Read the 3 FITS images (data, error map, mask map) that are input to the procedure MeasureFluxes. The procedure opens the files that are listed in the
MeasureFluxes parameter file, and generate the relevant default images in the event that one is not specified.
}
\usage{
readimage(outenv=parent.env(environment()),
quiet=FALSE, showtime=FALSE, env=NULL)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{quiet}{
logical; if TRUE, all output is suppressed
}
  \item{showtime}{
logical; if TRUE, execution times are printed along with output
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
%- maybe also 'usage' for other objects documented here.
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
NULL; the procedure outputs images into the 'outenv' environment, for use by 'MeasureFluxes'. These images/objects are:
  \item{im}{
  2D array containing the fits image data
  }
  \item{imm}{
  2D array containing the mask for the data image
  }
  \item{ime}{
  2D array containing the error map for the data image
  }
  \item{hdr_str}{
  List containing the FITS header information for the data image
  }
  \item{hdr_mask}{
  List containing the FITS header information for the mask image
  }
  \item{hdr_err}{
  List containing the FITS header information for the error map image
  }
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
\link{MeasureFluxes}, \link{readparfile}, \link{createparfile}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
