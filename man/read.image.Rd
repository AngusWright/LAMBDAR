\name{read.images}
\alias{read.images}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Read Input Images
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Read the FITS images (data, error map, mask map, and/or weight map) that are input to the procedure measure.fluxes.
The procedure opens the files that are listed in the measure.fluxes parameter file, and generate the relevant
default images in the event that one is not specified.
}
\usage{
read.images(outenv=parent.env(environment()),
quiet=FALSE, showtime=FALSE, env=NULL)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in `measure.fluxes'
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
%- maybe also `usage' for other objects documented here.
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
NULL; the procedure outputs images into the `outenv' environment, for use by `measure.fluxes'. These images/objects may include:
  \item{im}{
  2D array containing the fits image data
  }
  \item{imm}{
  2D array containing the mask for the data image
  }
  \item{ime}{
  2D array containing the error map for the data image
  }
  \item{imwt}{
  2D array containing the weight map for the data image
  }
  \item{data.hdr}{
  List containing the FITS header information for the data image
  }
  \item{mask.hdr}{
  List containing the FITS header information for the mask image
  }
  \item{error.hdr}{
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
\link{measure.fluxes}, \link{read.par.file}, \link{create.par.file}
}
%\examples{
%
%}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Imaging }
