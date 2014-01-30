\name{MeasureFluxes}
\alias{MeasureFluxes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure Fluxes
}
\description{
Measure galactic/object fluxes from an arbitrary FITS image, when provided an input catalogue of object locations and aperture parameters.
Detailed descriptions of the many functions called by this procedure can be found in their relevant help files.
The procedure has many high-level options available; all options are designated an input parameter file.
}
\usage{
MeasureFluxes(parfile, quiet = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parfile}{
Parameter file path and name
}
  \item{quiet}{
logical; if TRUE all output to screen is supressed
}
  \item{...}{
Any arguments that are intended to be passed to the 'createparfile' routine
}
}
%%\details{
%%}
\value{
Measure Fluxes returns a list, contining some of the results which are sent to file.
  \item{SFAflux}{
  The flux of each aperture, without consideration of deblending
}
  \item{SFAerror}{
  The error on the un-deblended flux
}
  \item{SFAfluxWeight}{
  The weight associated with the un-deblended flux
}
  \item{SFAerrorWeight}{
  The weight associated with the un-deblended flux error
}
  \item{DFAflux}{
  The final deblended flux of each aperture
}
  \item{DFAerror}{
  The error associated with the final deblended flux for each aperture
}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
createparfile
}
\examples{
# To generate a default parameter file for use:
MeasureFluxes("--makepar")

# To generate a parameter file for use, with customised values:
##  i.e. Changing the input Catalogue name
MeasureFluxes("--makepar", Catalogue="MyCat.csv")

##  i.e. Changing the input Catalogue name, and changing number of processors
MeasureFluxes("--makepar", Catalogue="MyCat.csv", nProcessors=10)

# To execute the program using the parameter file "Lambdar.par":
MeasureFluxes("Lambdar.par")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
