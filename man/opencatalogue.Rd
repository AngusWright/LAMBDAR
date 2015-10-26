\name{opencatalogue}
\alias{opencatalogue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Open and Read Catalogue
}
\description{
__This function is not designed to be called directly, but rather is an internal function to MeasureFluxes__

Open and read an input catalogue for the MeasureFluxes procedure.
}
\usage{
opencatalogue(outenv = parent.env(environment()), env = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
\details{
  The Catalogue is expected to contain a series of columns that will vary depending on the input parameters selected, and
  whose exact labels can be specified in the parameter file. For all runs of LAMBDAR, it is expected that the user catalogue
  specifies columns containing:
  Object IDs: Can be strings or numeric
  Object RAs: Expected to be in Decimal Degrees
  Object Decs: Expected to be in Decimal Degrees

  In addition to these, if the user has not specified that all objects should be forced to be point sources then the following
  additional columns must also be specified:
  Object Aperture SemiMajor Axis: Expected to be in arcseconds
  Object Aperture SemiMinor Axis: Expected to be in arcseconds
  Object Aperture orientation: Expected to be in decimal degrees, and in N90E0. However can be specified in N0E90 if the AngularOffset parameter
  is set to 1 in the parameter file.

  If the user has asked for Contaminants to be deblended and removed from the final catalogue (by specifying the RemoveContams parameter as 1)
  then it is expected that there is a boolean 1/0 column provided, containing flags for which objects are contaminants [1], or main objects [0].

  Finally, the user may optionally provide a column containing object weightings in conjunction with any other series of parameters. The weightings
  can be specified in terms of flux, magnitude, or scale, as determined by the users specification of the 'FluxWgtType' parameter.
}
\value{
NULL; The procedure assigns variables into the parameter environment 'outenv'.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
