\name{make_data_mask}
\alias{make_data_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Data Stamps
}
\description{
__This function is not designed to be called directly, but rather is an internal function to MeasureFluxes__

Create a series of length (m) list of (n,n) arrays, containing the data cutouts around all input catalogue objects.
}
\usage{
make_data_mask(outenv=parent.env(environment()),env=NULL)
}
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
%- maybe also 'usage' for other objects documented here.
\details{
This function creates cutouts of the input image, input mask, and error map. This is done so that parallelisation is as efficient
as possible.
}
\value{
numeric [[m]][n,n] list; stamps of the input image
Further numeric [[m]][n,n] list may be generated, but are assigned to the output environment directly, and are not returned.
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
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
