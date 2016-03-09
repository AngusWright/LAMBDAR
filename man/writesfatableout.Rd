\name{write.flux.measurements.table}
\alias{write.flux.measurements.table}
\alias{write.deblend.fraction.table}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Write measure.fluxes Results
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Output a CSV table, to file, containing the measure.fluxes results
}
\usage{
write.flux.measurements.table(env=NULL, filename)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the enovironment containing the parameters and results of measure.fluxes.
}
  \item{filename}{
string; path to, and filename of, output filename
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
NULL; Output is to file
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
\link{measure.fluxes}, \link{flux.measurements}
}
\examples{

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Output }
