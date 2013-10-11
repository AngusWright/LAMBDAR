\name{readparfile}
\alias{readparfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read Parameter File
}
\description{
Reads the Parameter file required for the execution of the MeasureFluxes procedure
}
\usage{
readparfile(parfile, starttime = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parfile}{
string; relative path to, and filename of, the parameter file
}
  \item{starttime}{
numeric; CPU start time of the calling procedure, in seconds. Used for output of elapsed times; if not supplied, current CPU time is used.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
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
MeasureFluxes
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
