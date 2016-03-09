\name{read.par.file}
\alias{read.par.file}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Read Parameter File
}
\description{
__This function is not designed to be called directly, but rather is an internal function to measure.fluxes__

Reads the Parameter file required for the execution of the measure.fluxes procedure. This parameter file can be
generated using the `create.par.file' routine, or by calling create.par.file through `measure.fluxes'.
}
\usage{
read.par.file(par.file=NA, start.time=NA,
quiet=FALSE, env=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{par.file}{
string; relative path to, and filename of, the parameter file
}
  \item{start.time}{
numeric; CPU start time of the calling procedure, in seconds. Used for output of elapsed times; if not supplied, current CPU time is used.
}
  \item{quiet}{
logical; if TRUE, all output is suppressed.
}
  \item{env}{
environment; the environment to which all parameters will be assigned. if NULL, a new environment is generated.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
NULL; all output is to the parameter environment `env'
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
\link{measure.fluxes}, \link{create.par.file}
}
\examples{

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Parameters }
