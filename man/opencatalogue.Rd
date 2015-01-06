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
opencatalogue(env=NULL,outenv=NULL)
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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
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
