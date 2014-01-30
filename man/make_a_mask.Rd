\name{make_a_mask}
\alias{make_a_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a Full Mask
}
\description{
Create a full image mask of all smaller masks/apertures which are input. The full mask is of dimension fullmaskdim.
Pixel locations for the placement of smaller masks on large mask are determined from the
input variables x, y.
}
\usage{
make_a_mask(env=NULL, masks, fullmaskdim, outenv=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
environment; the environment containing the parameters used in 'MeasureFluxes'
}
  \item{masks}{
numeric length(m) list of (n,n) arrays; length (m) list of (n,n) stamps, each containing a mask.
}
  \item{fullmaskdim}{
numeric (2) vector; fullmaskdim[1] is the x-dimension of the full mask. fullmaskdim[2] is the y-dimension. fullmaskdim values must be > submask width.
}
  \item{outenv}{
environment; the environment where the returned values are to be stored. If not provided, outenv defaults to the input env.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
array; full mask, with dimesions fullmaskdim, of smaller masks placed in their appropriate locations.
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
