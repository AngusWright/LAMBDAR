\name{make_a_mask}
\alias{make_a_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Full Mask
}
\description{
Create a full image mask of all smaller masks which are input. The full mask is of dimension fullmaskdim.
Pixel locations for the placement of smaller masks on large mask are determined from the
input variables x, y.
}
\usage{
make_a_mask(masks, fullmaskdim, x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{masks}{
numeric length(m) list of (n,n) arrays; length (m) list of (n,n) stamps, each containing a mask.
}
  \item{fullmaskdim}{
numeric (2) vector; fullmaskdim[1] is the x-dimension of the full mask. fullmaskdim[2] is the y-dimension. fullmaskdim values must be > submask width.
}
  \item{xcen}{
numeric (m) vector; central x positions of masks to be placed onto full mask (pixel number)
}
  \item{ycen}{
numeric (m) vector; central y positions of masks to be placed onto full mask (pixel number)
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
MeasureFluxes, fluxmeasurements
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
