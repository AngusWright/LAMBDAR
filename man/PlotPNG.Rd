\name{PlotPNG}
\alias{PlotPNG}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Open PNG plot device
}
\description{
Opens a Cairo PNG plot device if available, otherwise opens a standard R png plot device. 
}
\usage{
PlotPNG(...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{...}{Arguments to pass to either CairoPNG() or png() }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
device identifier 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ AIfA}
 {\email{awright@astro.uni-bonn.de}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

#Do not Run {{{
# PlotPNG()
# plot(rnorm(1E4),rnorm(1E4))
# dev.off()
#}}}

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Plot }% __ONLY ONE__ keyword per line
