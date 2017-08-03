\name{as.time}
\alias{as.time}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Convert time (sec) into human readable string (day hr min sec) 
}
\description{
Converts an elapsed time (in seconds) into a human readable string. String 
is output in units decreasing from the largest integer unit of day/hour/min/sec. 
If the time is greater than 1 day, seconds are not printed. 
}
\usage{
as.time(sec,digits)
}
%- maybe also `usage' for other objects documented here.
\arguments{
\item{sec}{Elapsed time in seconds }
\item{digits}{Number of decimal places to use when reporting smalled time unit }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
string: time elapsed in Day, Hr, Min, Sec. 
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

elapsed.time <- Sys.time()
for (i in 1:6) { 
  print(as.time(10^i,digits=3))
}


}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Plot }% __ONLY ONE__ keyword per line
