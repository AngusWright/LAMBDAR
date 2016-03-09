\name{lastnchar}
\alias{lastnchar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
String Manipulations
}
\description{
String functions; lastnchar returns the last `n' characters of the input string
}
\usage{
lastnchar(x, n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{
string; length > n
}
  \item{n}{
integer; number of characters to be returned
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
string; length (n)
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
\link{MeasureFluxes}, \link{readparfile}
}
\examples{

#Show the last 6 characters of a random string
string<-"Methinks it is like a weasel"
print(lastnchar(string,6))

#Is the last character of the filename a '/'?
path<-getwd()
#Using getwd alone: No
lastnchar(path,1)=="/"
#Using getwd and file.path: Yes
lastnchar(file.path(path,""),1)=="/"

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Parameters }
