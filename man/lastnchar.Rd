\name{last.n.char}
\alias{last.n.char}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
String Manipulations
}
\description{
String functions; last.n.char returns the last `n' characters of the input string
}
\usage{
last.n.char(x, n)
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
\link{measure.fluxes}, \link{read.par.file}
}
\examples{

#Show the last 6 characters of a random string
string<-"Methinks it is like a weasel"
print(last.n.char(string,6))

#Is the last character of the filename a '/'?
path<-getwd()
#Using getwd alone: No
last.n.char(path,1)=="/"
#Using getwd and file.path: Yes
last.n.char(file.path(path,""),1)=="/"

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Parameters }
