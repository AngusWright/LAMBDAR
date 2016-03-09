\name{lsos}
\alias{lsos}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
R Object List
}
\description{
A better object list
}
\usage{
lsos(pattern="", envir=NULL, pos = 1, order.by="Size",
decreasing=TRUE, head=TRUE, n=10,noShow="",exact=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{pattern}{
character; the pattern to search for in the object list
}
  \item{envir}{
environment; the environment to list
}
  \item{pos}{
integer; environment position
}
  \item{order.by}{
character; field name to order object list by. Options are "Type", "Size", "PrettySize", "Rows", "Columns", "Value".
}
  \item{decreasing}{
logical; if TRUE, objects are listed in decreasing order
}
  \item{head}{
logical; if TRUE, only the first `n' objects are listed
}
  \item{n}{
integer; number of items to list when `head' is specified
}
  \item{noShow}{
character vector; fields to not show. Can be any combination of "type","size","pretty","dim","vals".
}
  \item{exact}{
logical; only show items that match 'pattern' exactly.
}
}
\details{
A simple procedure to list objects matching `pattern' in enviroment at `pos', including
details on objects type, size, and dimensions.
}
\value{
data frame;
\item{Type}{ Object class/type}
\item{Size}{ Size of the object (bits)}
\item{PrettySize}{ Size of the object in nearest major unit }
\item{Rows}{ Number of rows in the object}
\item{Columns}{ Number of columns in the object}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Procedure was authored by Dirk Eddelbuettel and Tony Breyal in the StackOverflow Forum:
http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
It has been modified by A.H.Wright.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Add some things to the an environment
param.env<-new.env(parent=environment())
assign("x",runif(1),envir=param.env)
assign("y",runif(1E5),envir=param.env)
assign("z",runif(1E2),envir=param.env)
assign("string","Methinks it is like a weasel",envir=param.env)
assign("im",matrix(runif(1E6),1000,1000),envir=param.env)
assign("bigmat",matrix(runif(1E8),10000,10000),envir=param.env)

#Show the parameter space variables
lsos(envir=param.env)

#Only show the type and pretty size
lsos(envir=param.env, noShow=c("size","dim","vals"))

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Parameters }
