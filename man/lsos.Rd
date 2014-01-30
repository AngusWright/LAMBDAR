\name{lsos}
\alias{lsos}
\alias{ls.objects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
R Object List
}
\description{
A better object list
}
\usage{
lsos(..., n = 10)
ls.objects(pattern="", pos = 1, order.by, decreasing=FALSE, head=FALSE, n=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
arguments to be sent to 'ls.objects'
}
  \item{n}{
integer; number of items to list
}
}
\details{
A simple procedure to list objects matching 'pattern' in enviroment at 'pos', including
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
Prodedure was authored by Dirk Eddelbuettel and Tony Breyal in the StackOverflow Forum:
http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., n = 10)
{
    .ls.objects(..., order.by = "Size", decreasing = TRUE, head = TRUE,
        n = n)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
