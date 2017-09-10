\name{get.nth.var}
\alias{get.nth.var}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Parameter Assignment to Environment
}
\description{
Assign the Nth variable from a parameter list to a designated environment
}
\usage{
 get.nth.var(varList,n,inenv,outenv,lim)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{varList}{
character vector; list of parameters to assign to output environment.
}
  \item{n}{
integer; which element to assign from each parameter vector
}
  \item{inenv}{
environment; the parameter space containing the full parameter vectors
}
  \item{outenv}{
environment; the parameter space to which we want to output the variables. _Cannot_ be the same as inenv
}
  \item{lim}{
numeric; the maximum number of loops to be performed
}
}
\details{
The function takes a parameter space, inenv, and for each parameter designated in varList assigns its
'n'th component to the outenv parameter space. When assigning parameters, the function behaves as follows:
If a parameter is of length 1, that value is returned, regardless of `n'. If a parameter of length `lim', then
the `n'th value is returned. If the parameter has length such that length(par)\%\%'lim'==0, then the `n'\%\%'lim'th value is
returned.
}
\value{
NULL; all parameters are assigned to the outenv environment
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Load LAMBDAR
library(LAMBDAR)
#Create the 'parameters' environment
parameter.env<-new.env(parent=environment())
#Assign a bunch of values there, of different lengths
assign("par1",1:12,env=parameter.env)
assign("par2",c("a","b","c"),env=parameter.env)
assign("par3",11:22,env=parameter.env)
assign("par4","xyz",env=parameter.env)
assign("par5",function(x) print(x),env=parameter.env)

#Assign the 2nd value of parameters 1, 2, and 5 to this environment
get.nth.var(c("par1","par2","par4"),2,inenv=parameter.env,
outenv=environment(),lim=12)
#Show Parameters
lsos(envir=environment())

#Assign the 5th value of all parameters to this environment
get.nth.var(ls(envir=parameter.env),2,inenv=parameter.env,
outenv=environment(),lim=12)
#Show Parameters
lsos(envir=environment(),noShow=c('size','pretty'))


}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Parameters }
