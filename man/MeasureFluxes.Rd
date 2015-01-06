\name{MeasureFluxes}
\alias{MeasureFluxes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure Fluxes
}
\description{
Measure Photometric fluxes from an arbitrary FITS image, when provided an appropriate input parameter file.
Detailed descriptions of the many functions called by this procedure can be found in their relevant help files.
The procedure has many high-level options available, all of which are designated in the input parameter file.
A default parameter file can be created by running MeasureFluxes("--makepar"). For more in depth examples, see below.
}
\usage{
MeasureFluxes(parfile, quiet = FALSE, MPIBackend=FALSE, doReturn=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parfile}{
Path to, and name of, the Input Parameter file (ie. '/path/to/parfile.par').
}
  \item{quiet}{
logical; if TRUE all output to screen is supressed
}
  \item{...}{
Any arguments that are intended to be passed to the 'createparfile' routine; see examples below.
}
}
\details{
MeasureFluxes is the priciple routine in the LAMBDAR package. The LAMBDAR package is intended for use
in astronomy, and other related fields, to determine photometry from input FITS photometric images of
arbitrary resolution.

The program has many optional functions, and can be run in a standard threaded mode or in a Cluster-enabled
MPI mode. Examples of how to run the program utilising the MPI backend are supplied in the examples below.
The many options available to the program are outined in the CreateParFile function's manual page, as well
as being explained in greater detail in related manual pages.
}
\value{
MeasureFluxes will return either NULL (if doReturn=FALSE) or a list contining a subset the photometric results
(if doReturn=TRUE). By Default, MeasureFluxes will return NULL and all outputs are to file. If requested
the returned list contains the following:
  \item{SFAflux}{
  The flux of each aperture, without consideration of deblending, in image units.
}
  \item{SFAerror}{
  The error on the un-deblended flux, in image units.
}
  \item{DFAflux}{
  The final deblended flux of each aperture, in image units.
}
  \item{DFAerror}{
  The error associated with the final deblended flux for each aperture, in image units.
}
Additionally, if the Diagnostic flag is 1 in the parameter file, doReturn=TRUE will also return the entire aperture stamps lists:

  \item{SA_Stamps}{
  The list of aperture stamps genreated by the make_sa_mask function
}
  \item{SFA_Stamps}{
  The list of aperture stamps genreated by the make_sfa_mask function without fluxweighting
}
  \item{WSFA_Stamps}{
  The list of aperture stamps genreated by the make_sfa_mask function with fluxweighting active
}
  \item{DFA_Stamps}{
  The list of aperture stamps genreated by the make_sfa_mask function, and after deblending
}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Angus H Wright{ ICRAR}
 {\email{angus.wright@icrar.org}}
}
%\note{
%}
\note{
All parameters and data used by MeasureFluxes are stored in environments created within the function.
At completion of the routine, these environments are all discarded, meaning that there is no modification
to the runtime environment (as per CRAN requirements).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{createparfile}
}
\examples{
# To generate a default parameter file for use:
MeasureFluxes("--makepar")

# When generating the parameter file, all parameters can be made custom
# by including the parameter name and its desired value in the makepar call.
# For Example, to generate a parameter file for use, with customised Catalogue
# name value:
MeasureFluxes("--makepar", Catalogue="MyCat.csv")

# Or changing the input Catalogue name, and changing number of foreach/doParallel processors
MeasureFluxes("--makepar", Catalogue="MyCat.csv", nProcessors=10)

# To execute the program using the parameter file "Lambdar.par":
MeasureFluxes("Lambdar.par")

# To execute the program using the MPI backend (provided that you have a functional
# System installation of MPI, as well as the relevant R packages Rmpi, and doMPI)
# you would execute:

mpi -np <NUM_PROCESSORS> R --slave -e "MeasureFluxes('/path/to/parfile.par',MPIBackend=TRUE)"

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
