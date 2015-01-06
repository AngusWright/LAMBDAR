\name{CreateSim}
\alias{CreateSim}
\alias{create.sim.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Simulated Image
}
\description{
Create Simulate Photometric Image from an arbitrary FITS image. This function takes the input
image and catalogue from the MeasureFluxes Parameter File and uses them to generate a simulated
image. Objects that are point sources are modelled with a perfect analytical gaussian, with FWHM
as specified by the user or as determined empirically from the input PSF.
create.sim.image is a subroutine within CreateSim and is not designed to be called directly.
}
\usage{
CreateSim(parfile, ObsParm=NA, noNoise=FALSE, convolveNoise=TRUE, padGalaxies=TRUE, colourCorr=0.0, quiet=FALSE)
create.sim.image(ObsParm, noNoise=FALSE, convolveNoise=TRUE, padGals=TRUE, outenv=parent.env(environment()), env=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parfile}{
Path to, and name of, the Input Parameter file (ie. '/path/to/parfile.par'). This is of the same
as the MeasureFluxes Parameter File, however only some parameters are relevant.
}
  \item{ObsParm}{
List containing Observation Parameters needed in creating simulated galaxies; Exposure Time (sec; exp=),
Collecting Area (m^2; area=), Effective wavelength of filter (Å; lamEff=), and Effective width of filter (Å; WEff=).
For example, SDSS r-band has parameters ObsParm=list(exp=53.9,area=4.91,lamEff=1111.2,Weff=6165).
}
  \item{noNoise}{
logical; if TRUE then noise is not added to the simulated image
}
  \item{convolveNoise}{
logical; if TRUE then noise is convolved with the PSF prior to addition
}
  \item{padGalaxies}{
logical; if TRUE then the simulated image will have galaxies added beyond those in the input catalogue, such that
the galactic number counts follow the power law log(N)=0.38*m_r-4.78.
}
  \item{colourCorr}{
numeric; The colour correction need to convert SDSS r-band magnitude to magnitude in the band being analysed.
This is needed when padding galaxies because our power law is defined in the r-band. Naturally, this is a 0th
order correction, and should be noted when determining the accuracy of simulations that are well beyond the
SDSS r-band.
}
  \item{quiet}{
logical; if TRUE all output to screen is supressed
}
}
\details{
CreateSim is an additional routine included in the LAMBDAR package for the purposes of testing.
The LAMBDAR package is intended for use in astronomy, and other related fields, to determine photometry from
input FITS photometric images of arbitrary resolution. As such, testing of the program is paramount. This
function allows the user to generate simulated images from input galactic apertures and fluxes. Using the
known input, the user may then test the MeasureFluxes routine with all its various options permutted, thus
testing the programs accuracy. Results of this type of rigorous testing are demonstrated in Wright et al. 2015 (in prep).
}
\value{
NULL; CreateSim will output a simulated image to "sim_image.fits", and a catalogue of the created objects
called "SimFlux.csv".
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{createparfile}
}
\examples{

#Create simulated image using parameter file '/path/to/parfile.par' and SDSS r-band observation parameters:
CreateSim('/path/to/parfile.par',ObsParm=list(exp=53.9,area=4.91,lamEff=1111.2,Weff=6165),convolveNoise=T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
