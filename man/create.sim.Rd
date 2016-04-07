\name{create.sim}
\alias{create.sim}
\alias{create.sim.image}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{
Create Simulated Image
}
\description{
Create Simulate Photometric Image from an arbitrary FITS image. This function takes the input
image and catalogue from the measure.fluxes Parameter File and uses them to generate a simulated
image. Objects that are point sources are modelled with a perfect analytical gaussian, with FWHM
as specified by the user or as determined empirically from the input PSF.
create.sim.image is a subroutine within create.sim and is not designed to be called directly.
}
\usage{
create.sim(par.file=NA, ObsParm=NA, noNoise=FALSE,convolveNoise=TRUE,
padGalaxies=TRUE, colourCorr=0, quiet=FALSE, confuse=FALSE)
create.sim.image(ObsParm, noNoise=FALSE, convolveNoise=TRUE,
padGals=TRUE, col.corr=0,outenv=parent.env(environment()),
confuse=FALSE, env=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{par.file}{
Path to, and name of, the Input Parameter file (ie. `/path/to/par.file.par'). This is of the same
as the measure.fluxes Parameter File, however only some parameters are relevant.
}
  \item{ObsParm}{
Data frame containing Observation Parameters needed in creating simulated galaxies; Exposure Time (sec; exp=),
Collecting Area (meters squared; area=), Effective wavelength of filter (Angstrom; lamEff=), and Effective width of filter (Angstrom; WEff=).
For example, SDSS r-band has parameters
ObsParm=data.frame(exp=53.9, area=4.91, lamEff=1111.2, Weff=6165).
}
  \item{noNoise}{
logical; if TRUE then noise is not added to the simulated image
}
  \item{convolveNoise}{
logical; if TRUE then noise is convolved with the PSF prior to addition
}
  \item{padGalaxies}{
logical; if TRUE then the simulated image will have galaxies added beyond those in the input catalogue, such that
the galactic number counts follow the polynomial log(N)=-0.01476*m_r^2+1.025*m_r-11.76.
}
  \item{colourCorr}{
numeric; The colour correction need to convert SDSS r-band magnitude to magnitude in the band being analysed.
This is needed when padding galaxies because our power law is defined in the r-band. Naturally, this is a 0th
order correction, and should be noted when determining the accuracy of simulations that are well beyond the
SDSS r-band.
}
  \item{confuse}{
logical; do you want to confuse the simulation such that the source density is greater than 1 object per FWHM^2
}
  \item{col.corr}{
same as the colourCorr argument above
}
  \item{padGals}{
same as the padGalaxies argument above
}
  \item{outenv}{
environment; parameter space to output variables to
}
  \item{quiet}{
logical; if TRUE all output to screen is suppressed
}
  \item{env}{
environment;; parameter space containing the input measure.fluxes/LAMBDAR core variables
}
}
\details{
create.sim is an additional routine included in the LAMBDAR package for the purposes of testing.
The LAMBDAR package is intended for use in astronomy, and other related fields, to determine photometry from
input FITS photometric images of arbitrary resolution. As such, testing of the program is paramount. This
function allows the user to generate simulated images from input galactic apertures and fluxes. Using the
known input, the user may then test the measure.fluxes routine with all its various options permuted, thus
testing the programs accuracy. Results of this type of rigorous testing are demonstrated in Wright et al. 2015 (in prep).
}
\value{
NULL; create.sim will output a simulated image to "sim_image.fits", and a catalogue of the created objects
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
\link{create.par.file}\link{lsos}\link{write.fits}\link{measure.fluxes}
}
\examples{

#Create simulated image using an SDSS r-band image and observation parameters:
#Load a Sample Image and Catalogue
data("SDSS.sample",envir=environment())
data("ApCat.sample",envir=environment())
lsos(envir=environment())
#Write Sample Image to File
write.fits(file="SampleImage.fits",SDSS.sample)
#Write Sample Catalogue to File
write.csv(file="SampleCat.csv",ApCat.sample,row.names=FALSE,quote=FALSE)
#Generate Parameter File
create.par.file(file="Lambdar_Sample.par",DataMap="SampleImage.fits",
Catalogue="SampleCat.csv",IterateFluxes=1,nIterations=2,Gauss_FWHM_AS=2,
PSFConvolve=1)

#Measure Photometry
measure.fluxes("Lambdar_Sample.par")

#Generate Simulation Parameter File
create.par.file(file="Lambdar_Simulation.par",DataMap="SampleImage.fits",
Catalogue="LAMBDAR_Results.csv",Gauss_FWHM_AS=2,SimGauss_AS=1.1,
FluxWgtColumnLabel="DFAFlux_Jy",OutputDirectory="Simulation/",PSFConvolve=1)

#Create Simulation from Catalogued objects only
create.sim('Lambdar_Simulation.par',ObsParm=data.frame(exp=53.9,area=4.91,
lamEff=1111.2,Weff=6165),convolveNoise=TRUE,padGalaxies=TRUE)

#Read the Input Image
inp.im<-read.fits("SampleImage.fits")
#Read the Simulated Image
sim.im<-read.fits("Simulation/sim_image.fits")

#Plot the Input and Simulation Images
layout(matrix(1:4,2,2,byrow=TRUE))
#Get Object Positions
xy<-ad.to.xy(ApCat.sample$RAdeg,ApCat.sample$DECdeg,
read.astrometry("SampleImage.fits"))
aspp<-abs(read.astrometry("SampleImage.fits")$CD[1,1])*3600
#Plot Input & Aperture Ellipses
image(x=1:length(inp.im$dat[[1]][,1]),y=1:length(inp.im$dat[[1]][1,]),
inp.im$dat[[1]],col=grey.colors(1E3),zlim=quantile(inp.im$dat[[1]],
c(0,0.999)),asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)")
invisible(Vectorize(function(x,y,a,b,pa) {
   lines(ellipse(x,y,a,b,pa=pa),col='red',lty=2)
   points(x,y,pch=3,cex=0.5,col='red')
})(xy[,1],xy[,2],ApCat.sample$radmajasec/aspp,ApCat.sample$radminasec/aspp,
ApCat.sample$rotW2N))
#Plot Simulation & Aperture Ellipses
image(x=1:length(sim.im$dat[[1]][,1]),y=1:length(sim.im$dat[[1]][1,]),
sim.im$dat[[1]],col=grey.colors(1E3),zlim=quantile(inp.im$dat[[1]],
c(0,0.999)),asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)")
invisible(Vectorize(function(x,y,a,b,pa) {
   lines(ellipse(x,y,a,b,pa=pa),col='red',lty=2)
   points(x,y,pch=3,cex=0.5,col='red')
})(xy[,1],xy[,2],ApCat.sample$radmajasec/aspp,ApCat.sample$radminasec/aspp,
ApCat.sample$rotW2N))
#Plot Zoomed Sections of the Image:
#-> Input & Aperture Ellipses
image(x=1:length(inp.im$dat[[1]][,1]),y=1:length(inp.im$dat[[1]][1,]),
inp.im$dat[[1]],col=grey.colors(1E3),zlim=quantile(inp.im$dat[[1]],
c(0,0.999)),asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)",
xlim=c(800,1200),ylim=c(800,1200))
points(xy[,1],xy[,2],pch=3,cex=0.5,col='red')
#-> Simulation & Aperture Ellipses
image(x=1:length(sim.im$dat[[1]][,1]),y=1:length(sim.im$dat[[1]][1,]),
sim.im$dat[[1]],col=grey.colors(1E3),zlim=quantile(inp.im$dat[[1]],
c(0,0.999)),asp=1,useRaster=TRUE,xlab='X (pix)',ylab="Y (pix)",
xlim=c(800,1200),ylim=c(800,1200))
points(xy,pch=3,cex=0.5,col='red')

#Measure Simulation Fluxes
create.par.file(file="Lambdar_SimAnalysis.par",
DataMap="Simulation/sim_image.fits",Catalogue="Simulation/SimFlux.csv",
Gauss_FWHM_AS=2,OutputDirectory="SimResults/",PSFConvolve=1,BlankCor=1)
measure.fluxes('Lambdar_SimAnalysis.par')

#Read simulation results
sim.out.cat<-read.csv("SimResults/LAMBDAR_Results.csv",stringsAsFactors=FALSE)
sim.in.cat.full<-read.csv("Simulation/SimFlux.csv",stringsAsFactors=FALSE)
sim.in.cat<-sim.in.cat.full[which(sim.in.cat.full$CATAID \%in\% sim.out.cat$CATAID),]

#Plot Flux Comparisons
layout(c(1,2))
magplot(sim.in.cat$Sim.Flux.Jy, sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy,
pch=20, log='xy', ylim=c(0.5,2),majorn=10,minorn=10,ylab='Measured / Input Flux',
xlab='Input Flux (Jy)')
#Add Uncertainties
suppressWarnings(magerr(x=sim.in.cat$Sim.Flux.Jy,
y=sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy,
ylo=abs(sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy)*
sim.out.cat$DFAErr_Jy/sim.out.cat$DFAFlux_Jy))
abline(h=1,col='grey')
abline(h=0.96,col='blue',lty=2)
legend('topright',legend=c("Expected Fraction\n(2.5Kron w/ exponential profile)"),
col='blue',lty=2)

#Improve Flux Measurement by Iterating
create.par.file(file="Lambdar_SimAnalysis.par",
DataMap="Simulation/sim_image.fits",Catalogue="Simulation/SimFlux.csv",
Gauss_FWHM_AS=2,OutputDirectory="SimResults/",IterateFluxes=1,
PSFConvolve=1,nIterations=5,BlankCor=1)
measure.fluxes('Lambdar_SimAnalysis.par')

#Read simulation results
sim.out.cat<-read.csv("SimResults/LAMBDAR_Results.csv",stringsAsFactors=FALSE)
sim.in.cat<-sim.in.cat.full[which(sim.in.cat.full$CATAID \%in\% sim.out.cat$CATAID),]

#Plot Flux Comparisons
magplot(sim.in.cat$Sim.Flux.Jy, sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy,
pch=20, log='xy', ylim=c(0.5,2),majorn=10,minorn=10,
ylab='Measured / Input Flux',xlab='Input Flux (Jy)')
#Add Uncertainties
suppressWarnings(magerr(x=sim.in.cat$Sim.Flux.Jy,
y=sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy,
ylo=abs(sim.out.cat$DFAFlux_Jy / sim.in.cat$Sim.Flux.Jy)*
sim.out.cat$DFAErr_Jy/sim.out.cat$DFAFlux_Jy))
abline(h=1,col='grey')
abline(h=0.96,col='blue',lty=2)
legend('topright',legend=c("Expected Fraction\n(2.5Kron w/ exponential profile)"),
col='blue',lty=2)

}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{ Simulation }
