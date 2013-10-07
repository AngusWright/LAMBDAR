MeasureFluxes <-
function(parfile=NA, quiet=FALSE, ...){ #Begin
  #
  # Proceedure measures GAMA object fluxes from an arbitrary fits image
  #
  
  #For Setup, warnings are handled internally - print nothing 
  options(warn=-1)
  
# Create the LAMBDAR_WorkEnv environment, which will contain all 
  # Data for the LAMBDAR Proceedure
  #assign(".w.env",environment())  
#
  ##Set the LAMBDAR_Work Environment as default
  #environment(readparfile)<-environment()
  #environment(crop_im)<-.w.env
  #environment(readimage)<-.w.env
  #environment(opencatalogue)<-.w.env
  #environment(fluxmeasurements)<-.w.env
  #environment(make_sa_mask)<-.w.env
  #environment(make_a_mask)<-.w.env
  #environment(make_sfa_mask)<-.w.env
  #environment(make_deblended_weightmap)<-.w.env
  #environment(sourcesubtraction)<-.w.env

  #Set quiet variable
  assign("quiet", quiet)

  #Check for appropriate calling syntax
  if (is.na(parfile)) {
  stop(paste("Parameter file not supplied.\n",
             "Calling Syntax:\n",
             "       MeasureFluxes(<ParameterFile>,<QuietFlag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run MeasureFluxes('--makepar').", sep=""))
  }

  #If requested, produce the default .par file and end
  if (parfile == "--makepar") {
    if (!quiet) { cat("Outputting Default Parameter file to './Lambdar_default.par'\n") }
    createparfile(...)
    if (!quiet) { cat("Program Complete\n") }
    return()
  }

  #Check that required packages are installed
  if (require(foreach, quietly=TRUE)==FALSE) { stop("Required Package 'foreach' not installed") }
  if (require(doParallel, quietly=TRUE)==FALSE) { stop("Required Package 'doParallel' not installed") }
  if (require(astro, quietly=TRUE)==FALSE) { stop("Required Package 'astro' not installed") }

  #Set start timer & print opening
  starttime<-proc.time()[3]
  if (!quiet) { cat(paste('START: ',starttime,'==================================\n')) }

  #Setup Parameter Space (read .par file)
  readparfile(parfile,starttime,quiet,env=environment())

  #print(ls())
  #print("Global:")
  #print(ls(.GlobalEnv))

  #From here on, produce warnings as they occur
  options(warn=1)

  #If wanted, crop image prior to read
  if (cropimage) {
    if (verbose) { message(paste("Cropping Input Image: Outputting to", imfitsoutname)) }
    crop_im(ra0=ra0, dec0=dec0, pathroot=pathroot, inpim=datamap, cutrad=cutrad, fitsoutname=imfitsoutname)
    if (verbose) { message(paste("Using", imfitsoutname, "as data image")) }
    datamap<-imfitsoutname
    if (maskmap != "NONE") {
      if (verbose) { message(paste("Cropping Input Mask Map: Outputting to", immfitsoutname)) }
      crop_im(ra0=ra0, dec0=dec0, pathroot=pathroot, inpim=maskmap, cutrad=cutrad, fitsoutname=immfitsoutname)
      if (verbose) { message(paste("Using", immfitsoutname, "as mask map")) }
      maskmap<-immfitsoutname
    }
    if (errormap != "NONE") {
      if (verbose) { message(paste("Cropping Input Error Map: Outputting to", imefitsoutname)) }
      crop_im(ra0=ra0, dec0=dec0, pathroot=pathroot, inpim=errormap, cutrad=cutrad, fitsoutname=imefitsoutname)
      if (verbose) { message(paste("Using", imefitsoutname, "as error map")) }
      errormap<-imefitsoutname
    }
  }

  # Create the image environment, which will contain all 
  # Image Data for the LAMBDAR Proceedure - This segregates the 
  # large image arrays from the regular parameter space, and 
  # stops unnecessary memory usage in foreach loops
  image.env<-new.env(parent=environment())

  #Read in Data, Mask map, & Error map
  readimage(environment(),quiet,showtime,image.env)
  astr_struc<-image.env$astr_struc
  rm(astr_struc, envir=image.env)

  #If needed, read ZP Magnitude from Image header
  if ((makeMagnitudes) & (magZP==-999)){
    magZP<-read.fitskey(magZPlabel,paste(pathroot,datamap,sep="")) 
    if (!is.finite(magZP)) { sink(type='message'); stop("Zero Point Magnitude determination failed") }
  }
  
  #Read source catalogue
  opencatalogue(environment())

  #Determine which galaxies are inside the input image
  #and have physical aperture parameters
  if (!quiet) { cat("   Determining Correct Galaxy Sample ") }

  gamapos<-ad2xy(ra_g,dec_g,astr_struc)
  x_g<-gamapos[,1]
  y_g<-gamapos[,2]

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("X_G:", length(x_g),"\nY_G:",length(y_g)))
    message(paste("min/max X_G:", min(x_g), max(x_g),"\nmin/max Y_G:",min(x_g), max(y_g)))
  }

  #Discard any apertures that lie completely outside of the image ± 1 pixel 
  catlen<-length(x_g)
  insidemask<-!((x_g <= 0) | (x_g >= length(image.env$im[,1])+1) | (y_g <= 0) | (y_g >= length(image.env$im[1,])+1)) 
  if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Single Apertures are inside the image.") }  # Nothing inside the image
  x_g<-x_g[which(insidemask)]
  y_g<-y_g[which(insidemask)]
  id_g<-id_g[which(insidemask)]
  ra_g<-ra_g[which(insidemask)]
  dec_g<-dec_g[which(insidemask)]
  theta_g<-theta_g[which(insidemask)]
  a_g<-a_g[which(insidemask)]
  b_g<-b_g[which(insidemask)]
  if (filtcontam) { contams<-contams[which(insidemask)] }
  if (verbose) { message(paste("There are",length(x_g),"supplied objects inside the image (",
                                round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied were outside the image )")) }
  #Discard any apertures that have nonphysical aperture axis values
  catlen<-length(x_g)
  insidemask<-!((a_g < 0)|(b_g < 0))
  if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Apertures remaining have physical axis values.") } 
  x_g<-x_g[which(insidemask)]
  y_g<-y_g[which(insidemask)]
  id_g<-id_g[which(insidemask)]
  ra_g<-ra_g[which(insidemask)]
  dec_g<-dec_g[which(insidemask)]
  theta_g<-theta_g[which(insidemask)]
  a_g<-a_g[which(insidemask)]
  b_g<-b_g[which(insidemask)]
  if (filtcontam) { contams<-contams[which(insidemask)] }
  if (verbose) { message(paste("There are",length(x_g),"supplied objects with physical aperture values (",
                                round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied had unphysical values )")) }

  if (!quiet) { cat(" - Done\n") }
  
  #Set object Astrometry
  if (!quiet) { cat("   Setting Astrometry   ") }

  #Get pixel resolution from astrometry
  if (all(is.finite(astr_struc$CDELT))) {
    #Using CDELT keywords
    asperpix<-max(abs(astr_struc$CDELT))*3600. 
  } else if (all(is.finite(astr_struc$CD))) { 
    #Using CD matrix keywords
    asperpix<-max(astr_struc$CD[1,1],astr_struc$CD[2,2])*3600. 
  } else {
    #ERROR: no pixel width keywords
    sink(type=c("output","message")) ; stop("Data image header does not contain CD or CDELT keywords")
  }
  
  #Set apertures with NA/NULL aperture axis or minoraxis<apertruediag to point-sources
  diag_arcsec<-abs(asperpix)*sqrt(2)
  message("Forcing",length(which((b_g<diag_arcsec)|!is.finite(a_g))),"apertures to be point sources")
  a_g[which((b_g<diag_arcsec)|!is.finite(a_g))]<-0
  theta_g[which((b_g<diag_arcsec)|!is.finite(theta_g))]<-0
  b_g[which((b_g<diag_arcsec)|!is.finite(b_g))]<-0

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("X_G:", length(x_g),"\nY_G:",length(y_g)))
    message(paste("min/max X_G:", min(x_g), max(x_g),"\nmin/max Y_G:",min(x_g), max(y_g)))
  }

  #Fix beam area in pixels
  beamarea_pix<-beamarea_SOM_as/(asperpix)^2.

  #Finished setting astrometry
  if (!quiet) { 
    cat(" - Done\n") 
    cat('} Initialisation Complete ') 
  }
  if (showtime) { cat(paste(' (  Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'  )\n')) }

  #-----Diagnostic-----#
  if (diagnostic) {
  message(paste('Arcsec per pixel in map: ', asperpix))
  message(paste('Beam area (from observers manual) converted into pixel units: ', beamarea_pix))
  } 
  #Enter flux measurements
  results<-fluxmeasurements(environment())
  #fluxmeasurements(id_g,ra_g,dec_g,x_g,y_g,theta_g,a_g,b_g,im,hdr_str,ime,imm,pathroot,pathout,
  #psfmap,lambda,field,conf,fluxcorr,fluxweight,beamarea_pix,asperpix,Jybm,sourcemask,nopsf,
  #filtcontam,contams)

   #Program Completed - Print closing and remove sink
   if (!quiet) { cat(paste('-----------------------------------------------------\nProgram Complete\n'))
                 cat(paste('Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n')) 
                 message(paste('-----------------------------------------------------\nProgram Complete\n'))
                 message(paste('Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n')) }
   sink(type="message")
   return=results

}
