CreateSim <-
function(parfile=NA, ObsParm=NA, noNoise=FALSE, convolveNoise=TRUE, padGalaxies=TRUE, colourCorr=0.0, quiet=FALSE, confuse=FALSE){
#Proceedure measures GAMA object fluxes from an arbitrary fits image

  #For Setup, warnings are handled internally - print nothing {{{
  options(warn=-1)
  #}}}

  #Set function Environments {{{
  environment(opencatalogue)<-environment()
  environment(readimage)<-environment()
  environment(readparfile)<-environment()
  environment(create.sim.image)<-environment()
  #}}}

  #Check for appropriate calling syntax {{{
  if (is.na(parfile)) {
  stop(paste("Parameter file not supplied.\n",
             "Calling Syntax:\n",
             "       MeasureFluxes(<ParameterFile>,<Observation Parameters>,<Quiet Flag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<Observation Parameters> = Data frame containing the details of the observation of the input image: \n",
             "            $exp = Effective Exposure Time (s)\n",
             "            $area = Telescope Collecting Area (m^2)\n",
             "            $lamEff = Filter Effective Central Wavelength (Ang)\n",
             "            $Weff = Filter Effective Width (Ang)\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run MeasureFluxes('--makepar').", sep=""))
  }
  if (is.na(ObsParm)|!is.data.frame(ObsParm)) {
  warning(paste("Observation Parameters not supplied, or supplied incorrectly.\n",
             "Calling Syntax:\n",
             "       MeasureFluxes(<ParameterFile>,<Observation Parameters>,<Quiet Flag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<Observation Parameters> = Data frame containing the details of the observation of the input image: \n",
             "            $exp = Effective Exposure Time (s)\n",
             "            $area = Telescope Collecting Area (m^2)\n",
             "            $lamEff = Filter Effective Central Wavelength (Ang)\n",
             "            $Weff = Filter Effective Width (Ang)\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run MeasureFluxes('--makepar').", sep=""))
  warning("Using SDSSIII r-band default parameters for observation paramters\n")
  ObsParm=data.frame(exp=53.9,area=pi*(2.5/2)^2,lamEff=1111.2,Weff=6165)
  }#}}}

  #If requested, resume LAMBDAR run from last loop state {{{
  if (parfile == "--resume") {
    if (file.exists(".LambdarParameters.Rdata")&file.exists(".LambdarActiveLoop.txt")) {
      if (!quiet) { cat("Resuming from Previous Loop State\n") }
      load(".LambdarParameters.Rdata")
      lstart<-as.numeric(read.csv(".LambdarActiveLoop.txt"))[1]
      resume<-TRUE
    } else {
      stop("Resume Requested, but one/both of the required resume files .LambdarParameters.Rdata & .LambdarLoopState.txt are missing")
    }
  } else {
    lstart<-1
    resume<-FALSE
  }#}}}

  #Set start timer & print opening {{{
  starttime<-proc.time()[3]
  #}}}

  memLim<-Inf

  #Setup Parameter Space (read .par file) {{{
  if (!resume) {
    param.env<-new.env(parent=environment())
    readparfile(parfile,starttime,quiet,env=param.env)
  }
  parameterList<-ls(envir=param.env)
  #}}}

  #From here on, produce warnings as they occur {{{
  options(warn=1)
  #}}}

  #If needed, register the parallel backend {{{
  registerDoParallel(cores=param.env$ncores)
  if (!quiet) { cat("   Program running with ",getDoParWorkers()," workers/threads.\n") }
  #}}}

  #Initialise Loop Counter {{{
  results<-{}
  nloops<-length(param.env$datamap)
  if (!quiet) { cat("   There are ",nloops," files to analyse:\n") }
  #}}}

  #If doing multiple loops; save the parameters in case we need to reset mid-run {{{
  if ((nloops>1)&(!resume)) {
    save(file=".LambdarParameters.Rdata", param.env)
  }
  #}}}

  #Loop through files supplied {{{
  for (f in lstart:nloops) {
    #Set restart value {{{
    write.csv(file=".LambdarActiveLoop.txt", c(f), row.names=FALSE)
    #}}}

    #Initialise Timer and Get Parameters for this run {{{
    loopstarttime<-proc.time()[3]
    getNthVar(parameterList,n=f,inenv=param.env,outenv=environment(),lim=nloops)
    dir.create(file.path(pathroot,pathwork, pathout), showWarnings = FALSE)
    #}}}

    #Send Message output to logfile {{{
    sinkfile<-file(file.path(pathroot,pathwork,pathout,logfile),open="wt")
    sink(sinkfile, type="message")
    on.exit(sink(type="message"), add=TRUE)
    #Print any warnings {{{
    if ((!quiet)&(!(is.null(warnings())))) {
    cat('\n')
    print(warnings())
    cat('\n')
    }
    #}}}
    #}}}

    #If wanted, crop image prior to read {{{
    if (cropimage) {
      #Data Image {{{
      if (verbose) { message(paste("Cropping Input Image: Outputting to", imfitsoutname)) }
      crop_im(ra0=ra0, dec0=dec0, pathroot=file.path(pathroot,pathwork), inpim=datamap, cutrad=cutrad, fitsoutname=imfitsoutname)
      if (verbose) { message(paste("Using", imfitsoutname, "as data image")) }
      datamap<-imfitsoutname
      #}}}
    }#}}}

    #Create Seperate Images Environment {{{
    #Details {{{
    # Create the image environment, which will contain all
    # Image Data for the LAMBDAR Proceedure - This segregates the
    # large image arrays from the regular parameter space, and
    # stops unnecessary memory usage in foreach loops
    #}}}
    image.env<-new.env(parent=environment())
    #}}}

    #Read in Data, Mask map, & Error map {{{
    readimage(env=NULL,quiet,showtime,outenv=image.env)
    #Move Astrometry list from image env to main env {{{
    astr_struc<-image.env$astr_struc
    rm(astr_struc, envir=image.env)
    #}}}
    #}}}

    #Read source catalogue {{{
    opencatalogue(outenv=environment())
    #}}}

    #If wanted, Set Minimum Aperture Radius {{{
    if(MinApRad>0){
      a_g[a_g<MinApRad]=MinApRad
      b_g[b_g<MinApRad]=MinApRad
    }
    #}}}

    #Determine Initial Galaxy Sample {{{
    #Details {{{
    #Determine which galaxies are inside the input image,
    #and which galaxies have physical aperture parameters.
    #If either of these are not true, then the aperture is
    #discarded.}}}
    if (!quiet) { cat("   Determining Correct Galaxy Sample ") }

    #Get object locations in pixel space {{{
    gamapos<-ad2xy(ra_g,dec_g,astr_struc)
    x_g<-gamapos[,1]
    y_g<-gamapos[,2]
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
      message(paste("X_G:", length(x_g),"\nY_G:",length(y_g)))
      message(paste("min/max X_G:", min(x_g), max(x_g),"\nmin/max Y_G:",min(x_g), max(y_g)))
    }#}}}

    #Discard any apertures that lie completely outside of the image ± 1 pixel {{{
    catlen<-length(x_g)
    insidemask<-!((x_g <= 0) | (x_g >= length(image.env$im[,1])+1) | (y_g <= 0) | (y_g >= length(image.env$im[1,])+1))
    #Check that something is inside the image {{{
    if (length(which(insidemask==TRUE))==0) {
      warning("No Single Apertures are inside the image.")
      #Notify & Close Logfile {{{
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      message(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loopstarttime, digits=3),'\n'))
      message(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n'))
      sink(type='message')
      close(sinkfile)
      if (!quiet) {
        cat(paste('- Done\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
        cat(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loopstarttime, digits=3),'\n'))
        cat(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n'))
        if (f!=nloops) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
      }
      next
      #}}}
    }
    #}}}
    #Remove object catalogue entries {{{
    x_g<-x_g[which(insidemask)]
    y_g<-y_g[which(insidemask)]
    id_g<-id_g[which(insidemask)]
    ra_g<-ra_g[which(insidemask)]
    dec_g<-dec_g[which(insidemask)]
    theta_g<-theta_g[which(insidemask)]
    a_g<-a_g[which(insidemask)]
    b_g<-b_g[which(insidemask)]
    if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
    if (filtcontam) { contams<-contams[which(insidemask)] }
    #}}}
    #Notify how many objects remain {{{
    if (verbose) { message(paste("There are",length(x_g),"supplied objects inside the image (",
                                  round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied were outside the image )")) }
    #}}}
    #}}}

    #Discard any apertures that have nonphysical aperture axis values {{{
    catlen<-length(x_g)
    insidemask<-!((a_g < 0)|(b_g < 0))
    #Check that something is inside the image {{{
    if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Apertures remaining have physical axis values.") }
    #}}}
    #Remove object catalogue entries {{{
    x_g<-x_g[which(insidemask)]
    y_g<-y_g[which(insidemask)]
    id_g<-id_g[which(insidemask)]
    ra_g<-ra_g[which(insidemask)]
    dec_g<-dec_g[which(insidemask)]
    theta_g<-theta_g[which(insidemask)]
    a_g<-a_g[which(insidemask)]
    b_g<-b_g[which(insidemask)]
    if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
    if (filtcontam) { contams<-contams[which(insidemask)] }
    chunkSize=length(id_g)/getDoParWorkers()
    mpiopts<-list(chunkSize=chunkSize)
    #}}}
    #Notify how many objects remain {{{
    if (verbose) { message(paste("There are",length(x_g),"supplied objects with physical aperture values (",
                                  round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied had unphysical values )")) }
    #}}}
    #}}}

    if (!quiet) { cat(" - Done\n") }
    #}}}

    #Set object Astrometry {{{
    if (!quiet) { cat("   Setting Astrometry   ") }

    #Get pixel resolution from astrometry {{{
    if (all(is.finite(astr_struc$CDELT))) {
      #Using CDELT keywords {{{
      asperpix<-max(abs(astr_struc$CDELT))*3600.
      #}}}
    } else if (all(is.finite(astr_struc$CD[1:2,1:2]))) {
      #Using CD matrix keywords {{{
      asperpix<-max(astr_struc$CD[1,1],astr_struc$CD[2,2])*3600.
      #}}}
    } else {
      #Otherwise; Error - unknown Astrometry Resolution Keyword {{{
      sink(type=c("output","message")) ; stop("Data image header does not contain CD or CDELT keywords")
      #}}}
    }#}}}

    #Set apertures with NA/NULL aperture axis or minoraxis<aperturediag to point-sources {{{
    #Details {{{
    #If an aperture has no provided aperture parameters
    #or if the aperture width is smaller than the pixel
    #diagonal, the object should be treated as a point source.
    #This is because if an aperture is too small/thin, then
    #positional information of the aperture inside the pixel
    #is completely lost, and convolution with a PSF
    #will be inaccurate because it cannot be accurately
    #convolved. }}}
    #Get pixel diagonal size in arcsec {{{
    diag_arcsec<-abs(asperpix)*sqrt(2)/2
    #}}}
    #Make needed changes, and notify {{{
    message("Forcing ",length(which((b_g<diag_arcsec)|!is.finite(a_g)))," apertures to be point sources")
    a_g[which((b_g<diag_arcsec)|!is.finite(a_g))]<-0
    theta_g[which((b_g<diag_arcsec)|!is.finite(theta_g))]<-0
    b_g[which((b_g<diag_arcsec)|!is.finite(b_g))]<-0
    #}}}
    #}}}

    if (!quiet) { cat(" - Done\n") }
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
      message(paste("X_G:", length(x_g),"\nY_G:",length(y_g)))
      message(paste("min/max X_G:", min(x_g), max(x_g),"\nmin/max Y_G:",min(x_g), max(y_g)))
    }#}}}

    #Finished setting astrometry. Notify {{{
    if (!quiet) {
      cat("   } - Done\n")
      cat('} Initialisation Complete ')
    }
    if (showtime) { cat(paste(' (  Loop Time Elapsed (s): ',round(proc.time()[3]-loopstarttime, digits=3),'  )\n')) }
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
    message(paste('Arcsec per pixel in map: ', asperpix))
    message(paste('Beam area (from observers manual) converted into pixel units: ', beamarea_pix))
    }#}}}

    #Send Parameter to logfile {{{
    sink(sinkfile, type="output")
    cat("Parameters used in this run:\n")
    print(lsos(envir=environment(), head=FALSE))
    cat("Images used in this run:\n")
    print(lsos(envir=image.env, head=FALSE))
    sink(type="output")
    #}}}

    #Create Simulated Image from Catalogue {{{
    simcat<-create.sim.image(ObsParm=ObsParm,noNoise=noNoise,convolveNoise=convolveNoise,padGals=padGalaxies,col.corr=colourCorr,confuse=confuse)
    #}}}

    #Notify & Close Logfile {{{
    message(paste('-----------------------------------------------------\nDatamap Complete\n'))
    message(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loopstarttime, digits=3),'\n'))
    message(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n'))
    sink(type='message')
    close(sinkfile)
    #}}}

    #Loop Completed - Print Update, Loop {{{
    if (!quiet) {
      cat(paste('-----------------------------------------------------\nDatamap Complete\n'))
      cat(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loopstarttime, digits=3),'\n'))
      cat(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-starttime, digits=3),'\n'))
      if (f!=nloops) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
    }
    #}}}
  }#}}}

  #Remove State Preservation Files {{{
  if (nloops>1) {
    system("rm -f .LambdarActiveLoop.txt .LambdarParameters.Rdata")
  }#}}}

  #Program Completed - Print closing {{{
  if (!quiet) {
    cat(paste('-----------------------------------------------------\nProgram Complete\n'))
  }
  return=NULL
  #}}}

}
