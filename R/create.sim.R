create.sim <-
function(par.file=NA, ObsParm=NA, injection=FALSE, n.inject=100, noNoise=FALSE, convolveNoise=TRUE, padGalaxies=TRUE, colourCorr=0.0, quiet=FALSE, confuse=FALSE){
#Proceedure measures object fluxes from an arbitrary fits image

  #For Setup, warnings are handled internally - print nothing {{{
  options(warn=-1)
  #}}}

  #Set function Environments {{{
  environment(open.catalogue)<-environment()
  environment(read.images)<-environment()
  environment(read.par.file)<-environment()
  environment(create.sim.image)<-environment()
  environment(create.injection.image)<-environment()
  for (nam in ls.deb("package:LAMBDAR",simple=TRUE)) { 
    if (nam%in%ls(envir=environment())) { 
      debug(get(nam,envir=environment()))
    }
  }
  #}}}

  #Check for appropriate calling syntax {{{
  if (is.na(par.file)) {
  stop(paste("Parameter file not supplied.\n",
             "Calling Syntax:\n",
             "       measure.fluxes(<ParameterFile>,<Observation Parameters>,<Quiet Flag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<Observation Parameters> = Data frame containing the details of the observation of the input image: \n",
             "            $exp = Effective Exposure Time (s)\n",
             "            $area = Telescope Collecting Area (m^2)\n",
             "            $lamEff = Filter Effective Central Wavelength (Ang)\n",
             "            $Weff = Filter Effective Width (Ang)\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run measure.fluxes('--makepar').", sep=""))
  }
  if (is.na(ObsParm)|!is.data.frame(ObsParm)) {
  warning(paste("Observation Parameters not supplied, or supplied incorrectly.\n",
             "Calling Syntax:\n",
             "       measure.fluxes(<ParameterFile>,<Observation Parameters>,<Quiet Flag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<Observation Parameters> = Data frame containing the details of the observation of the input image: \n",
             "            $exp = Effective Exposure Time (s)\n",
             "            $area = Telescope Collecting Area (m^2)\n",
             "            $lamEff = Filter Effective Central Wavelength (Ang)\n",
             "            $Weff = Filter Effective Width (Ang)\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run measure.fluxes('--makepar').", sep=""))
  warning("Using SDSSIII r-band default parameters for observation paramters\n")
  ObsParm=data.frame(exp=53.9,area=pi*(2.5/2)^2,lamEff=1111.2,Weff=6165)
  }#}}}

  #If requested, resume LAMBDAR run from last loop state {{{
  if (par.file == "--resume") {
    if (file.exists(".LambdarParameters.Rdata")&file.exists(".LambdarActiveLoop.txt")) {
      if (!quiet) { cat("Resuming from Previous Loop State\n") }
      load(".LambdarParameters.Rdata")
      loop.start<-as.numeric(read.csv(".LambdarActiveLoop.txt"))[1]
      resume<-TRUE
    } else {
      stop("Resume Requested, but one/both of the required resume files .LambdarParameters.Rdata & .LambdarLoopState.txt are missing")
    }
  } else {
    loop.start<-1
    resume<-FALSE
  }#}}}

  #Set start timer & print opening {{{
  start.time<-proc.time()[3]
  #}}}

  mem.lim<-Inf

  #Setup Parameter Space (read .par file) {{{
  if (!resume) {
    param.env<-new.env(parent=environment())
    read.par.file(par.file,start.time,quiet,env=param.env)
  }
  parameter.list<-ls(envir=param.env)
  #}}}

  #From here on, produce warnings as they occur {{{
  options(warn=1)
  #}}}

  #Check if magnitudes are active {{{
  if (any(!param.env$magnitudes) | any(is.na(param.env$mag.zp))) {
    stop("magnitudes must be TRUE for sim-image generation; mag.zp is used in all cases to output image as Jy/pix")
  }
  #}}}

  #If needed, register the parallel backend {{{
  registerDoParallel(cores=param.env$num.cores)
  if (!quiet) { cat("   Program running with ",getDoParWorkers()," workers/threads.\n") }
  #}}}

  #Initialise Loop Counter {{{
  results<-{}
  loop.total<-length(param.env$data.map)
  if (!quiet) { cat("   There are ",loop.total," files to analyse:\n") }
  #}}}

  #If doing multiple loops; save the parameters in case we need to reset mid-run {{{
  if ((loop.total>1)&(!resume)) {
    save(file=".LambdarParameters.Rdata", param.env)
  }
  #}}}

  #Loop through files supplied {{{
  for (f in loop.start:loop.total) {
    #Set restart value {{{
    write.csv(file=".LambdarActiveLoop.txt", c(f), row.names=FALSE)
    #}}}

    #Initialise Timer and Get Parameters for this run {{{
    loop.start.time<-proc.time()[3]
    get.nth.var(parameter.list,n=f,inenv=param.env,outenv=environment(),lim=loop.total)
    dir.create(file.path(path.root,path.work, path.out), showWarnings = FALSE)
    #}}}

    #Send Message output to logfile {{{
    sink.file<-file(file.path(path.root,path.work,path.out,logfile),open="wt")
    sink(sink.file, type="message")
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
    if (crop.image) {
      #Data Image {{{
      if (verbose) { message(paste("Cropping Input Image: Outputting to", data.fits.output.filename)) }
      crop.fits.image(ra0=ra0, dec0=dec0, path.root=file.path(path.root,path.work), inpim=data.map, crop.radius=crop.radius, fitsoutname=data.fits.output.filename)
      if (verbose) { message(paste("Using", data.fits.output.filename, "as data image")) }
      data.map<-data.fits.output.filename
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
    read.images(env=NULL,quiet,showtime,outenv=image.env)
    #Move Astrometry list from image env to main env {{{
    astr.struc<-image.env$astr.struc
    saturation<-image.env$saturation
    gain<-image.env$gain
    rm(astr.struc, envir=image.env)
    #}}}
    #}}}

    #Read source catalogue {{{
    open.catalogue(outenv=environment())
    #}}}

    #If wanted, Set Minimum Aperture Radius {{{
    if(min.ap.rad>0){
      cat.a[cat.a<min.ap.rad]=min.ap.rad
      cat.b[cat.b<min.ap.rad]=min.ap.rad
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
    gama.pos<-ad.to.xy(cat.ra,cat.dec,astr.struc)
    cat.x<-gama.pos[,1]
    cat.y<-gama.pos[,2]
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
      message(paste("X.pix:", length(cat.x),"\nY.pix:",length(cat.y)))
      message(paste("min/max X.pix:", min(cat.x), max(cat.x),"\nmin/max Y.pix:",min(cat.x), max(cat.y)))
    }#}}}

    #Discard any apertures that lie completely outside of the image ± 1 pixel {{{
    cat.len<-length(cat.x)
    inside.mask<-!((cat.x <= 0) | (cat.x >= length(image.env$im[,1])+1) | (cat.y <= 0) | (cat.y >= length(image.env$im[1,])+1))
    #Check that something is inside the image {{{
    if (length(which(inside.mask==TRUE))==0) {
      warning("No Single Apertures are inside the image.")
      #Notify & Close Logfile {{{
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      message(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'\n'))
      message(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-start.time, digits=3),'\n'))
      sink(type='message')
      close(sink.file)
      if (!quiet) {
        cat(paste('- Done\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
        cat(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'\n'))
        cat(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-start.time, digits=3),'\n'))
        if (f!=loop.total) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
      }
      next
      #}}}
    }
    #}}}
    #Remove object catalogue entries {{{
    cat.x<-cat.x[which(inside.mask)]
    cat.y<-cat.y[which(inside.mask)]
    cat.id<-cat.id[which(inside.mask)]
    cat.ra<-cat.ra[which(inside.mask)]
    cat.dec<-cat.dec[which(inside.mask)]
    cat.theta<-cat.theta[which(inside.mask)]
    cat.a<-cat.a[which(inside.mask)]
    cat.b<-cat.b[which(inside.mask)]
    if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
    if (filt.contam) { contams<-contams[which(inside.mask)] }
    #}}}
    #Notify how many objects remain {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects inside the image (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of supplied were outside the image )")) }
    #}}}
    #}}}

    #Discard any apertures that have nonphysical aperture axis values {{{
    cat.len<-length(cat.x)
    inside.mask<-!((cat.a < 0)|(cat.b < 0))
    #Check that something is inside the image {{{
    if (length(which(inside.mask==TRUE))==0) { sink(type="message") ; stop("No Apertures remaining have physical axis values.") }
    #}}}
    #Remove object catalogue entries {{{
    cat.x<-cat.x[which(inside.mask)]
    cat.y<-cat.y[which(inside.mask)]
    cat.id<-cat.id[which(inside.mask)]
    cat.ra<-cat.ra[which(inside.mask)]
    cat.dec<-cat.dec[which(inside.mask)]
    cat.theta<-cat.theta[which(inside.mask)]
    cat.a<-cat.a[which(inside.mask)]
    cat.b<-cat.b[which(inside.mask)]
    if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
    if (filt.contam) { contams<-contams[which(inside.mask)] }
    chunk.size=length(cat.id)/getDoParWorkers()
    mpi.opts<-list(chunkSize=chunk.size)
    message("Number of objects per thread:",chunk.size)
    #}}}
    #Notify how many objects remain {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects with physical aperture values (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of supplied had unphysical values )")) }
    #}}}
    #}}}

    if (!quiet) { cat(" - Done\n") }
    #}}}

    #Set object Astrometry {{{
    if (!quiet) { cat("   Setting Astrometry   ") }

    #Get pixel resolution from astrometry {{{
    if (all(is.finite(astr.struc$CDELT))) {
      #Using CDELT keywords {{{
      arcsec.per.pix<-max(abs(astr.struc$CDELT))*3600.
      #}}}
    } else if (all(is.finite(astr.struc$CD[1:2,1:2]))) {
      #Using CD matrix keywords {{{
      arcsec.per.pix<-max(astr.struc$CD[1,1],astr.struc$CD[2,2])*3600.
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
    diag.arcsec<-abs(arcsec.per.pix)*sqrt(2)/2
    #}}}
    #Make needed changes, and notify {{{
    message("Forcing ",length(which((cat.b<diag.arcsec)|!is.finite(cat.a)))," apertures to be point sources")
    cat.a[which((cat.b<diag.arcsec)|!is.finite(cat.a))]<-0
    cat.theta[which((cat.b<diag.arcsec)|!is.finite(cat.theta))]<-0
    cat.b[which((cat.b<diag.arcsec)|!is.finite(cat.b))]<-0
    #}}}
    #}}}

    if (!quiet) { cat(" - Done\n") }
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
      message(paste("X.pix:", length(cat.x),"\nY.pix:",length(cat.y)))
      message(paste("min/max X.pix:", min(cat.x), max(cat.x),"\nmin/max Y.pix:",min(cat.x), max(cat.y)))
    }#}}}

    #Finished setting astrometry. Notify {{{
    if (!quiet) {
      cat("   } - Done\n")
      cat('} Initialisation Complete ')
    }
    if (showtime) { cat(paste(' (  Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'  )\n')) }
    #}}}

    #-----Diagnostic-----# {{{
    if (diagnostic) {
    message(paste('Arcsec per pixel in map: ', arcsec.per.pix))
    message(paste('Beam area (from observers manual) converted into pixel units: ', beam.area.pix))
    }#}}}

    #Send Parameter to logfile /*fold*/ {{{
    sink(sink.file, type="output")
    cat("Parameters used in this run:\n")
    print(lsos(envir=environment(), head=FALSE))
    cat("Images used in this run:\n")
    print(lsos(envir=image.env, head=FALSE))
    sink(type="output")
    #/*fend*/ }}}

    #Create Simulated Image from Catalogue {{{
    if (injection) { 
      simcat<-create.injection.image(ObsParm=ObsParm,col.corr=colourCorr)
    } else { 
      simcat<-create.sim.image(ObsParm=ObsParm,noNoise=noNoise,convolveNoise=convolveNoise,padGals=padGalaxies,col.corr=colourCorr,confuse=confuse)
    } 
    #}}}

    #Notify & Close Logfile {{{
    message(paste('-----------------------------------------------------\nDatamap Complete\n'))
    message(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'\n'))
    message(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-start.time, digits=3),'\n'))
    sink(type='message')
    close(sink.file)
    #}}}

    #Loop Completed - Print Update, Loop {{{
    if (!quiet) {
      cat(paste('-----------------------------------------------------\nDatamap Complete\n'))
      cat(paste('Loop Completed: Indiv. Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'\n'))
      cat(paste('                      Total Time Elapsed (s): ',round(proc.time()[3]-start.time, digits=3),'\n'))
      if (f!=loop.total) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
    }
    #}}}
  }#}}}

  #Remove State Preservation Files {{{
  if (loop.total>1) {
    system("rm -f .LambdarActiveLoop.txt .LambdarParameters.Rdata")
  }#}}}

  #Program Completed - Print closing {{{
  if (!quiet) {
    cat(paste('-----------------------------------------------------\nProgram Complete\n'))
  }
  return=NULL
  #}}}

}
