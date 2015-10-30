MeasureFluxes <-
function(parfile=NA, quiet=FALSE, MPIBackend=FALSE, doReturn=FALSE, ...){
#Proceedure measures object fluxes from an arbitrary fits image

  #If needed, register the parallel backend {{{
  if (MPIBackend) {
    #Check for MPI libraries
    if (!(any(grepl("doMPI", search()))&&any(grepl("Rmpi", search())))) {
      #MPI libraries are *not* already loaded
      loadlibs<-try(require("Rmpi"),silent=TRUE)
      if (class(loadlibs)=="try-error") {
        stop(cat("MPI backend requested, but Rmpi package could not be loaded.",
                 "\nPackage may be in a special library location not known to LAMBDAR.",
                 "\nTry loading it before running MeasureFluxes"))
      }
      loadlibs<-try(require("doMPI"),silent=TRUE)
      if (class(loadlibs)=="try-error") {
        stop(cat("MPI backend requested, but doMPI package could not be loaded.",
                 "\nPackage may be in a special library location not known to LAMBDAR.",
                 "\nTry loading it before running MeasureFluxes"))
      }
    }
    #Create Cluster
    #message("The Size of the MPI Universe is ",mpi.universe.size()," with Comm Size ",mpi.comm.size(0))
    suppressMessages(cl<-startMPIcluster())
    registerDoMPI(cl)
    on.exit(closeCluster(cl))
    on.exit(mpi.quit(),add=TRUE)
    if (!quiet) { cat("The Size of the MPI Universe is ",mpi.universe.size()," with Comm Size ",mpi.comm.size(0)) }
  }
  #}}}

  #For Setup, warnings are handled internally - print nothing {{{
  options(warn=-1)
  #}}}

  #Set function Environments {{{
  environment(opencatalogue)<-environment()
  environment(readimage)<-environment()
  environment(readparfile)<-environment()
  environment(fluxmeasurements)<-environment()
  #}}}

  #Check for appropriate calling syntax {{{
  if (is.na(parfile)) {
  stop(paste("Parameter file not supplied.\n",
             "Calling Syntax:\n",
             "       MeasureFluxes(<ParameterFile>,<QuietFlag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run MeasureFluxes('--makepar').", sep=""))
  }#}}}

  #If requested, produce the default .par file and end {{{
  if (parfile == "--makepar") {
    if (!quiet) { cat("Outputting Default Parameter file to './Lambdar_default.par'\n") }
    createparfile(...)
    if (!quiet) { cat("Program Complete\n") }
    return()
  }#}}}

  #If parameter file is NULL, run program using specified parameters and (all other) default parameters file and end {{{
  if (parfile == "--run") {
    if (!quiet) { cat("Running LAMBDAR with parameters specified here only\n") }
    system(paste("echo -e '",paste(names(list(...)),list(...),sep=' ',collapse='\n'),"' > .tmpLam.par"))
    parfile='.tmpLam.par'
  }#}}}

  #If requested, resume LAMBDAR run from last loop state {{{
  if (parfile == "--resume") {
    if (file.exists(".LambdarParameters.Rdata")&file.exists(".LambdarActiveLoop.txt")) {
      if (!quiet) { cat("Resuming from Previous Loop State ") }
      load(".LambdarParameters.Rdata")
      lstart<-as.numeric(read.csv(".LambdarActiveLoop.txt"))[1]
      if (!quiet) { cat(paste0("(loop ",lstart,")\n")) }
      resume<-TRUE
    } else {
      stop("Resume Requested, but one/both of the required resume files .LambdarParameters.Rdata & .LambdarActiveLoop.txt are missing")
    }
  } else {
    lstart<-1
    resume<-FALSE
  }#}}}

  #Set start timer & print opening {{{
  starttime<-proc.time()[3]
  #}}}

  #Get System Memory Limit (platform dependant) before any assignments {{{
  #Can't do this when using MPI backend, as system calls can cause data corruption
  if (!MPIBackend) {
    sysName=Sys.info()[1]
    if (sysName=="Linux") {
      #Linux Systems {{{
      #Memory Limit returned in kBytes. Convert to bits.
      memTot<-as.numeric(system("awk '/MemTotal:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
      memAct<-as.numeric(system("awk '/Active:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
      memLim<-memTot-memAct
      if (!is.finite(memLim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        memLim<-Inf
      }
      #}}}
    } else if (sysName=="Darwin") {
      #Mac Systems {{{
      memLim<-system("top -l 1 | grep PhysMem | awk '{print $6}'", intern=TRUE)
      #Determine unit and convert to Bits
           if (grepl('G',memLim)) { memLim<-as.numeric(strsplit(memLim,'G'))*1E9*8 } #Gigabytes
      else if (grepl('M',memLim)) { memLim<-as.numeric(strsplit(memLim,'M'))*1E6*8 } #Megabytes
      else if (grepl('k',memLim)) { memLim<-as.numeric(strsplit(memLim,'k'))*1E3*8 } #Kilobytes
      else                        { memLim<-as.numeric(memLim)*8 }                   #bytes
      if (!is.finite(memLim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        memLim<-Inf
      }
      #}}}
    } else if (sysName=="Windows") {
      #Windows Machines {{{
      #Memory Limit returned in Bytes. Convert to Bits
      memLim<-as.numeric(memory.limit())*8
      if (!is.finite(memLim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        memLim<-Inf
      }
      #}}}
    } else {
      #Any Others, warn & set to Inf {{{
      warning("Unknown Operating System Name. Cannot determine Memory Limits.\nUsing Infinity")
      message("Unknown Operating System Name. Cannot determine Memory Limits.\nUsing Infinity")
      memLim<-Inf
      #}}}
    }
  }#}}}

  #Setup Parameter Space (read .par file) {{{
  if (!resume) {
    param.env<-new.env(parent=environment())
    parwarning<-readparfile(parfile,starttime,quiet,env=param.env)
  }
  parameterList<-ls(envir=param.env)
  #}}}

  #From here on, produce warnings as they occur {{{
  options(warn=1)
  #}}}

  #If needed, register the parallel backend {{{
  if (!MPIBackend) {
    #Straight Register of cores
    registerDoParallel(cores=param.env$ncores)
    if (!quiet) { cat("   Program running with ",getDoParWorkers()," workers/threads.\n") }
  }
  #}}}

  #Initialise Loop Counter {{{
  results<-{}
  nloops<-max(sapply(ls(envir=param.env), function(x) length(param.env[[x]])))
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
    sinkfile<-try(file(file.path(pathroot,pathwork,pathout,logfile),open="wt"),silent=TRUE)
    if (any(class(sinkfile)=='try-error')) {
      stop("Unable to open Log File. Likely the path is incorrect or you do not have permission to write there:\n   ",file.path(pathroot,pathwork,pathout,logfile))
    }
    sink(sinkfile, type="message")
    on.exit(sink(type="message"), add=TRUE)
    #Print any warnings {{{
    if (!exists("parwarning")) { parwarning<-NULL }
    if ((!(is.null(parwarning)))) {
      message("Warnings in Parameter File Read:")
      message(parwarning)
    }
    #}}}
    #}}}

    #If wanted, crop image prior to read {{{
    if (cropimage) {
      #Data Image {{{
      if (datamap!=imfitsoutname) {
        if (verbose) { message(paste("Cropping Input Image: Outputting to", imfitsoutname)) }
        if (!file.exists(file.path(pathroot,pathwork,datamap))) { sink(type='message'); stop("Data Image does not exist at location:",file.path(pathroot,pathwork,datamap)) }
        if (extn>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond extn 1)") }
        crop_im(ra0=ra0, dec0=dec0, pathroot=file.path(pathroot,pathwork), inpim=datamap, cutrad=cutrad, fitsoutname=imfitsoutname,extn=1)
        if (verbose) { message(paste("Using", imfitsoutname, "as data image")) }
        datamap<-imfitsoutname
      } else {
        message("Crop input and output are the same file. Crop will fail, so it is skipped")
        warning("Crop input and output are the same file. Crop will fail, so it is skipped")
      }
      #If RA and/or DEC are -999, get RA and DEC of the data cutout {{{
      if (ra0==-999 | dec0==-999) {
        astrtmp<-read.astr(file.path(pathroot,pathwork,datamap),hdu=extn)
        pos<-xy2ad(astrtmp$NAXIS[1]/2,astrtmp$NAXIS[2]/2,astrtmp)
        if (ra0==-999) {
          ra0<-pos[1]
          message(paste("Setting RA0 to centre of datamap:",round(ra0,4)))
        }
        if (dec0==-999) {
          dec0<-pos[2]
          message(paste("Setting DEC0 to centre of datamap:",round(dec0,4)))
        }
      }
      #}}}
      #}}}
      #Mask Image {{{
      if (maskmap != "NONE") {
        if (maskmap!=immfitsoutname) {
          if (verbose) { message(paste("Cropping Input Mask Map: Outputting to", immfitsoutname)) }
          if (!file.exists(file.path(pathroot,pathwork,maskmap))) { sink(type='message'); stop("Mask Image does not exist at location:",file.path(pathroot,pathwork,maskmap)) }
          if (extnmask>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond extn 1)") }
          crop_im(ra0=ra0, dec0=dec0, pathroot=file.path(pathroot,pathwork), inpim=maskmap, cutrad=cutrad, fitsoutname=immfitsoutname,extn=1)
          if (verbose) { message(paste("Using", immfitsoutname, "as mask map")) }
          maskmap<-immfitsoutname
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }
      #}}}
      #Or Weights Image {{{
      if (wgtmap != "NONE") {
        if (!file.exists(file.path(pathroot,pathwork,wgtmap))) {
          wgtmap<-"NONE"
        } else if (wgtmap!=imwgtfitsoutname) {
          if (verbose) { message(paste("Cropping Input Weight Map: Outputting to", imwgtfitsoutname)) }
          if (!file.exists(file.path(pathroot,pathwork,wgtmap))) { sink(type='message'); stop("Weight Image does not exist at location:",file.path(pathroot,pathwork,wgtmap)) }
          if (extnwgt>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond extn 1)") }
          crop_im(ra0=ra0, dec0=dec0, pathroot=file.path(pathroot,pathwork), inpim=wgtmap, cutrad=cutrad, fitsoutname=imwgtfitsoutname,extn=1)
          if (verbose) { message(paste("Using", imwgtfitsoutname, "as weight map")) }
          wgtmap<-imwgtfitsoutname
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }#}}}
      #Error Image {{{
      if ((errormap != "NONE")&(is.na(as.numeric(errormap)))) {
        if (errormap!=imefitsoutname) {
          if (verbose) { message(paste("Cropping Input Error Map: Outputting to", imefitsoutname)) }
          if (!file.exists(file.path(pathroot,pathwork,errormap))) { sink(type='message'); stop("Error Image does not exist at location:",file.path(pathroot,pathwork,errormap)) }
          if (extnerr>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond extn 1)") }
          crop_im(ra0=ra0, dec0=dec0, pathroot=file.path(pathroot,pathwork), inpim=errormap, cutrad=cutrad, fitsoutname=imefitsoutname,extn=1)
          if (verbose) { message(paste("Using", imefitsoutname, "as error map")) }
          errormap<-imefitsoutname
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }#}}}
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
    saturation<-image.env$saturation
    rm(astr_struc, envir=image.env)
    #}}}
    #}}}

    #If needed, read ZP Magnitude from Image header {{{
    if ((Magnitudes) & (magZP==-999)){
      magZP<-as.numeric(read.fitskey(magZPlabel,paste(pathroot,pathwork,datamap,sep=""),hdu=extn))
      #If Failed, do not output Magnitudes {{{
      if (!is.finite(magZP)) {
        message("Zero Point Magnitude determination failed - Not outputting Magnitudes")
        warning("Zero Point Magnitude determination failed")
        Magnitudes<-FALSE
      }#}}}
    }#}}}

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
      warning("No catalogue entries are centred within the limits of the input image.")
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
    #if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Single Apertures are inside the image.") }
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
    insidemask<-insidemask[which(insidemask)]
    #}}}
    if (filtcontam) {
      insidemask<-insidemask & !contams
      #Check that there is a science target inside the image {{{
      if (length(which(insidemask==TRUE))==0) {
        warning("No science targets are centred within the limits of the input image.")
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
          if (f!=nloops) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap (',f+1,' of ',nloops,')\nInitialising Workspace {\n',sep="")) }
        }
        next
        #}}}
      }
      #}}}
    }
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
    message("Number of objects per thread:",chunkSize)
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
    if (all(is.finite(astr_struc$CDELT[1:2]))) {
      #Using CDELT keywords {{{
      asperpix<-max(abs(astr_struc$CDELT[1:2]))*3600.
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

    #If wanted, check memory-safe {{{
    if ((!MPIBackend)&&(memSafe & is.finite(memLim))) {
      #Check that computation is able to be performed within memory limits {{{
      if (!quiet) { cat("   Checking Memory Usage & Limits {") }
      #Details {{{
      #Perform a rudimentary check that the image isn't too
      #large for the Available system memory, using
      # approx aperture memory:
      #     apmem = napertures*apsizeinBits*nLists
      # where
      #     nLists=12
      #     apsizeinbits = (90th percentile{(semimajor.axis)}*2/asperpix)^2*bitsperpixel
      #
      # approx image memory during calculations:
      #     immem = nimages*imsizeinBytes*nThreads
      #
      # current memory usage from parameter/image storage
      #     currmem = memory available at initialisation - memory currently used.
      #}}}
      #Get System Memory Usage (platform dependant) at Current Time {{{
      if (sysName=="Linux") {
        #Linux Systems {{{
        #Memory returned in kBytes. Convert to bits.
        memTot<-as.numeric(system("awk '/MemTotal:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
        memAct<-as.numeric(system("awk '/Active:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
        memCur<-memTot-memAct
        if ((length(memCur)==0)||!is.finite(memCur)) {
          warning("Memory determination failed. Setting to Inf.")
          memCur<-Inf
        }
        #}}}
      } else if (sysName=="Darwin") {
        #Mac Systems {{{
        memCur<-system("top -l 1 | grep PhysMem | awk '{print $6}'", intern=TRUE)
        #Determine unit and convert to Bits
             if (grepl('G',memCur)) { memCur<-as.numeric(strsplit(memCur,'G'))*1E9*8 } #Gigabytes
        else if (grepl('M',memCur)) { memCur<-as.numeric(strsplit(memCur,'M'))*1E6*8 } #Megabytes
        else if (grepl('k',memCur)) { memCur<-as.numeric(strsplit(memCur,'k'))*1E3*8 } #Kilobytes
        else                        { memCur<-as.numeric(memCur)*8 }                   #bytes
        if (!is.finite(memCur)) {
          warning("Memory determination failed. Setting to Inf.")
          memCur<-Inf
        }
        #}}}
      } else if (sysName=="Windows") {
        #Windows Machines {{{
        #Memory returned in Bytes. Convert to Bits
        memCur<-as.numeric(memory.limit())*8
        if (!is.finite(memCur)) {
          warning("Memory determination failed. Setting to Inf.")
          memCur<-Inf
        }
        #}}}
      } else {
        #Any Others, warn & set to Inf {{{
        warning("Unknown Operating System Name. Cannot determine Memory Usages.\nUsing Infinity")
        message("Unknown Operating System Name. Cannot determine Memory Usages.\nUsing Infinity")
        memCur<-Inf
        #}}}
      }
      memCur<-memLim-memCur
      #}}}
      #Aperture Memory requirements {{{
      catlen<-length(id_g)
      #Use 90th Quantile of aperture semimajor axes {{{
      aprad.quant<-quantile(a_g[which(a_g>0)],0.9)
      if (is.na(aprad.quant)) {
        #All apertures are point sources - make a default width @ 10 pix
        aprad.quant<-10
      }
      #}}}
      apsizeinbits<-(2*aprad.quant/asperpix)^2*64*2
      apmem<-catlen*apsizeinbits*12
      #}}}
      #Image Memory requirements {{{
      nimage=4
      if (length(image.env$imm)>1){ nimage=nimage+1 }
      if (length(image.env$ime)>1){ nimage=nimage+1 }
      if (sourcemask){ nimage=nimage+1 }
      immem<-nimage*lsos(pattern="im",envir=image.env)[1,'Size']
      #}}}
      #Check Memory Allocation is less than available free memory {{{
      if ((apmem+immem+memCur) >= memLim) {
          #If this is too much, less threads cannot help. Stop {{{
          cat("Failed\n")
          sink(type='message')
          stop(paste("This computation may exceed the available memory on this machine.",
                     "The required memory (",(apmem+immem+memCur)*1E-9/8,"Gb) is greater than that which is available to the system (",(memLim)*1E-9/8,"Gb).",
                     "\nHowever, using the in-built crop function, seperating the image into smaller chuncks will enable computation.",
                     "\nThe memory usage is",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                     "\n",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n"))
          #}}}
        ##If not, try and continue by reducing the number of threads {{{
        #warning("Memory Required for this computation could possibly exceed system RAM\n")
        #cat(paste("\n        Memory Required for this computation could possibly exceed system RAM\n",
        #             "       Attempting to limiti memory usage by lowering number of threads - "))
        ##Memory Required on 1 thread {{{
        #immem<-3*lsos(pattern="im",envir=image.env)[1,'Size']
        ##}}}
        #if ((apmem+immem+memCur) >= memLim) {
        #  #If this is too much, less threads cannot help. Stop {{{
        #  cat("Failed\n")
        #  sink(type='message')
        #  stop(paste("This computation is not possible on this machine,",
        #             "as the required memory (",(apmem+immem+memCur)*1E-9/8,"Gb) is greater than that which is available to the system (",(memLim)*1E-9/8,"Gb).",
        #             "\nHowever, using the in-built crop function, seperating the image into smaller chuncks will enable computation.",
        #             "\nThe memory usage is",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
        #             "\n",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n"))
        #  #}}}
        #} else {
        ##Otherwise, determine the maximum number of threads that won't fail {{{
        #  for( i in 2:ncores) {
        #    immem<-3*lsos(pattern="im",envir=image.env)[1,'Size']*i
        #    if ((apmem+immem) >= memLim) {
        #      ncores<-i-1
        #      immem<-3*lsos(pattern="im",envir=image.env)[1,'Size']*ncores
        #      break
        #    }
        #  }
        #  #}}}
        #  #Update number of threads to be used, and notify {{{
        #  cat(paste("Success. Continuing computations with",ncores,"threads\n"))
        #  registerDoParallel(cores=ncores)
        #  #}}}
        #}#}}}
      }#}}}
      #Notify Memory Usage Status {{{
      if (!quiet) {
        cat(paste("\n       This computation will require approximately ",round((apmem+immem)*1E-9/8,digits=3)," Gb of memory.\n",
                     "      The memory usage is ",round(apmem*1E-9/8,digits=3),"Gb for apertures/results/storage/etc, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                  "\n       ",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n",
                     "       Total memory available for the system is ",(memLim)*1E-9/8," Gb.\n",sep=""))
      }
      message(paste("This computation will require approximately",round((apmem+immem)*1E-9/8,digits=3),"Gb of memory.\n",
                     "      The memory usage is ",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                  "\n       ",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n",
                    "Total memory available for the system is",(memLim)*1E-9/8,"Gb."))
      #}}}
      #}}}
    }
    #}}}

    #Get beam area in pixels {{{
    beamarea_pix<-beamarea_SOM_as/(asperpix)^2.
    #}}}

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

    #Run Flux Measurements {{{
    loopresults<-fluxmeasurements()
    if (doReturn) {
      results<-c(results, loopresults)
    } else {
      loopresults<-NULL
    }
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
      if (f!=nloops) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap (',f+1,' of ',nloops,')\nInitialising Workspace {\n',sep="")) }
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
  if (doReturn) {
    return=results
  } else {
    return="There are no results being returned because doReturn=FALSE"
  }
  #}}}

}
