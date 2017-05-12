measure.fluxes <-
function(par.file=NA, quiet=FALSE, mpi.backend=FALSE, do.return=FALSE, stop.on.missing=TRUE, ...){
#Proceedure measures object fluxes from an arbitrary fits image

  #If needed, register the parallel backend /*fold*/ {{{
  if (mpi.backend) {
    #Check for MPI libraries
    if (!(any(grepl("doMPI", search()))&&any(grepl("Rmpi", search())))) {
      #MPI libraries are *not* already loaded
      loadlibs<-try(require("Rmpi"),silent=TRUE)
      if (class(loadlibs)=="try-error") {
        stop(cat("MPI backend requested, but Rmpi package could not be loaded.",
                 "\nPackage may be in a special library location not known to LAMBDAR.",
                 "\nTry loading it before running measure.fluxes"))
      }
      loadlibs<-try(require("doMPI"),silent=TRUE)
      if (class(loadlibs)=="try-error") {
        stop(cat("MPI backend requested, but doMPI package could not be loaded.",
                 "\nPackage may be in a special library location not known to LAMBDAR.",
                 "\nTry loading it before running measure.fluxes"))
      }
    }
    #Create Cluster
    suppressMessages(cl<-startMPIcluster())
    registerDoMPI(cl)
    on.exit(closeCluster(cl))
    on.exit(mpi.quit(),add=TRUE)
    if (!quiet) { cat("The Size of the MPI Universe is ",mpi.universe.size()," with Comm Size ",mpi.comm.size(0)) }
  }
  #/*fend*/ }}}

  #For Setup, warnings are handled internally - print nothing /*fold*/ {{{
  options(warn=-1)
  #/*fend*/ }}}

  #Set function Environments /*fold*/ {{{
  environment(open.catalogue)<-environment()
  environment(read.images)<-environment()
  environment(read.par.file)<-environment()
  environment(flux.measurements)<-environment()
  #/*fend*/ }}}

  #Check for appropriate calling syntax /*fold*/ {{{
  if (is.na(par.file)) {
  stop(paste("Parameter file not supplied.\n",
             "Calling Syntax:\n",
             "       measure.fluxes(<ParameterFile>,<QuietFlag>)\n",
             "<ParameterFile> = Path to, and Filename of, the LAMBDAR .par file\n",
             "<QuietFlag> = TRUE/FALSE\n\n",
             "To create the default parameter file, run measure.fluxes('--makepar').", sep=""))
  }#/*fend*/ }}}

  #If requested, produce the default .par file and end /*fold*/ {{{
  if (par.file == "--makepar") {
    if (!quiet) { cat("Outputting Default Parameter file to './Lambdar_default.par'\n") }
    create.par.file(...)
    if (!quiet) { cat("Program Complete\n") }
    return()
  }#/*fend*/ }}}

  #If parameter file is NULL, run program using specified parameters and (all other) default parameters file and end /*fold*/ {{{
  if (par.file == "--run") {
    if (!quiet) { cat("Running LAMBDAR with parameters specified here only\n") }
    system(paste("echo -e '",paste(names(list(...)),list(...),sep=' ',collapse='\n'),"' > .tmpLam.par"))
    par.file='.tmpLam.par'
  }#/*fend*/ }}}

  #If requested, resume LAMBDAR run from last loop state /*fold*/ {{{
  if (par.file == "--resume") {
    if (file.exists(".LambdarParameters.Rdata")&file.exists(".LambdarActiveLoop.txt")) {
      if (!quiet) { cat("Resuming from Previous Loop State ") }
      load(".LambdarParameters.Rdata")
      loop.start<-as.numeric(read.csv(".LambdarActiveLoop.txt"))[1]
      if (!quiet) { cat(paste0("(loop ",loop.start,")\n")) }
      resume<-TRUE
    } else {
      stop("Resume Requested, but one/both of the required resume files .LambdarParameters.Rdata & .LambdarActiveLoop.txt are missing")
    }
  } else {
    loop.start<-1
    resume<-FALSE
  }#/*fend*/ }}}

  #Set start timer & print opening /*fold*/ {{{
  start.time<-proc.time()[3]
  #/*fend*/ }}}

  #Get System Memory Limit (platform dependant) before any assignments /*fold*/ {{{
  #Can't do this when using MPI backend, as system calls can cause data corruption
  if (!mpi.backend) {
    sys.name=Sys.info()[1]
    if (sys.name=="Linux") {
      #Linux Systems /*fold*/ {{{
      #Memory Limit returned in kBytes. Convert to bits.
      memTot<-as.numeric(system("awk '/MemTotal:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
      memAct<-as.numeric(system("awk '/Active:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
      mem.lim<-memTot-memAct
      if (!is.finite(mem.lim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        mem.lim<-Inf
      }
      #/*fend*/ }}}
    } else if (sys.name=="Darwin") {
      #Mac Systems /*fold*/ {{{
      mem.lim<-system("top -l 1 | grep PhysMem | awk '{print $6}'", intern=TRUE)
      #Determine unit and convert to Bits
           if (grepl('G',mem.lim)) { mem.lim<-as.numeric(strsplit(mem.lim,'G'))*1E9*8 } #Gigabytes
      else if (grepl('M',mem.lim)) { mem.lim<-as.numeric(strsplit(mem.lim,'M'))*1E6*8 } #Megabytes
      else if (grepl('k',mem.lim)) { mem.lim<-as.numeric(strsplit(mem.lim,'k'))*1E3*8 } #Kilobytes
      else                        { mem.lim<-as.numeric(mem.lim)*8 }                   #bytes
      if (!is.finite(mem.lim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        mem.lim<-Inf
      }
      #/*fend*/ }}}
    } else if (sys.name=="Windows") {
      #Windows Machines /*fold*/ {{{
      #Memory Limit returned in Bytes. Convert to Bits
      mem.lim<-as.numeric(memory.limit())*8
      if (!is.finite(mem.lim)) {
        warning("Memory Limit determination failed. Setting to Inf.")
        mem.lim<-Inf
      }
      #/*fend*/ }}}
    } else {
      #Any Others, warn & set to Inf /*fold*/ {{{
      warning("Unknown Operating System Name. Cannot determine Memory Limits.\nUsing Infinity")
      message("Unknown Operating System Name. Cannot determine Memory Limits.\nUsing Infinity")
      mem.lim<-Inf
      #/*fend*/ }}}
    }
  }#/*fend*/ }}}

  #Setup Parameter Space (read .par file) /*fold*/ {{{
  if (!resume) {
    param.env<-new.env(parent=environment())
    param.warnings<-read.par.file(par.file,start.time,quiet,env=param.env)
  }
  parameter.list<-ls(envir=param.env)
  #/*fend*/ }}}

  #From here on, produce warnings as they occur /*fold*/ {{{
  options(warn=1)
  #/*fend*/ }}}

  #Initialise Loop Counter /*fold*/ {{{
  results<-{}
  loop.total<-max(sapply(ls(envir=param.env), function(x) length(param.env[[x]])))
  if (!quiet) { cat("   There are ",loop.total," files to analyse:\n") }
  #/*fend*/ }}}

  #If doing multiple loops; save the parameters in case we need to reset mid-run /*fold*/ {{{
  if ((loop.total>1)&(!resume)) {
    save(file=".LambdarParameters.Rdata", param.env)
  }
  #/*fend*/ }}}

  #If doing multiple loops; check if we can save time on table reads /*fold*/ {{{
  if (loop.total>1) {
    if (length(param.env$catalogue)==1) {
      #There is only one table; save time on every loop except the first
      reuse.table<-c(FALSE,rep(TRUE,loop.total-1))
    } else if (any(param.env$catalogue[-length(param.env$catalogue)]==param.env$catalogue[-1])) {
      #There are some loops where the table is the same as the one before it; save time on those loops (and not the first)
      reuse.table<-c(FALSE,param.env$catalogue[-1]==param.env$catalogue[-length(param.env$catalogue)])
    } else {
      #There are no sequential duplicate tables
      reuse.table<-rep(FALSE,loop.total)
    }
    if (resume) {
      reuse.table[loop.start]<-FALSE
    }
  } else {
    reuse.table<-FALSE
  }
  #/*fend*/ }}}

  #Loop through files supplied /*fold*/ {{{
  for (f in loop.start:loop.total) {
    #Set restart value /*fold*/ {{{
    write.csv(file=".LambdarActiveLoop.txt", c(f ), row.names=FALSE)
    #/*fend*/ }}}

    #Initialise Timer and Get Parameters for this run /*fold*/ {{{
    loop.start.time<-proc.time()[3]
    get.nth.var(parameter.list,n=f,inenv=param.env,outenv=environment(),lim=loop.total)
    dir.create(file.path(path.root,path.work, path.out), showWarnings = FALSE,recursive=TRUE)
    #/*fend*/ }}}

    #If needed, register the parallel backend /*fold*/ {{{
    if (!mpi.backend) {
      #Register Number of cores
      registerDoParallel(cores=num.cores)
      if (!quiet) { cat("   Program running with ",getDoParWorkers()," workers/threads.\n") }
    }
    #/*fend*/ }}}

    #Send Message output to logfile /*fold*/ {{{
    sink.file<-try(file(file.path(path.root,path.work,path.out,logfile),open="wt"),silent=TRUE)
    if (any(class(sink.file)=='try-error')) {
      stop("Unable to open Log File. Likely the path is incorrect or you do not have permission to write there:\n   ",file.path(path.root,path.work,path.out,logfile))
    }
    sink(sink.file, type="message")
    on.exit(sink(type="message"), add=TRUE)
    #Print Header in logfile /*fold*/ {{{
    message(paste('------------------------------------------------------\n'))
    message(paste('   LAMBDAR     version ',packageVersion("LAMBDAR"),': Log File\n'))
    message(paste('   created on ',Sys.time(),' by ',system("whoami",intern=TRUE)))
    message(paste('------------------------------------------------------\n'))
    #/*fend*/ }}}
    #Print any warnings /*fold*/ {{{
    if (!exists("param.warnings")) { param.warnings<-NULL }
    if ((!(is.null(param.warnings)))) {
      message("Warnings in Parameter File Read:")
      message(param.warnings)
    }
    #/*fend*/ }}}
    #/*fend*/ }}}

    #If wanted, crop image prior to read /*fold*/ {{{
    if (crop.image) {
      #Data Image /*fold*/ {{{
      if (data.map!=data.fits.output.filename) {
        if (verbose) { message(paste("Cropping Input Image: Outputting to", data.fits.output.filename)) }
        if (!file.exists(file.path(path.root,path.work,data.map))) { sink(type='message'); stop("Data Image does not exist at location:",file.path(path.root,path.work,data.map)) }
        if (data.extn>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond data.extn 1)") }
        crop.fits.image(ra0=ra0, dec0=dec0, path.root=file.path(path.root,path.work), inpim=data.map, crop.radius=crop.radius, fitsoutname=data.fits.output.filename,data.extn=1)
        if (verbose) { message(paste("Using", data.fits.output.filename, "as data image")) }
        data.map<-data.fits.output.filename
      } else {
        message("Crop input and output are the same file. Crop will fail, so it is skipped")
        warning("Crop input and output are the same file. Crop will fail, so it is skipped")
      }
      #If RA and/or DEC are -999, get RA and DEC of the data cutout /*fold*/ {{{
      if (ra0==-999 | dec0==-999) {
        astrtmp<-read.astrometry(file.path(path.root,path.work,data.map),hdu=1)
        pos<-xy.to.ad(astrtmp$NAXIS[1]/2,astrtmp$NAXIS[2]/2,astrtmp)
        if (ra0==-999) {
          ra0<-pos[1]
          message(paste("Setting RA0 to centre of data.map:",round(ra0,4)))
        }
        if (dec0==-999) {
          dec0<-pos[2]
          message(paste("Setting DEC0 to centre of data.map:",round(dec0,4)))
        }
      }
      #/*fend*/ }}}
      #/*fend*/ }}}
      #Mask Image /*fold*/ {{{
      if (mask.map != "NONE") {
        if (mask.map!=mask.fits.output.filename) {
          if (verbose) { message(paste("Cropping Input Mask Map: Outputting to", mask.fits.output.filename)) }
          if (!file.exists(file.path(path.root,path.work,mask.map))) { sink(type='message'); stop("Mask Image does not exist at location:",file.path(path.root,path.work,mask.map)) }
          if (data.mask.extn>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond data.extn 1)") }
          crop.fits.image(ra0=ra0, dec0=dec0, path.root=file.path(path.root,path.work), inpim=mask.map, crop.radius=crop.radius, fitsoutname=mask.fits.output.filename,data.extn=1)
          if (verbose) { message(paste("Using", mask.fits.output.filename, "as mask map")) }
          mask.map<-mask.fits.output.filename
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }
      #/*fend*/ }}}
      #Or Weights Image /*fold*/ {{{
      if (weight.map != "NONE") {
        if (!file.exists(file.path(path.root,path.work,weight.map))) {
          weight.map<-"NONE"
        } else if (weight.map!=weight.fits.output.filename) {
          if (verbose) { message(paste("Cropping Input Weight Map: Outputting to", weight.fits.output.filename)) }
          if (!file.exists(file.path(path.root,path.work,weight.map))) { sink(type='message'); stop("Weight Image does not exist at location:",file.path(path.root,path.work,weight.map)) }
          if (data.weight.extn>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond data.extn 1)") }
          crop.fits.image(ra0=ra0, dec0=dec0, path.root=file.path(path.root,path.work), inpim=weight.map, crop.radius=crop.radius, fitsoutname=weight.fits.output.filename,data.extn=1)
          if (verbose) { message(paste("Using", weight.fits.output.filename, "as weight map")) }
          weight.map<-weight.fits.output.filename
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }#/*fend*/ }}}
      #Error Image /*fold*/ {{{
      if ((error.map != "NONE")&(is.na(as.numeric(error.map)))) {
        if (error.map!=error.fits.output.filename) {
          if (verbose) { message(paste("Cropping Input Error Map: Outputting to", error.fits.output.filename)) }
          if (!file.exists(file.path(path.root,path.work,error.map))) { sink(type='message'); stop("Error Image does not exist at location:",file.path(path.root,path.work,error.map)) }
          if (data.error.extn>1) { sink(type='message'); stop("Internal Crop routine is unable to crop multi-extension FITS (beyond data.extn 1)") }
          crop.fits.image(ra0=ra0, dec0=dec0, path.root=file.path(path.root,path.work), inpim=error.map, crop.radius=crop.radius, fitsoutname=error.fits.output.filename,data.extn=1)
          if (verbose) { message(paste("Using", error.fits.output.filename, "as error map")) }
          error.map<-error.fits.output.filename
        } else {
          message("Crop input and output are the same file. Crop will fail, so it is skipped")
          warning("Crop input and output are the same file. Crop will fail, so it is skipped")
        }
      }#/*fend*/ }}}
    }#/*fend*/ }}}

    #Create Seperate Images Environment /*fold*/ {{{
    #Details /*fold*/ {{{
    # Create the image environment, which will contain all
    # Image Data for the LAMBDAR Proceedure - This segregates the
    # large image arrays from the regular parameter space, and
    # stops unnecessary memory usage in foreach loops
    #/*fend*/ }}}
    image.env<-new.env(parent=environment())
    #/*fend*/ }}}

    #Read in Data, Mask map, & Error map /*fold*/ {{{
    if (stop.on.missing) {
      read.images(env=NULL,quiet,showtime,outenv=image.env)
    } else {
      status<-try(read.images(env=NULL,quiet,showtime,outenv=image.env))
      if (class(status)=='try-error') {
        if (!quiet) { cat("Failed! Files missing\n") }
        sink(type="message")
        close(sink.file)
        #If we were supposed to read & reuse this table, we can't (this table never gets read)
        if (!reuse.table[f] && (f!=loop.total) && reuse.table[f+1]) { reuse.table[f+1]<-FALSE }
        next
      }
    }
    #Move Astrometry list from image env to main env /*fold*/ {{{
    astr.struc<-image.env$astr.struc
    saturation<-image.env$saturation
    rm(astr.struc, envir=image.env)
    #/*fend*/ }}}
    #/*fend*/ }}}

    #Read source catalogue /*fold*/ {{{
    if (reuse.table[f] && f!=loop.start) {
      if (!quiet) { cat(paste("   Restoring Previous Catalogue",catalogue,"   ")) }
      #Restore the previous catalogue columns
      cat.id<-saved.table[,cata.lab]
      cat.ra<-saved.table[,ra.lab]
      cat.dec<-saved.table[,dec.lab]
      cat.a<-saved.table[,semimaj.lab]
      cat.b<-saved.table[,semimin.lab]
      cat.theta<-saved.table[,theta.lab]
      flux.weight<-saved.table[,flux.weight.lab]
      contams<-saved.table[,contam.lab]
      if (!quiet) { cat("-- Done\n") }
    } else {
      #Read the catalogue
      open.catalogue(outenv=environment(),save.table=(f!=loop.total && reuse.table[f+1]))
    }
    #/*fend*/ }}}

    #If needed, read ZP Magnitude from Image header /*fold*/ {{{
    if ((magnitudes) & (mag.zp==-999)){
      mag.zp<-as.numeric(read.fitskey(mag.zp.label,paste(path.root,path.work,data.map,sep=""),hdu=data.extn))
      #If Failed, do not output magnitudes /*fold*/ {{{
      if (!is.finite(mag.zp)) {
        message("Zero Point Magnitude determination failed - Not outputting magnitudes")
        warning("Zero Point Magnitude determination failed")
        magnitudes<-FALSE
      }#/*fend*/ }}}
    }#/*fend*/ }}}

    #If wanted, Set Minimum Aperture Radius /*fold*/ {{{
    if(min.ap.rad>0){
      cat.a[cat.a<min.ap.rad]=min.ap.rad
      cat.b[cat.b<min.ap.rad]=min.ap.rad
    }
    #/*fend*/ }}}

    #Determine Initial Galaxy Sample /*fold*/ {{{
    #Details /*fold*/ {{{
    #Determine which galaxies are inside the input image,
    #and which galaxies have physical aperture parameters.
    #If either of these are not true, then the aperture is
    #discarded./*fend*/ }}}
    if (!quiet) { cat("   Determining Correct Galaxy Sample ") }

    #Get object locations in pixel space /*fold*/ {{{
    gama.pos<-ad.to.xy(cat.ra,cat.dec,astr.struc)
    cat.x<-gama.pos[,1]
    cat.y<-gama.pos[,2]
    #/*fend*/ }}}

    #-----Diagnostic-----# /*fold*/ {{{
    if (diagnostic) {
      message(paste("X.pix:", length(cat.x),"\nY.pix:",length(cat.y)))
      message(paste("min/max X.pix:", min(cat.x), max(cat.x),"\nmin/max Y.pix:",min(cat.x), max(cat.y)))
    }#/*fend*/ }}}

    #Discard any apertures that lie completely outside of the image ± 1 pixel /*fold*/ {{{
    cat.len<-length(cat.x)
    inside.mask<-!((cat.x <= 0) | (cat.x >= length(image.env$im[,1])+1) | (cat.y <= 0) | (cat.y >= length(image.env$im[1,])+1))
    #Check that something is inside the image /*fold*/ {{{
    if (length(which(inside.mask==TRUE))==0) {
      warning("No catalogue entries are centred within the limits of the input image.")
      #Notify & Close Logfile /*fold*/ {{{
      looptimestring<-as.time(proc.time()[3]-loop.start.time, digits=0)
      totaltimestring<-as.time(proc.time()[3]-start.time, digits=0)
      remaintimestring<-as.time((proc.time()[3]-start.time)/(f-loop.start+1)*(loop.total-f),digits=0)
      message(paste0('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      message(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
      message(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
      message(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
      sink(type='message')
      close(sink.file)
      if (!quiet) {
        cat(paste('- Done\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      message(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
      message(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
      message(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
        if (f !=loop.total) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
      }
      next
      #/*fend*/ }}}
    }
    #/*fend*/ }}}
    #Remove object catalogue entries /*fold*/ {{{
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
    if (exists('groups')) { groups<-groups[which(inside.mask)] }
    inside.mask<-inside.mask[which(inside.mask)]
    #/*fend*/ }}}
    if (filt.contam) {
      inside.mask<-inside.mask & !contams
      #Check that there is a science target inside the image /*fold*/ {{{
      if (length(which(inside.mask==TRUE))==0) {
        warning("No science targets are centred within the limits of the input image.")
        #Notify & Close Logfile /*fold*/ {{{
        looptimestring<-as.time(proc.time()[3]-loop.start.time, digits=0)
        totaltimestring<-as.time(proc.time()[3]-start.time, digits=0)
        remaintimestring<-as.time((proc.time()[3]-start.time)/(f-loop.start+1)*(loop.total-f),digits=0)
        message(paste0('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
        message(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
        message(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
        message(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
        sink(type='message')
        close(sink.file)
        if (!quiet) {
          cat(paste('- Done\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
        message(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
        message(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
        message(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
          if (f !=loop.total) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap\nInitialising Workspace {\n')) }
        }
        next
        #/*fend*/ }}}
      }
      #/*fend*/ }}}
    }
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects inside the image (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of supplied were outside the image )")) }
    #/*fend*/ }}}
    #/*fend*/ }}}

    #Discard any apertures that have nonphysical aperture axis values /*fold*/ {{{
    cat.len<-length(cat.x)
    inside.mask<-!((cat.a < 0)|(cat.b < 0))
    #Check that something is inside the image /*fold*/ {{{
    if (length(which(inside.mask==TRUE))==0) { sink(type="message") ; stop("No Apertures remaining have physical axis values.") }
    #/*fend*/ }}}
    #Remove object catalogue entries /*fold*/ {{{
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
    if (exists('groups')) { groups<-groups[which(inside.mask)] }
    inside.mask<-inside.mask[which(inside.mask)]
    if (num.cores < 0) { 
      registerDoParallel(cores=(min(floor(length(cat.x)/5000),abs(num.cores))))
    }
    chunk.size=length(cat.id)/getDoParWorkers()
    mpi.opts<-list(chunkSize=chunk.size)
    message("Number of objects per thread:",chunk.size)
    #/*fend*/ }}}
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects with physical aperture values (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of supplied had unphysical values )")) }
    #/*fend*/ }}}
    #/*fend*/ }}}

    if (!quiet) { cat(" - Done\n") }
    #/*fend*/ }}}

    #Set object Astrometry /*fold*/ {{{
    if (!quiet) { cat("   Setting Astrometry   ") }

    #Get pixel resolution from astrometry /*fold*/ {{{
    if (all(is.finite(astr.struc$CDELT[1:2]))) {
      #Using CDELT keywords /*fold*/ {{{
      arcsec.per.pix<-max(abs(astr.struc$CDELT[1:2]))*3600.
      #/*fend*/ }}}
    } else if (all(is.finite(astr.struc$CD[1:2,1:2]))) {
      #Using CD matrix keywords /*fold*/ {{{
      arcsec.per.pix<-max(astr.struc$CD[1,1],astr.struc$CD[2,2])*3600.
      #/*fend*/ }}}
    } else {
      #Otherwise; Error - unknown Astrometry Resolution Keyword /*fold*/ {{{
      sink(type=c("output","message")) ; stop("Data image header does not contain CD or CDELT keywords")
      #/*fend*/ }}}
    }#/*fend*/ }}}

    #Check that the minor and major axes are the correct way around /*fold*/ {{{
    if (all(cat.a <= cat.b) & !all(cat.a == cat.b)) {
      warning("All major axes are smaller than all minor axes! Are the labels in the Parameter File the correct wayy around?")
      message(paste("#################\n",
                    "WARNING:  All of the major axes are smaller than all of the minor axes. They are probably the wrong way around!\n",
                    "This is technically not allowed, because it can cause issues. So they have been swapped, and the angles rotated by 90deg.\n",
                    "That will allow us to continue anyway, but you're probably going to get garbage fits if you have accidentally supplied the axes in the wrong order (check the parameter file!)...\n#################\n"))
      tmp<-cat.a
      cat.b<-cat.a
      cat.a<-tmp
      cat.theta<-cat.theta-90
    }
    #/*fend*/ }}}

    #Set apertures with NA/NULL aperture axis or minoraxis<aperturediag to point-sources /*fold*/ {{{
    #Details /*fold*/ {{{
    #If an aperture has no provided aperture parameters
    #or if the aperture width is smaller than the pixel
    #diagonal, the object should be treated as a point source.
    #This is because if an aperture is too small/thin, then
    #positional information of the aperture inside the pixel
    #is completely lost, and convolution with a PSF
    #will be inaccurate because it cannot be accurately
    #convolved. /*fend*/ }}}
    #Get pixel diagonal size in arcsec /*fold*/ {{{
    diag.arcsec<-abs(arcsec.per.pix)*sqrt(2)/2
    #/*fend*/ }}}
    #Make needed changes, and notify /*fold*/ {{{
    message("Forcing ",length(which((cat.b<diag.arcsec)|!is.finite(cat.a)))," apertures to be point sources")
    cat.a[which((cat.b<diag.arcsec)|!is.finite(cat.a))]<-0
    cat.theta[which((cat.b<diag.arcsec)|!is.finite(cat.theta))]<-0
    cat.b[which((cat.b<diag.arcsec)|!is.finite(cat.b))]<-0
    #/*fend*/ }}}
    #/*fend*/ }}}

    if (!quiet) { cat(" - Done\n") }
    #/*fend*/ }}}

    #-----Diagnostic-----# /*fold*/ {{{
    if (diagnostic) {
      message(paste("X.pix:", length(cat.x),"\nY.pix:",length(cat.y)))
      message(paste("min/max X.pix:", min(cat.x), max(cat.x),"\nmin/max Y.pix:",min(cat.x), max(cat.y)))
    }#/*fend*/ }}}

    #If wanted, check memory-safe /*fold*/ {{{
    if ((!mpi.backend)&&(mem.safe & is.finite(mem.lim))) {
      #Check that computation is able to be performed within memory limits /*fold*/ {{{
      if (!quiet) { cat("   Checking Memory Usage & Limits {") }
      #Details /*fold*/ {{{
      #Perform a rudimentary check that the image isn't too
      #large for the Available system memory, using
      # approx aperture memory:
      #     apmem = napertures*apsizeinBits*nLists
      # where
      #     nLists=12
      #     apsizeinbits = (90th percentile{(semimajor.axis)}*2/arcsec.per.pix)^2*bitsperpixel
      #
      # approx image memory during calculations:
      #     immem = nimages*imsizeinBytes*nThreads
      #
      # current memory usage from parameter/image storage
      #     currmem = memory available at initialisation - memory currently used.
      #/*fend*/ }}}
      #Get System Memory Usage (platform dependant) at Current Time /*fold*/ {{{
      if (sys.name=="Linux") {
        #Linux Systems /*fold*/ {{{
        #Memory returned in kBytes. Convert to bits.
        memTot<-as.numeric(system("awk '/MemTotal:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
        memAct<-as.numeric(system("awk '/Active:/ {print $2}' /proc/meminfo", intern=TRUE))*1E3*8
        memCur<-memTot-memAct
        if ((length(memCur)==0)||!is.finite(memCur)) {
          warning("Memory determination failed. Setting to Inf.")
          memCur<-Inf
        }
        #/*fend*/ }}}
      } else if (sys.name=="Darwin") {
        #Mac Systems /*fold*/ {{{
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
        #/*fend*/ }}}
      } else if (sys.name=="Windows") {
        #Windows Machines /*fold*/ {{{
        #Memory returned in Bytes. Convert to Bits
        memCur<-as.numeric(memory.limit())*8
        if (!is.finite(memCur)) {
          warning("Memory determination failed. Setting to Inf.")
          memCur<-Inf
        }
        #/*fend*/ }}}
      } else {
        #Any Others, warn & set to Inf /*fold*/ {{{
        warning("Unknown Operating System Name. Cannot determine Memory Usages.\nUsing Infinity")
        message("Unknown Operating System Name. Cannot determine Memory Usages.\nUsing Infinity")
        memCur<-Inf
        #/*fend*/ }}}
      }
      memCur<-mem.lim-memCur
      #/*fend*/ }}}
      #Aperture Memory requirements /*fold*/ {{{
      cat.len<-length(cat.id)
      #Use 90th Quantile of aperture semimajor axes /*fold*/ {{{
      aprad.quant<-quantile(cat.a[which(cat.a>0)],0.9)
      if (is.na(aprad.quant)) {
        #All apertures are point sources - make a default width @ 10 pix
        aprad.quant<-10
      }
      #/*fend*/ }}}
      apsizeinbits<-(2*aprad.quant/arcsec.per.pix)^2*64*2
      apmem<-cat.len*apsizeinbits*12
      #/*fend*/ }}}
      #Image Memory requirements /*fold*/ {{{
      nimage=4
      if (length(image.env$imm)>1){ nimage=nimage+1 }
      if (length(image.env$ime)>1){ nimage=nimage+1 }
      if (sourcemask){ nimage=nimage+1 }
      immem<-nimage*lsos(pattern="im",envir=image.env)[1,'Size']
      #/*fend*/ }}}
      #Check Memory Allocation is less than available free memory /*fold*/ {{{
      if ((apmem+immem+memCur) >= mem.lim) {
          #If this is too much, less threads cannot help. Stop /*fold*/ {{{
          cat("Failed\n")
          sink(type='message')
          if (force.safe) {
          sink(type='message')
          stop(paste("This computation may exceed the available memory on this machine.",
                     "The required memory (",(apmem+immem+memCur)*1E-9/8,"Gb) is greater than that which is available to the system (",(mem.lim)*1E-9/8,"Gb).",
                     "\nHowever, using the in-built crop function, seperating the image into smaller chuncks will enable computation.",
                     "\nThe memory usage is",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                     "\n",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n"))
          } else {
          message(paste("This computation may exceed the available memory on this machine.",
                     "The required memory (",(apmem+immem+memCur)*1E-9/8,"Gb) is greater than that which is available to the system (",(mem.lim)*1E-9/8,"Gb).",
                     "\nHowever, using the in-built crop function, seperating the image into smaller chuncks will enable computation.",
                     "\nThe memory usage is",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                     "\n",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n"))
          }
          #/*fend*/ }}}
      }#/*fend*/ }}}
      #Notify Memory Usage Status /*fold*/ {{{
      if (!quiet) {
        cat(paste("\n       This computation will require approximately ",round((apmem+immem)*1E-9/8,digits=3)," Gb of memory.\n",
                     "      The memory usage is ",round(apmem*1E-9/8,digits=3),"Gb for apertures/results/storage/etc, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                  "\n       ",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n",
                     "       Total memory available for the system is ",(mem.lim)*1E-9/8," Gb.\n",sep=""))
      }
      message(paste("This computation will require approximately",round((apmem+immem)*1E-9/8,digits=3),"Gb of memory.\n",
                     "      The memory usage is ",round(apmem*1E-9/8,digits=3),"Gb for apertures, ",round(immem*1E-9/8,digits=3),"Gb for images.",
                  "\n       ",round(memCur*1E-9/8,digits=3),"Gb was assigned during Parameter & Image initialisation.\n",
                    "Total memory available for the system is",(mem.lim)*1E-9/8,"Gb."))
      #/*fend*/ }}}
      #/*fend*/ }}}
    }
    #/*fend*/ }}}

    #Get beam area in pixels /*fold*/ {{{
    beam.area.pix<-beam.area.input.as/(arcsec.per.pix)^2.
    #/*fend*/ }}}

    #Finished setting astrometry. Notify /*fold*/ {{{
    if (!quiet) {
      cat("   } - Done\n")
      cat('} Initialisation Complete ')
    }
    if (showtime) { cat(paste(' (  Loop Time Elapsed (s): ',round(proc.time()[3]-loop.start.time, digits=3),'  )\n')) }
    #/*fend*/ }}}

    #-----Diagnostic-----# /*fold*/ {{{
    if (diagnostic) {
    message(paste('Arcsec per pixel in map: ', arcsec.per.pix))
    message(paste('Beam area (from observers manual) converted into pixel units: ', beam.area.pix))
    }#/*fend*/ }}}

    #Send Parameter to logfile /*fold*/ {{{
    sink(sink.file, type="output")
    cat("Parameters used in this run:\n")
    print(lsos(envir=environment(), head=FALSE))
    cat("Images used in this run:\n")
    print(lsos(envir=image.env, head=FALSE))
    sink(type="output")
    #/*fend*/ }}}

    #Run Flux Measurements /*fold*/ {{{
    loopresults<-flux.measurements()
    if (do.return) {
      results<-c(results, loopresults)
    } else {
      loopresults<-NULL
    }
    #/*fend*/ }}}


    #Notify & Close Logfile /*fold*/ {{{
    looptimestring<-as.time(proc.time()[3]-loop.start.time, digits=0)
    totaltimestring<-as.time(proc.time()[3]-start.time, digits=0)
    remaintimestring<-as.time((proc.time()[3]-start.time)/(f-loop.start+1)*(loop.total-f),digits=0)
    message(paste0('\n-----------------------------------------------------\nDatamap Complete\n'))
    message(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
    message(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
    message(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
    sink(type='message')
    close(sink.file)
    #/*fend*/ }}}

    #Loop Completed - Print Update, Loop /*fold*/ {{{
    if (!quiet) {
      cat(paste0('-----------------------------------------------------\nDatamap Complete\n'))
      cat(paste0('Loop Completed: Indiv. Loop Time Elapsed: ',looptimestring,'\n'))
      cat(paste0('                      Total Time Elapsed: ',totaltimestring,'\n'))
      cat(paste0('                  Approx. Time Remaining: ',remaintimestring,'\n'))
      if (f !=loop.total) { cat(paste('-----------------------------------------------------\nLooping to Next DataMap (',f+1,' of ',loop.total,')\nInitialising Workspace {\n',sep="")) }
    }
    #/*fend*/ }}}
  }#/*fend*/ }}}

  #Remove State Preservation Files /*fold*/ {{{
  if (loop.total>1) {
    system("rm -f .LambdarActiveLoop.txt .LambdarParameters.Rdata")
  }#/*fend*/ }}}

  #Program Completed - Print closing /*fold*/ {{{
  if (!quiet) {
    cat(paste('-----------------------------------------------------\nProgram Complete\n'))
  }
  if (do.return) {
    return=results
  } else {
    return="There are no results being returned because do.return=FALSE"
  }
  #/*fend*/ }}}

}
#As Time functions /*fold*/ {{{
as.time<-function(sec,digits) {
  if (sec > 60*60*24) {
    day<-floor(sec/(60*60*24))
    hr<-floor((sec-day*(60*60*24))/(60*60))
    min<-round((sec-day*(60*60*24)-hr*(60*60))/60,digits=digits)
    timestr<-paste(day,'day',hr,'hr',min,'min')
  } else if (sec > 60*60) {
    hr<-floor(sec/(60*60))
    min<-floor((sec-hr*(60*60))/60)
    sec<-round(sec-hr*(60*60)-min*60,digits=digits)
    timestr<-paste(hr,'hr',min,'min',sec,'sec')
  } else if (sec > 60 ) {
    min<-floor(sec/60)
    sec<-round(sec-min*60,digits=digits)
    timestr<-paste(min,'min',sec,'sec')
  } else {
    sec<-round(sec,digits=digits)
    timestr<-paste(sec,'sec')
  }
  return=timestr
}
#/*fend*/ }}}
