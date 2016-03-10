read.images <-
function(outenv=parent.env(environment()), quiet=FALSE, showtime=FALSE, env=NULL){

  # Load Parameter Space {{{
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  if(is.null(outenv)&!is.null(env)) { outenv<-env }
  else if (is.null(outenv)) {
    warning("Output Environment cannot be NULL; using parent env")
    outenv<-parent.env(environment())
  }
  #}}}

  #Read Data, Mask Map, and Error Map {{{

  if (!quiet) { cat(paste("   Reading Data from Image",data.map,"   ")) }
  if (showtime) { timer<-proc.time() }

  #Check data file exists {{{
  if (!file.exists(paste(path.root,path.work,data.map,sep=""))) {
    sink(type='message')
    stop("Data Image does not exist at location specified:",paste(path.root,path.work,data.map,sep=""))
  }
  #}}}

  #Setup Astrometry Structure {{{
  astr.struc<-read.astrometry(paste(path.root,path.work,data.map,sep=""),hdu=data.extn)
  # Check WCS {{{
  if (any(is.na(astr.struc$CTYPE[c(1:2)]))) {
  } else if (all(grepl("TAN", astr.struc$CTYPE[c(1:2)]))) {
    #WCS is TAN Gnomonic; Continue without issue
  } else if (all(grepl("SIN", astr.struc$CTYPE[c(1:2)]))) {
    #WCS is SIN Orthographic; if no rotation, continue without error
    if (!(astr.struc$CD[1,2]==0 && astr.struc$CD[2,1]==0)) {
      #Error; WCS has rotation
      sink(type='message')
      stop("Fits file has WCS 'SIN' with Rotation; currently only TAN Gnomonic & SIN without rotation are compatible.")
    }
  } else {
    #Error; WCS not compatible
    sink(type='message')
    stop("Fits file has an incompatible WCS; currently only TAN Gnomonic & SIN without rotation are compatible.")
  }
  #}}}
  #}}}

  #Test Read of Data Image for errors {{{
  im_fits<-try(read.fits(paste(path.root,path.work,data.map,sep=""),hdu=data.extn, comments=FALSE),silent=TRUE)
  if (class(im_fits)=="try-error") {
    #Stop on Error
    sink(type='message')
    geterrmessage()
    stop("Data Image File read failed")
  }
  #}}}

  #Read Data Image {{{
  hdr=im_fits$hdr[[1]][which(im_fits$hdr[[1]][,"key"]!="COMMENT"),]
  data.hdr<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
  im<-im_fits$dat[[1]]
  #}}}

  #Read Saturation Level {{{
  if (!is.finite(saturation)) {
    saturation<-try(as.numeric(read.fitskey(satur.label,paste(path.root,path.work,data.map,sep=""),hdu=data.extn)),silent=TRUE)
    if (class(saturation)=='try-error' | is.na(saturation)) {
      saturation<-Inf
    }
  }
  message(paste0("Using Saturation level: ",saturation))
  #}}}

  #Finished Setting Astrometry {{{
  if (showtime) { cat(paste(" - Done (Image Read took",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }
  #}}}

  #Read Weightmap if needed {{{
  if ((weight.map!="NONE")&(mask.map=="NONE"|error.map=="NONE")) {
    if (!quiet) { cat(paste("   Reading Data from Weight Map",weight.map,"   ")) }
    #Try read weight map {{{
    imwt_fits<-try(read.fits(paste(path.root,path.work,weight.map,sep=""),hdu=data.weight.extn,comments=FALSE),silent=TRUE)
    if (class(imwt_fits)=="try-error") {
      #Stop on Error
      sink(type='message')
      geterrmessage()
      stop("Weight Map read failed")
    }
    #}}}
    #Read Weightmap {{{
    hdr=imwt_fits$hdr[[1]][which(imwt_fits$hdr[[1]][,"key"]!="COMMENT"),]
    imwt<-imwt_fits$dat[[1]]
    if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
    timer<-proc.time()
    } else if (!quiet) { cat(" - Done\n") }
    #}}}
  }#}}}

  #Read mask map {{{
  if (mask.map=='NONE') {
    #No Mask Present {{{
    if (weight.map=='NONE') {
      #No Weight map for mask generation; Generate transparent mask {{{
      if (!quiet) { cat(paste("   Generating Mask Map ")) }
      #If no mask, set mask to 1: transparent will be used everywhere
      imm<-1
      mask.hdr<-NULL
      #}}}
    } else {
      #Generate Mask {{{
      if (!quiet) { cat(paste("   Generating Mask Map from Weight Map ")) }
      #Make mask same dimensions as weightmap {{{
      imm<-array(1,dim=dim(imwt))
      #}}}
      #Make mask 0 where weightmap == weightmap zero point {{{
      imm[which(imwt==wgt.zp)]<-0
      mask.hdr<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
      #}}}
      #}}}
    }#}}}
  } else {
    #Mask Present, Read {{{
    if (!quiet) { cat(paste("   Reading Data from MaskMap",mask.map,"   ")) }
    #Test Read of Mask Map for errors {{{
    imm_fits<-try(read.fits(paste(path.root,path.work,mask.map,sep=""),hdu=data.mask.extn, comments=FALSE),silent=TRUE)
    if (class(imm_fits)=="try-error") {
      #Stop on Error
      sink(type='message')
      geterrmessage()
      stop("Mask Map read failed")
    }
    #}}}
    #Read Mask Map {{{
    hdr=imm_fits$hdr[[1]][which(imm_fits$hdr[[1]][,"key"]!="COMMENT"),]
    imm<-imm_fits$dat[[1]]
    #}}}
    #Remove NA/NaN/Inf's {{{
    imm[which(!is.finite(imm))]<-0.0
    mask.hdr<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
    #}}}
    #}}}
  }
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }
  #}}}

  #Read error map {{{
  if (error.map=='NONE') {
    #If no error map, generate the sigma-map using the provided gain and the image_map {{{
    if (!quiet) { cat(paste("   Generating Error Map ")) }
    if (weight.map=="NONE") {
      #Try reading gain from header {{{
      gain<-try(as.numeric(read.fitskey(gain.label,file=paste(path.root,path.work,data.map,sep=""),hdu=data.extn)),silent=TRUE)
      if ((class(gain)=="try-error")|is.na(gain)){
        #No Gain; No Weightmap; SNR map is poisson-like {{{
        message(paste0("No Gain supplied or able to be read from header; making Poisson Error Map"))
        ime<-sqrt(abs(im))
        #}}}
      } else {
        #Max Gain; No Weightmap; Sigma map is sqrt(abs(image/gain)) {{{
        message(paste0("Using Gain level from header: ",gain))
        ime<-sqrt(abs(im/(gain)))
        #}}}
      }
      #}}}
    } else {
      #Use Weight Map in generating Error Map {{{
      message(paste0("Generating Errormap from Weightmap"))
      gain<-try(as.numeric(read.fitskey(gain.label,file=paste(path.root,path.work,data.map,sep=""),hdu=data.extn)),silent=TRUE)
      if ((class(gain)=="try-error")|is.na(gain)){
        #No Max gain; Weightmap present; sigma map is 1/sqrt(weight.map) {{{
        message(paste0("No Gain supplied or able to be read from header; Using Weight-map as absolute 1/Var(x) to generate the Error Map."))
        #sigma map = 1/sqrt(weight.map)
        ime<-1/sqrt(imwt)
        #Determine Implied Gain: gain=abs(im)/(sigma)^2
        gain<-abs(im)/ime^2
        ind<-which(imwt!=wgt.zp)
        message(paste0("Implied Median (Max) Equivalent Gain level from Weight Map: ",median(gain[ind],na.rm=TRUE), "(",max(gain[ind],na.rm=TRUE),")"))
        gain<-max(gain[ind],na.rm=TRUE)
        #}}}
      } else {
        #Max Gain; Weightmap present; sigma map is 1/sqrt(image/gain) {{{
        message(paste0("Using Max Equiv. Gain from header to scale Weight Map when generating the Error Map."))
        message(paste0("Max Equivalent Gain level from Header: ",gain))
        #convert from wt ~ 1/var(x) to ~ sig(x)
        imwt<-1/sqrt(imwt)
        ind<-which(is.finite(imwt))
        #use relative sigma to get absolute pixel varying gain
        imwt<-imwt/max(imwt[ind],na.rm=TRUE)*gain
        #sigma map = sqrt(abs(im)/gain)
        ime<-sqrt(abs(im/(imwt)))
        #}}}
      }
      #}}}
    }
    error.hdr<-data.hdr
    #}}}
  } else {
    #Check if Gain value or Error map File {{{
    if (!is.na(suppressWarnings(as.numeric(error.map)))) {
      #Gain Value {{{
      gain<-error.map
      if (!quiet) { cat(paste("   Using Supplied Single Gain value of",error.map," for errors   ")) }
      #Sigma map is sqrt(abs(im)/gain)
      ime<-sqrt(abs(im/as.numeric(error.map)))
      error.hdr<-data.hdr
      #Remove NA/NaN/Inf's {{{
      ime[which(!is.finite(ime))]<-0.0
      #}}}
      #}}}
    } else {
      #Sigma Map File {{{
      if (!quiet) { cat(paste("   Reading Data from ErrorMap",error.map,"   ")) }
      gain<-NA
      #Test Read of Error Map for errors {{{
      ime_fits<-try(read.fits(paste(path.root,path.work,error.map,sep=""),hdu=data.error.extn, comments=FALSE),silent=TRUE)
      if (class(ime_fits)=="try-error") {
        #Stop on Error
        sink(type='message')
        geterrmessage()
        stop("Error Map File read failed: Provided Entry is neither a file, nor NONE, nor a numeric Gain")
      }
      #}}}
      #Read Sigma Map {{{
      hdr=ime_fits$hdr[[1]][which(ime_fits$hdr[[1]][,"key"]!="COMMENT"),]
      error.hdr<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
      ime<-ime_fits$dat[[1]]
      #}}}
      #Remove NA/NaN/Inf's {{{
      ime[which(!is.finite(ime))]<-0.0
      #}}}
      #If no mask, generate one using Errormap {{{
      if (length(imm)==1) {
        if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
        timer<-proc.time()
        } else if (!quiet) { cat(" - Done\n") }
        if (!quiet) { cat(paste("   Generating Mask Map from Weight Map ")) }
        # Use Errormap Footprint
        imm<-ime*0
        # Wherever there is non-zero error, make transparent.
        imm[which(ime>0)]<-1
      }
      #}}}
      #}}}
    }
    #}}}
  }
  #Scale error map by error.factor {{{
  ime=ime*error.factor
  #}}}
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }
  #}}}
  #}}}

  #Remove NA/NaN/Inf's {{{
  im[which(!is.finite(im))]<-0.0
  #}}}


  #Notify {{{
  if (verbose) {
    message(paste('ReadSDP: (Data,Error,Mask)'))
    message(paste(data.map))
    message(paste(error.map))
    message(paste(mask.map))
  }#}}}

  #Parse Parameter Space {{{
  if (!is.null(env)) { detach(env) }
  assign("imm"       , imm       , envir = outenv)
  assign("ime"       , ime       , envir = outenv)
  assign("im"        , im        , envir = outenv)
  assign("data.hdr"   , data.hdr   , envir = outenv)
  assign("error.hdr"   , error.hdr   , envir = outenv)
  assign("mask.hdr"  , mask.hdr  , envir = outenv)
  assign("astr.struc", astr.struc, envir = outenv)
  assign("saturation", saturation, envir = outenv)
  assign("gain"      , gain      , envir = outenv)
  return=NULL
  #}}}

}
