readimage <-
function(env=NULL, quiet=FALSE, showtime=FALSE, outenv=NULL){

  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))

  #Read Data, Mask Map, and Error Map

  if (!quiet) { cat(paste("   Reading Data from Image",datamap,"   ")) }
  if (showtime) { timer<-proc.time() }

  #Test Read of Data Image for errors
  error<-try(read.fits(paste(pathroot,datamap,sep=""),hdu=extn, comments=FALSE))
  if (class(error)=="try error") {
    #Stop on Error
    geterrmessage()
    stop("Data Image File read failed")
  }
  #Read Data Image
  im_fits<-read.fits(paste(pathroot,datamap,sep=""),hdu=extn, comments=FALSE)
  hdr=im_fits$hdr[[1]][which(im_fits$hdr[[1]][,"key"]!="COMMENT"),]
  hdr_str<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
  im<-im_fits$dat[[1]]
  #Remove NA/NaN/Inf's
  im[which(!is.finite(im))]<-0.0

  #Setup Astrometry Structure
  astr_struc<-read.astr(paste(pathroot,datamap,sep=""))

  #Finished Setting Astrometry
  if (showtime) { cat(paste(" - Done (Image Read took",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }

  #Read error map
  if (errormap=='NONE') {
    if (!quiet) { cat(paste("   Generating Error Map ")) }
    #If no error map, set map to 1: equal weighting will be used everywhere
    ime<-1
    hdr_err<-NULL
  } else {
    #If map present, read
    if (!quiet) { cat(paste("   Reading Data from ErrorMap",errormap,"   ")) }
    #Test Read of Error Map for errors
    error<-try(read.fits(paste(pathroot,errormap,sep=""),hdu=extnerr, comments=FALSE))
    if (class(error)=="try error") {
      #Stop on Error
      geterrmessage()
      stop("Error Map File read failed")
    }
    ime_fits=read.fits(paste(pathroot,errormap,sep=""),hdu=extnerr, comments=FALSE)
    hdr=ime_fits$hdr[[1]][which(ime_fits$hdr[[1]][,"key"]!="COMMENT"),]
    ime<-ime_fits$dat[[1]]
    #Remove NA/NaN/Inf's
    ime[which(!is.finite(ime))]<-0.0
    hdr_err<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
  }
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }

  #Read mask map
  if (maskmap=='NONE') {
    if (!quiet) { cat(paste("   Generating Mask Map ")) }
    #If no mask, set mask to 1: transparent will be used everywhere
    imm<-1
    hdr_mask<-NULL
  } else {
    #If mask present, read
    if (!quiet) { cat(paste("   Reading Data from MaskMap",maskmap,"   ")) }
    #Test Read of Mask Map for errors
    error<-try(read.fits(paste(pathroot,maskmap,sep=""),hdu=extnerr, comments=FALSE))
    if (class(error)=="try error") {
      #Stop on Error
      geterrmessage()
      stop("Mask Map read failed")
    }
    imm_fits=read.fits(paste(pathroot,maskmap,sep=""),hdu=extnmask, comments=FALSE)
    hdr=imm_fits$hdr[[1]][which(imm_fits$hdr[[1]][,"key"]!="COMMENT"),]
    imm<-imm_fits$dat[[1]]
    #Remove NA/NaN/Inf's
    imm[which(!is.finite(imm))]<-0.0
    hdr_mask<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
  }
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }


  #Scale error map by Efactor
  ime=ime*Efactor

  #-----Verbose Output-----#
  if (verbose) {
    message(paste('ReadSDP: (Data,Error,Mask)'))
    message(paste(datamap))
    message(paste(errormap))
    message(paste(maskmap))
  }

  #Parse Parameter Space
  detach(env)
  assign("imm"       , imm       , envir = outenv)
  assign("ime"       , ime       , envir = outenv)
  assign("im"        , im        , envir = outenv)
  assign("im_fits"   , im_fits   , envir = outenv)
  assign("hdr_str"   , hdr_str   , envir = outenv)
  assign("hdr_err"   , hdr_err   , envir = outenv)
  assign("hdr_mask"  , hdr_mask  , envir = outenv)
  assign("astr_struc", astr_struc, envir = outenv)

}
