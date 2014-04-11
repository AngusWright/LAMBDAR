readimage <-
function(env=NULL, quiet=FALSE, showtime=FALSE, outenv=NULL){

  # Load Parameter Space {{{
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #}}}

  #Read Data, Mask Map, and Error Map {{{

  if (!quiet) { cat(paste("   Reading Data from Image",datamap,"   ")) }
  if (showtime) { timer<-proc.time() }

  #Test Read of Data Image for errors {{{
  im_fits<-try(read.fits(paste(pathroot,datamap,sep=""),hdu=extn, comments=FALSE))
  if (class(im_fits)=="try-error") {
    #Stop on Error
    geterrmessage()
    stop("Data Image File read failed")
  }
  #}}}

  #Read Data Image {{{
  hdr=im_fits$hdr[[1]][which(im_fits$hdr[[1]][,"key"]!="COMMENT"),]
  hdr_str<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
  im<-im_fits$dat[[1]]
  #}}}

  #Remove NA/NaN/Inf's {{{
  im[which(!is.finite(im))]<-0.0
  #}}}

  #Setup Astrometry Structure {{{
  astr_struc<-read.astr(paste(pathroot,datamap,sep=""))
  #}}}

  #Finished Setting Astrometry {{{
  if (showtime) { cat(paste(" - Done (Image Read took",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }
  #}}}

  #Read Weightmap if needed {{{
  if ((wgtmap!="NONE")&(maskmap=="NONE"|errormap=="NONE")) {
    if (!quiet) { cat(paste("   Reading Data from Weight Map",wgtmap,"   ")) }
    #Try read weight map {{{
    imwt_fits<-try(read.fits(paste(pathroot,wgtmap,sep=""),hdu=extnwgt,comments=FALSE))
    if (class(imwt_fits)=="try-error") {
      #Stop on Error
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
  if (maskmap=='NONE') {
    #No Mask Present {{{
    if (wgtmap=='NONE') {
      #No Weight map for mask generation; Generate transparent mask {{{
      if (!quiet) { cat(paste("   Generating Mask Map ")) }
      #If no mask, set mask to 1: transparent will be used everywhere
      imm<-1
      hdr_mask<-NULL
      #}}}
    } else {
      #Generate Mask {{{
      if (!quiet) { cat(paste("   Generating Mask Map from Weight Map ")) }
      #Make mask same dimensions as weightmap {{{
      imm<-array(1,dim=dim(imwt))
      #}}}
      #Make mask 0 where weightmap == weightmap zero point {{{
      imm[which(imwt==wgtzp)]<-0
      hdr_mask<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
      #}}}
      #}}}
    }#}}}
  } else {
    #Mask Present, Read {{{
    if (!quiet) { cat(paste("   Reading Data from MaskMap",maskmap,"   ")) }
    #Test Read of Mask Map for errors {{{
    imm_fits<-try(read.fits(paste(pathroot,maskmap,sep=""),hdu=extnmask, comments=FALSE))
    if (class(imm_fits)=="try-error") {
      #Stop on Error
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
    hdr_mask<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
    #}}}
    #}}}
  }
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  } else if (!quiet) { cat(" - Done\n") }
  #}}}

  #Read error map {{{
  if (errormap=='NONE') {
    #If no error map, generate the sigma-map using the provided gain and the image_map {{{
    if (!quiet) { cat(paste("   Generating Error Map ")) }
    if (wgtmap=="NONE") {
      #No Weightmap; sigma map is poisson-like {{{
      ime<-sqrt(abs(im))
      #}}}
    } else {
      #Weightmap present; sigma map is image/image*gain {{{
      ime<-abs(im)/sqrt(abs(im*(1/imwt)))
      #}}}
    }
    hdr_err<-hdr_str
    #}}}
  } else {
    #Check if Gain value or Error map File {{{
    if (!is.na(as.numeric(errormap))) {
      #Gain Value {{{
      if (!quiet) { cat(paste("   Using Supplied Single Gain value of",errormap," for errors   ")) }
      ime<-abs(im)/sqrt(abs(im*as.numeric(errormap)))
      #}}}
    } else {
      #Sigma Map File {{{
      if (!quiet) { cat(paste("   Reading Data from ErrorMap",errormap,"   ")) }
      #Test Read of Error Map for errors {{{
      ime_fits<-try(read.fits(paste(pathroot,errormap,sep=""),hdu=extnerr, comments=FALSE))
      if (class(ime_fits)=="try-error") {
        #Stop on Error
        geterrmessage()
        stop("Error Map File read failed: Provided Entry is neither a file, nor NONE, nor a numeric Gain")
      }
      #}}}
      #Read Sigma Map {{{
      hdr=ime_fits$hdr[[1]][which(ime_fits$hdr[[1]][,"key"]!="COMMENT"),]
      hdr_err<-as.data.frame(hdr[,"value"], row.names=hdr[,"key"], stringsAsFactors=FALSE)
      ime<-ime_fits$dat[[1]]
      #}}}
      #Remove NA/NaN/Inf's {{{
      ime[which(!is.finite(ime))]<-0.0
      #}}}
      #}}}
    }
    #}}}
  }
  #Scale error map by Efactor {{{
  ime=ime*Efactor
  #}}}
  if (showtime) { cat(paste(" - Done (",round(proc.time()[3]-timer[3], digits=3),"sec )\n"))
  timer<-proc.time()
  } else if (!quiet) { cat(" - Done\n") }
  #}}}

  #Notify {{{
  if (verbose) {
    message(paste('ReadSDP: (Data,Error,Mask)'))
    message(paste(datamap))
    message(paste(errormap))
    message(paste(maskmap))
  }#}}}

  #Parse Parameter Space {{{
  detach(env)
  assign("imm"       , imm       , envir = outenv)
  assign("ime"       , ime       , envir = outenv)
  assign("im"        , im        , envir = outenv)
  assign("hdr_str"   , hdr_str   , envir = outenv)
  assign("hdr_err"   , hdr_err   , envir = outenv)
  assign("hdr_mask"  , hdr_mask  , envir = outenv)
  assign("astr_struc", astr_struc, envir = outenv)
  return=NULL
  #}}}

}
