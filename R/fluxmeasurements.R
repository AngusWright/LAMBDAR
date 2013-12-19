fluxmeasurements <-
function(env=NULL) {
#function(id_g,ra_g,dec_g,x_g,y_g,theta_g,a_g,b_g,im,hdr,ime,imm,pathroot,
#pathout,psfmap,lambda,field,conf,fluxcorr,fluxweight,beamarea_pix,asperpix,Jybm,
#sourcemask,nopsf,filtcontam, contams) {
  # Procedure measures the fluxes present in a supplied image
  # inside catalogued apertures.

  # Load Parameter Space
  if (is.null(env)) {
    stop("No Parameter Space Environment Specified in Call")
  }
  attach(env)

  if (!quiet) { cat('Begining Flux Measurements\n') }
  message('} Initialisation Complete\nBegining Flux Measurements')

  timer<-proc.time()
  if (!quiet) { cat('Getting PSF details') }
  message('Getting PSF details')
  # If we are filtering apertures with the PSF
  if (!(nopsf)) {
    #Calculate PSF - if one has not be supplied,
    #then a gaussian PSF will be created. Stampsize of PSF
    #should be = maximum stampsize: max aperture * stamp mult
    psf<-readpsf(environment(),paste(pathroot,psfmap,sep=""),asperpix,defbuff*max(a_g),stampmult,gauss_fwhm_as=gauss_fwhm_as,normalize=TRUE)
    if (verbose) { message(paste("Maxima of the PSF is at pixel", which(psf == max(psf)),"and has value",max(psf))) }
    #Normalise Beam Area
    beamarea_nn<-sumpsf
    beamarea_n<-as.single(sum(psf))
    #-----Diagnostic-----#
    if (diagnostic) { message(paste('Beam area before/after norm: ',beamarea_nn,beamarea_n)) }
  } else {
    #As we are not using PSF filtering, we set stamp size ourselves
    #and set the beamarea to unity
    beamarea_n<-1.0
    psffwhm<-0
  }
  if (beamarea_pix == 0) {
    # If possible, use the beamarea just determined at high resolution
    beamarea<-beamarea_n
  } else {
    # Otherwise just use the input beamarea
    beamarea<-beamarea_pix
  }
  #-----Diagnostic-----#
  if (diagnostic) { message(paste('Beam area adopted: ',beamarea)) }

  #Convert images from Jy/beam to Jy/pixel
  if (Jybm) {
    image.env$im<-image.env$im/beamarea
    if (length(image.env$ime)!=1) { image.env$ime<-image.env$ime/beamarea }
    conf<-conf/beamarea
  }
  if (showtime) { cat(paste(' - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )\n'))
    message(paste('Getting PSF details - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat(' - Done\n') }

  #---------------------------------------------------------------
  #PART ONE:  SINGLE APERTURE MASKS - perform on contams & sources
  #---------------------------------------------------------------

  #Convert precise (i.e. non-integer) pixel values into nearest
  #whole pixel values
  pix<-objpos2objpix(x_g, y_g)
  x_p<-pix["x_p",]
  y_p<-pix["y_p",]

  #Create an array of stamps containing the apertures for all objects
  timer=system.time(sa<-make_sa_mask(environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make SA Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }

  #Discard any apertures that were zero'd in the make process
  totsa<-foreach(i=1:length(sa), .combine='c') %dopar% { if (max(sa[[i]])>1){warning(paste("Max of Aperture",i,"is",max(sa[[i]]))) } ; sum(sa[[i]]) }
  insidemask<-(totsa > 0)
  if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Single Apertures are inside the Mask.") }  # Nothing inside the mask
  sa<-sa[which(insidemask)]
  stamplen<-stamplen[which(insidemask)]
  #This is to stop crashes when only 1 aperture is
  #present inside the image.
  if (length(which(insidemask))==1) {
    stamp_lims<-rbind(stamp_lims[which(insidemask),],stamp_lims[which(insidemask),])
  } else {
    stamp_lims<-stamp_lims[which(insidemask),]
  }
  x_p<-x_p[which(insidemask)]
  y_p<-y_p[which(insidemask)]
  x_g<-x_g[which(insidemask)]
  y_g<-y_g[which(insidemask)]
  id_g<-id_g[which(insidemask)]
  ra_g<-ra_g[which(insidemask)]
  dec_g<-dec_g[which(insidemask)]
  theta_g<-theta_g[which(insidemask)]
  a_g<-a_g[which(insidemask)]
  a_g_pix<-a_g_pix[which(insidemask)]
  b_g<-b_g[which(insidemask)]
  if (filtcontam) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]

  #Initialise arrays
  npos<-length(id_g)

  #Create an full mask of all apertures in their correct locations
  timer=system.time(image.env$aa<-make_a_mask(environment(), sa, dim(image.env$im)))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
   message(paste('Make A Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }

  #Do we want to output the All Apertures Mask?
  if (makeaamask) {
    if (!quiet) { cat(paste('Outputting All Apertures Mask to',aafilename,"   ")) }
    #Write All Apertures Mask to file
    writefitsout(paste(pathout,aafilename,sep=""),image.env$aa,image.env$hdr_str,nochange=TRUE)
    if (!quiet) { cat(" - Done\n") }
  }

  #---------------------------------------------------------------
  #PART TWO:  SINGLE FILTERED APERTURE MASKS
  #---------------------------------------------------------------

  #Do we want to filter the apertures using the PSF?
  if ((nopsf)&(length(which(fluxweight!=1))!=0)) {
    #If not convolving with psf, and if all the fluxweights are not unity,
    #then skip filtering, duplicate the arrays, and only weight the stamps/apertures
    if (verbose) { message("NoPSF: Convolved Apertures are identical to Simple Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    #Perform aperture weighting
    timer=system.time(wsfa<-make_sfa_mask(environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #Create a full mask of stamps/apertures
    timer=system.time(image.env$wfa<-make_a_mask(environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  } else if ((nopsf)&(length(which(fluxweight!=1))==0)) {
    #If not convolving with psf, and if all the fluxweights are unity,
    #then skip filtering and weighting, and simply duplicate the arrays
    if (verbose) { message("NoPSF: Convolved Apertures are identical to Simple Apertures")
                   message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    wsfa<-sfa
    image.env$wfa<-image.env$fa
  } else if ((!nopsf)&(length(which(fluxweight!=1))!=0)) {
    #If convolving with psf, and if all the fluxweights are not unity,
    #then colvolve and weight the stamps/apertures
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #Create a full mask of all convolved stamps/apertures
    timer=system.time(image.env$fa<-make_a_mask(environment(), sfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #If we have nonunity Fluxweights, then weight apertures
    #weight the stamps/apertures by the input fluxweights
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #Create a full mask of all the convolved weighted stamps/apertures
    timer=system.time(image.env$wfa<-make_a_mask(environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  } else if ((!nopsf)&(length(which(fluxweight!=1))==0)) {
    #If convolving with psf, and if all the fluxweights are unity,
    #then colvolve and weight the stamps/apertures
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #Create a full mask of all convolved stamps/apertures
    timer=system.time(image.env$fa<-make_a_mask(environment(), sfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #If we have unity Fluxweights, duplicate apertures
    if (verbose) { message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    wsfa<-sfa
    image.env$wfa<-image.env$fa
  }
  if (diagnostic) {
    if (verbose) { message("Checking Apertures for scaling errors") }
    foreach(i=1:length(sa), .combine=function(a,b){NULL},.inorder=FALSE) %dopar% {
      if (max(sa[[i]])>1){warning(paste("Max of Aperture",i,"is",max(sa[[i]]))) }
      if (max(sfa[[i]])>1){warning(paste("Max of Convolved Aperture",i,"is",max(sfa[[i]]))) }
      #if (max(wsfa[[i]])>1){warning(paste("Max of Weighted Convolved Aperture",i,"is",max(wsfa[[i]]))) }
    }
  }

  #Do we want to plot a sample of the apertures?
  if (plotsample) {
    #Set output name
    pdf(paste(pathout,"PSF_Samples.pdf",sep=""))
    #Set Layout
    par(mfrow=c(2,2))
    #Output a random 15 apertures to file
    for (i in (runif(1:15)*npos)) {
      image(sa[[i]], main="Single Aperture (SA)", asp=1, col=heat.colors(1000), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
      image(sfa[[i]], main="Single Convolved Apertrure (SFA)", asp=1, col=heat.colors(1000), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
    }
    #Close the file
    dev.off()
  }
  #Remove arrays that are no longer needed
  rm(sa, envir=env)
  gc()


  #Do we want to output the convolved apertures mask?
  if (makefamask) {
    if (!quiet) { cat(paste('Outputting All Convolved Apertures Mask to',fafilename,"   ")) }
    timer=system.time(writefitsout(paste(pathout,fafilename,sep=""),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
  #Do we want to make the SourceMask to file?
  if (sourcemask) {
    #Note:
    #> In the case of flux measurements: we want SM to be 1 only where the sources are, and within the mask region
    #> In the case of creating a source mask: we want SM to be 1 within the mask region in between sources only
    #
    #To do this we:
    #         a) set the sm array = image mask (0s outside image, 1s inside)
    #         b) set all nonzero locations in aa to 0 in sm array
    if (!quiet) { cat(paste('Outputting Source Mask to',smfilename,"   ")) }
    if (length(image.env$imm)!=1) { sm<-image.env$imm } else { sm<-array(1, dim=dim(image.env$im)) }
    sm[which(image.env$fa > height*(pnorm(-3)))]<-0
    #-----Diagnoistic-----#
    if (diagnostic) {
      message(paste("SourceMask Max/Min:",max(sm),min(sm)))
      message(paste("OLDMethod - SourceMask Max/Min:",max(1-image.env$fa),min(1-image.env$fa)))
    }
    #Are we outputting the sourcemask?
    if (!is.null(smfilename)){ writefitsout(paste(pathout,smfilename, sep=""),sm,image.env$hdr_str,nochange=TRUE) }
    if (!quiet) { cat(" - Done\n") }

    #If we want the SourceMask only, then end here
    if (sourcemaskonly) {
      if (!quiet) { cat("SourceMaskOnly Flag Set\n")  }
      return()
    }
    #Remove arrays that are no longer needed
    rm(fa, envir=image.env)
    #rm(sm, envir=env)
    #gc()
  }

  #---------------------------------------------------------------
  #PART TWO-A:  DEBLENDING
  #---------------------------------------------------------------

  #Perform Deblending of apertures
  timer=system.time(dbw<-make_deblended_weightmap(environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Deblended Weightmap - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }

  #Remove arrays that are no longer needed
  rm(wfa, envir=image.env)
  rm(wsfa, envir=env)
  gc()

  #Generate Deblended Flux Arrays
  if (!quiet) { cat("Generating Deblended Flux Arrays") }
  timer<-proc.time()

  ##-----Diagnostic-----#
  #if (diagnostic) {
  ##check dbw - all values of deblended array >=0
    #dbw_new<-foreach(i=1:npos) %dopar% {
      #if(length(which(zapsmall(dbw[[i]]) <0))>0)  {
        #warning("There are values of the deblended matrix < 0")
        #dbw[which(dbw[[i]] <0)]<-0
      #}
      #dbw[[i]]
    #}
  #dbw<-dbw_new
  #Remove array that is no longer needed
  #rm(dbw_new, envir=env)
  #}

  #Create the Deblended Flux Array
  #dfa<-foreach(i=1:npos) %dopar% { dbw[[i]]*fluxweight[i] }
  dfa<-foreach(i=1:npos) %dopar% { dbw[[i]]*sfa[[i]] }
  #dfa<-foreach(i=1:npos) %dopar% { dbw[[i]]*0.0 +1 }
  #dfa<-foreach(i=1:npos) %dopar% { wsfa[[i]] }

  if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Make DFA - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   - Done\n") }

  #Do we want to output the deblended convolved apertures mask?
  if (makedfamask) {
    if (!quiet) { cat(paste('Making All Deblended Convolved Apertures Mask - ')) }
    timer=system.time(image.env$adfa<-make_a_mask(environment(), dfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (!quiet) { cat(paste('Outputting All Deblended Convolved Apertures Mask to',dfafilename,"   ")) }
    timer=system.time(writefitsout(paste(pathout,dfafilename,sep=""),image.env$adfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }

  #Remove array that is no longer needed
  rm(dbw, envir=env)
  rm(adfa, envir=image.env)
  gc()

  #---------------------------------------------------------------
  #PART THREE: COMPUTE AND OUTPUT GALAXY-BY-GALAXY RESULTS
  #---------------------------------------------------------------

  # Make Images Available
  attach(image.env)

  if (!quiet) { cat("Computing Galaxy-by-Galaxy Results {\n") }

  #-----Diagnostic-----#
  if (diagnostic){
    if (length(which(is.na(im)))!=0) {        message(paste("Input Parameter 'im' contains NA elements")) }
    if (length(which(is.na(ime)))!=0) {        message(paste("Input Parameter 'ime' contains NA elements")) }
    if (length(which(is.na(sfa)))!=0) {        message(paste("Input Parameter 'sfa' contains NA elements")) }
    if (length(which(is.na(dfa)))!=0) {        message(paste("Input Parameter 'dfa' contains NA elements")) }
    if (length(which(is.na(fluxcorr)))!=0) {        message(paste("Input Parameter 'fluxcorr' contains NA elements")) }

  }
  if (!quiet) { cat("   Calculating Fluxes {  \n") }
  timer<-proc.time()

  #Compute Galaxy-by-Galaxy Results
#-----
  #Image Flux at central pixel
  if (verbose) { cat("      Image Flux at central pixel ") }
  pixflux<-foreach(i=1:npos, .inorder=TRUE)%dopar% {
    if ((!(is.finite(x_p[i]) & is.finite(y_p[i]))) | (x_p[i] <= 0) | (y_p[i] <= 0) |
        (x_p[i] > (length(im[,1]))) | (y_p[i] > (length(im[1,])))) { -99 } else { im[x_p[i],y_p[i]]  }
  }
  pixflux<-array(unlist(pixflux),dim=c(dim(pixflux[[1]]),length(pixflux)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the aperture
  if (verbose) { cat("      Integral of the aperture") }
  ssa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(sa[[i]]) }
  ssa<-array(unlist(ssa),dim=c(dim(ssa[[1]]),length(ssa)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the convolved aperture
  if (verbose) { cat("      Integral of the convolved aperture") }
  ssfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(sfa[[i]]) }
  ssfa<-array(unlist(ssfa),dim=c(dim(ssfa[[1]]),length(ssfa)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(convolved aperture)^2]
  if (verbose) { cat("      Integral of the [(convolved aperture)^2]") }
  ssfa2<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum((sfa[[i]])^2.)  }
  ssfa2<-array(unlist(ssfa2),dim=c(dim(ssfa2[[1]]),length(ssfa2)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the (convolved aperture * image)
  if (verbose) { cat("      Integral of the (convolved aperture * image)") }
  ssfad<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(sfa[[i]]*(im[xlo:xup,ylo:yup]))
  }
  ssfad<-array(unlist(ssfad),dim=c(dim(ssfad[[1]]),length(ssfad)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the deblended convolved aperture
  if (verbose) { cat("      Integral of the deblended convolved aperture") }
  sdfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(dfa[[i]])  }
  sdfa<-array(unlist(sdfa),dim=c(dim(sdfa[[1]]),length(sdfa)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(deblended convolved aperture)^2]
  if (verbose) { cat("      Integral of the [(deblended convolved aperture)^2]") }
  sdfa2<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum((dfa[[i]])^2.)  }
  sdfa2<-array(unlist(sdfa2),dim=c(dim(sdfa2[[1]]),length(sdfa2)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the (deblended convolved aperture * image)
  if (verbose) { cat("      Integral of the (deblended convolved aperture * image)") }
  sdfad<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(dfa[[i]]*(im[xlo:xup,ylo:yup]))
  }
  sdfad<-array(unlist(sdfad),dim=c(dim(sdfad[[1]]),length(sdfad)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the (deblended convolved aperture * convolved aperture)
  if (verbose) { cat("      Integral of the (deblended convolved aperture * convolved aperture)") }
  sdfasfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(dfa[[i]]*sfa[[i]]) }
  sdfasfa<-array(unlist(sdfasfa),dim=c(dim(sdfasfa[[1]]),length(sdfasfa)))
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the (convolved aperture * image error)
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime)!=1) {
    ssfae<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(sfa[[i]]*ime[xlo:xup,ylo:yup])
    }
    ssfae<-array(unlist(ssfae),dim=c(dim(ssfae[[1]]),length(ssfae)))
  } else if (ime==1) {
    ssfae<-ssfa
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the (convolved aperture * image error)
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime)!=1) {
    ssfae2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]*(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfae2<-array(unlist(ssfae2),dim=c(dim(ssfae2[[1]]),length(ssfae2)))
  } else if (ime==1) {
    ssfae2<-ssfa
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(convolved aperture * image error)^2]
  if (verbose) { cat("      Integral of the [(convolved aperture * image error)^2]") }
  if (length(ime)!=1) {
    ssfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((sfa[[i]]*ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2e2<-array(unlist(ssfa2e2),dim=c(dim(ssfa2e2[[1]]),length(ssfa2e2)))
  } else if (ime==1) {
    ssfa2e2<-ssfa2
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(deblended convolved aperture * image error)^2]
  if (verbose) { cat("      Integral of the [(deblended convolved aperture * image error)^2]") }
  if (length(ime)!=1) {
    sdfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((dfa[[i]]*ime[xlo:xup,ylo:yup])^2.)
    }
    sdfa2e2<-array(unlist(sdfa2e2),dim=c(dim(sdfa2e2[[1]]),length(sdfa2e2)))
  } else if (ime==1) {
    sdfa2e2<-sdfa2
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(convolved aperture * image) / {(image error)^2} ]
  if (verbose) { cat("      Integral of the [(convolved aperture * image) / {(image error)^2} ]") }
  if (length(ime)!=1) {
    ssfadw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]*(im[xlo:xup,ylo:yup])/((ime[xlo:xup,ylo:yup])^2.))
    }
    ssfadw<-array(unlist(ssfadw),dim=c(dim(ssfadw[[1]]),length(ssfadw)))
  } else if (ime==1) {
    ssfadw<-ssfad
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [convolved aperture / (image error)^2 ]
  if (verbose) { cat("      Integral of the [convolved aperture / (image error)^2 ]") }
  if (length(ime)!=1) {
    ssfaw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]/(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfaw<-array(unlist(ssfaw),dim=c(dim(ssfaw[[1]]),length(ssfaw)))
  } else if (ime==1) {
    ssfaw<-ssfa
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  #Integral of the [(convolved aperture / image error)^2]
  if (verbose) { cat("      Integral of the [(convolved aperture / image error)^2]") }
  if (length(ime)!=1) {
    ssfa2w<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((sfa[[i]]/ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2w<-array(unlist(ssfa2w),dim=c(dim(ssfa2w[[1]]),length(ssfa2w)))
  } else if (ime==1) {
    ssfa2w<-ssfa2
  } else {
    sink(type="message")
    stop("Image Error is neither an array, nor equal to unity")
  }
  if (verbose) { cat(" - Done\n") }
#-----
  if (showtime) { cat("   } - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Perform Calculations - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   } - Done\n") }

#
# NB: in the following calculations - Int(<exp>) denotes the integral of the expression <exp>
#

  #Convolved aperture flux = Int(fltAp*Im) * Int(fltAp) / Int(fltAp^2)
  sfaflux<-ssfad*ssfa/ssfa2

  #Deblended convolved aperture flux = Int(DBfltAp*Im) * Int(fltAp) / Int(fltAp^2)
  dfaflux<-sdfad*ssfa/ssfa2

  sfaerr<-sqrt((ssfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(ssfa)))

  dfaerr<-sqrt((sdfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(sdfa)))
  sfafluxw<-ssfadw*ssfa/ssfa2w
  sfaerrw<-sqrt(1/ssfaw*(ssfa/beamarea*ssfa/ssfa2)^2. + ((conf*beamarea)^2.*sqrt(ssfa)))

  #Do we want to do sky estimation/subtraction?
  if (doskyest) {
    #Get sky estimates
    if (verbose) { message("Perfoming Sky Subtraction"); cat("   Performing Sky Subtraction") }
    skyest<-skyback(ra_g,dec_g,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,origim=list(dat=list(im)),maskim=list(dat=list(sm)),astrom=astr_struc)
    skyflux<-skyest[,1]
    skyerr<-skyest[,2]
    #browser()
    #Subrtract Sky Flux
    skyflux<-skyflux*sdfa
    skyerr<-skyerr*sdfa
    dfaflux<-dfaflux-skyflux
    sfaflux<-sfaflux-skyflux
    dfaerr<-sqrt(dfaerr^2+skyerr^2)
    dfaerr<-sqrt(sfaerr^2+skyerr^2)
    if (verbose) { message(paste("   - Done\n")); cat("   - Done\n")}
  } else {
    skyflux<-array(0, dim=c(length(dfaflux)))
    skyerr<-array(0, dim=c(length(dfaflux)))
  }

  if (!quiet) { cat("   Performing Final Calculations   ") }

  #Apply Flux Correction to the Finalised Values
  sfaflux<-sfaflux*fluxcorr
  sfafluxw<-sfafluxw*fluxcorr
  dfaflux<-dfaflux*fluxcorr
  sfaerr<-sfaerr*fluxcorr
  sfaerrw<-sfaerrw*fluxcorr
  dfaerr<-dfaerr*fluxcorr

  #Calculate Magnitudes
  if (Magnitudes) {
    mags<--2.5*(log10(dfaflux)-log10(ABvegaflux))+magZP
  } else {
    mags<-rep(NA, length(id_g))
  }

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(ssa    )))/length(ssa    )*100,digits=2),"% of the ssa matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfa   )))/length(ssfa   )*100,digits=2),"% of the ssfa matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sdfa   )))/length(sdfa   )*100,digits=2),"% of the sdfa matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfa2  )))/length(ssfa2  )*100,digits=2),"% of the ssfa2 matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sdfa2  )))/length(sdfa2  )*100,digits=2),"% of the sdfa2 matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sdfasfa)))/length(sdfasfa)*100,digits=2),"% of the sdfasfa matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfad  )))/length(ssfad  )*100,digits=2),"% of the ssfad matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sdfad  )))/length(sdfad  )*100,digits=2),"% of the sdfad matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sfaflux)))/length(sfaflux)*100,digits=2),"% of the sfaflux matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sfaerr )))/length(sfaerr )*100,digits=2),"% of the sfaerr matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(dfaflux)))/length(dfaflux)*100,digits=2),"% of the dfaflux matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(dfaerr )))/length(dfaerr )*100,digits=2),"% of the dfaerr matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfae  )))/length(ssfae  )*100,digits=2),"% of the ssfae matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfae2 )))/length(ssfae2 )*100,digits=2),"% of the ssfae2 matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfa2e2)))/length(ssfa2e2)*100,digits=2),"% of the ssfa2e2 matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(sdfa2e2)))/length(sdfa2e2)*100,digits=2),"% of the sdfa2e2 matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfadw )))/length(ssfadw )*100,digits=2),"% of the ssfadw matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfaw  )))/length(ssfaw  )*100,digits=2),"% of the ssfaw matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(ssfa2w )))/length(ssfa2w )*100,digits=2),"% of the ssfa2w matrix are NA"))
    message(paste("After assignment",round(length(which(is.na(pixflux)))/length(pixflux)*100,digits=2),"% of the pixflux matrix are NA"))
  }

  #Check that final values of Deblended Convolved Apertures are not NA/NaN/Inf
  if (length(which(!is.finite(dfaerr[which(sdfa > 0)]))) > 0) {
    message(paste(length(!is.finite(dfaerr[which(sdfa > 0)])), "elements of dfaerr are not finite"))
    browser()
    sink(type="message")
    stop("NaN or Infs Produced in calculations")
  }

  if (!quiet) { cat(" - Done\n} Galaxy Results Complete\n") }

  #---------------------------------------------------------------
  #PART FOUR: OUTPUT
  #---------------------------------------------------------------

  #If map was input in Jy/bm we need to convert it back before output in SourceSubtraction
  if (Jybm) { ba=beamarea } else { ba=1. }
  # Do we want to overlay the elipses?
  if (overlay) {
    if (!quiet) { cat(paste("Writing Ellipse-Overlaid Image Map to",overlaymap,"   ")) }
    #Make overlay stamps (where aperture=FWHM, pix=1E3, else == 0)
    overlaystamps<-foreach(i=1:length(sfa)) %dopar% {
      tmp=sfa[[i]]
      tmp[which(tmp<(max(tmp)*2/5))]=0
      tmp[which(tmp>(max(tmp)*3/5))]=0
      tmp[which(tmp>0)]=-1E3
      tmp
    }
    #Perform Overlay
    timer=system.time(sourcesubtraction(im,overlaystamps,stamp_lims,1,paste(pathout,overlaymap, sep=""),hdr_str,ba,insidemask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Contam Subtraction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
  # Do we want to make the residual map?
  if (makeresidmap) {
    if (filtcontam) {
      if (!quiet) { cat(paste("Writing Contaminant-subtracted Map to",nocontammap,"   ")) }
      #Perform Source Subtraction
      timer=system.time(sourcesubtraction(im,sfa,stamp_lims,dfaflux,paste(pathout,nocontammap, sep=""),hdr_str,ba,contams))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Contam Subtraction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
    if (!quiet) { cat(paste("Writing Source-subtracted Map to",residmap,"   ")) }
    #Perform Source Subtraction
    timer=system.time(sourcesubtraction(im,sfa,stamp_lims,dfaflux,paste(pathout,residmap, sep=""),hdr_str,ba,insidemask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Source Subtraction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }

  #If Subtracting Contaminants, then remove them before output
  if ((writetab)&(filtcontam)) {
    id_g   <-id_g[which(contams==0)]
    ra_g   <-ra_g[which(contams==0)]
    dec_g  <-dec_g[which(contams==0)]
    theta_g<-theta_g[which(contams==0)]
    a_g    <-a_g[which(contams==0)]
    b_g    <-b_g[which(contams==0)]
    ssa    <-ssa[which(contams==0)]
    ssfa   <-ssfa[which(contams==0)]
    ssfa2  <-ssfa2[which(contams==0)]
    ssfad  <-ssfad[which(contams==0)]
    ssfa2e2<-ssfa2e2[which(contams==0)]
    sfaflux<-sfaflux[which(contams==0)]
    skyflux<-skyflux[which(contams==0)]
    skyerr <-skyerr[which(contams==0)]
    sfaerr <-sfaerr[which(contams==0)]
    sdfa   <-sdfa[which(contams==0)]
    sdfa2  <-sdfa2[which(contams==0)]
    sdfad  <-sdfad[which(contams==0)]
    sdfa2e2<-sdfa2e2[which(contams==0)]
    dfaflux<-dfaflux[which(contams==0)]
    dfaerr <-dfaerr[which(contams==0)]
    pixflux<-pixflux[which(contams==0)]
    mags<-mags[which(contams==0)]
  }

  #Do we want to output the Results Table?
  if (writetab) {
    if (!quiet) { cat(paste('Writing Results Table to ',tableoutname,'   ')) }
    #Output the results table
    timer=system.time(writesfatableout(environment(), paste(pathout,tableoutname, sep="")) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Write Results Table - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }

  #-----Diagnostic-----#
  if (diagnostic) { message(paste('Sum of PSF = ',beamarea)) }

  if (interact) {
     cat(paste("Launching Interactive Mode: To end, type 'c'\n"))
     sink(type='message')
     browser()
     sink(sinkfile, type='message')
  }
  if (!quiet) { cat('\n') }
  #Return Output Variables
  detach(env)
  return=list(SFAflux=sfaflux,SFAerror=sfaerr, SFAfluxWeight=sfafluxw, SFAerrorWeight=sfaerrw, DFAflux=dfaflux,DFAerror=dfaerr)
}
