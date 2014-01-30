fluxmeasurements <-
function(env=NULL) {
#Procedure measures the fluxes present in a supplied image
#inside catalogued apertures.

  #PART ZERO:  PSF DETERMINATION & GENERAL INITIALISATION {{{
  # Load Parameter Space {{{
  if (is.null(env)) {
    stop("No Parameter Space Environment Specified in Call")
  }
  attach(env)
  #}}}

  if (!quiet) { cat('Begining Flux Measurements\n') }
  message('} Initialisation Complete\nBegining Flux Measurements')

  #Get PSF details {{{
  timer<-proc.time()
  if (!quiet) { cat('Getting PSF details') }
  message('Getting PSF details')
  if (!(nopsf)) {
    #We are filtering apertures with the PSF {{{
    #Details {{{
    #Calculate PSF - if one has not be supplied,
    #then a gaussian PSF will be created. Stampsize of PSF
    #should be = maximum stampsize: max aperture * stamp mult }}}
    #Get PSF {{{
    psf<-readpsf(environment(),paste(pathroot,psfmap,sep=""),asperpix,defbuff*max(a_g),confidence,gauss_fwhm_as=gauss_fwhm_as)
    #}}}
    #Notify {{{
    if (verbose) { message(paste("Maxima of the PSF is at pixel", which(psf == max(psf)),"and has value",max(psf))) }
    #}}}
    #Normalise Beam Area {{{
    beamarea_nn<-sumpsf
    beamarea_n<-as.single(sum(psf))
    #}}}
    #-----Diagnostic-----##{{{
    if (diagnostic) { message(paste('Beam area before/after norm: ',beamarea_nn,beamarea_n)) }
    #}}}
    #}}}
  } else {
    #We are not using PSF filtering {{{
    #set relevant paramters manually
    beamarea_n<-1.0
    psfwidth<-0
    psf.clip<-0
    #}}}
  }
  #If possible, Update Beamarea {{{
  if (beamarea_pix == 0) {
    #Use the beamarea just determined {{{
    beamarea<-beamarea_n
    #}}}
  } else {
    # Otherwise just use the input beamarea {{{
    beamarea<-beamarea_pix
    #}}}
  }#}}}
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) { message(paste('Beam area adopted: ',beamarea)) }
  #}}}

  #Convert images from Jy/beam to Jy/pixel {{{
  if (Jybm) {
    image.env$im<-image.env$im/beamarea
    if (length(image.env$ime)!=1) { image.env$ime<-image.env$ime/beamarea }
    conf<-conf/beamarea
  } #}}}

  #Notify {{{
  if (showtime) { cat(paste(' - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )\n'))
    message(paste('Getting PSF details - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat(' - Done\n') }
  #}}}
  #}}}

  #PART ONE:  SINGLE APERTURE MASKS {{{
  #Details {{{
  #Create apertures for every object remaining in
  #the catalogue. }}}

  #Convert decimal pixel values into actual pixel values {{{
  x_p<-floor(x_g)
  y_p<-floor(y_g)
  #}}}

  #Create an array of stamps containing the apertures for all objects {{{
  timer=system.time(sa<-make_sa_mask(environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make SA Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}

  #Discard any apertures that were zero'd in the make process {{{
  totsa<-foreach(i=1:length(sa), .combine='c') %dopar% { if (max(sa[[i]])>1){warning(paste("Max of Aperture",i,"is",max(sa[[i]]))) } ; sum(sa[[i]]) }
  insidemask<-(totsa > 0)
  if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Single Apertures are inside the Mask.") }  # Nothing inside the mask
  sa<-sa[which(insidemask)]
  stamplen<-stamplen[which(insidemask)]
  #Check if there is more than 1 aperture left (array dimension errors arise if not) {{{
  if (length(which(insidemask))==1) {
    stamp_lims<-rbind(stamp_lims[which(insidemask),],stamp_lims[which(insidemask),])
  } else {
    stamp_lims<-stamp_lims[which(insidemask),]
  } #}}}
  #Remove object that failed aperture creation {{{
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
  #}}}
  #}}}

  #Re-Initialise object count {{{
  npos<-length(id_g)
  #}}}

  #Create a full mask of all apertures in their correct image-space locations {{{
  timer=system.time(image.env$aa<-make_a_mask(environment(), sa, dim(image.env$im)))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
   message(paste('Make A Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}

  #If wanted, output the All Apertures Mask {{{
  if (makeaamask) {
    if (!quiet) { cat(paste('Outputting All Apertures Mask to',aafilename,"   ")) }
    #Write All Apertures Mask to file
    writefitsout(paste(pathout,aafilename,sep=""),image.env$aa,image.env$hdr_str,nochange=TRUE)
    if (!quiet) { cat(" - Done\n") }
  }#}}}
  #}}}

  #PART TWO:  SINGLE FILTERED APERTURE MASKS {{{
  #Details {{{
  #If wanted, perform the convolution of the
  #apertures with the PSF. If not wanted,
  #duplicate the single apertures.
  #Also, perform flux weighting. }}}

  #If needed, do Convolutions & Fluxweightings {{{
  if ((nopsf)&(length(which(fluxweight!=1))!=0)) {
    #No Convolution, Need Fluxweighting {{{
    #Details {{{
    #If not convolving with psf, and if all the fluxweights are not unity,
    #then skip filtering, duplicate the arrays, and only weight the stamps/apertures
    #}}}
    #No Convolution, duplicate apertures {{{
    if (verbose) { message("NoPSF: Convolved Apertures are identical to Simple Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    #}}}
    #Perform aperture weighting {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of stamps/apertures {{{
    timer=system.time(image.env$wfa<-make_a_mask(environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #}}}
  } else if ((nopsf)&(length(which(fluxweight!=1))==0)) {
    #No Convolution, No Fluxweighting {{{
    #Details {{{
    #If not convolving with psf, and if all the fluxweights are unity,
    #then skip filtering and weighting, and simply duplicate the arrays }}}
    #Duplicate arrays {{{
    if (verbose) { message("NoPSF: Convolved Apertures are identical to Simple Apertures")
                   message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    wsfa<-sfa
    image.env$wfa<-image.env$fa
    #}}}
    #}}}
  } else if ((!nopsf)&(length(which(fluxweight!=1))!=0)) {
    #Convolving PSF, Need Fluxweighting {{{
    #Details {{{
    #If convolving with psf, and if all the fluxweights are not unity,
    #then colvolve and weight the stamps/apertures }}}
    #Perform Convolution {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all convolved stamps/apertures {{{
    timer=system.time(image.env$fa<-make_a_mask(environment(), sfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Perform aperture weighting {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all the convolved weighted stamps/apertures {{{
    timer=system.time(image.env$wfa<-make_a_mask(environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #}}}
  } else if ((!nopsf)&(length(which(fluxweight!=1))==0)) {
    #Convolve PSF, No Fluxweighting {{{
    #Details {{{
    #If convolving with psf, and if all the fluxweights are unity,
    #colvolve, and then duplicate the stamps/apertures }}}
    #Perform Convolution {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all convolved stamps/apertures {{{
    timer=system.time(image.env$fa<-make_a_mask(environment(), sfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #No Fluxweights, duplicate apertures {{{
    if (verbose) { message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    wsfa<-sfa
    image.env$wfa<-image.env$fa
    #}}}
    #}}}
  }#}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    if (verbose) { message("Checking Apertures for scaling errors") }
    foreach(i=1:length(sa), .combine=function(a,b){NULL},.inorder=FALSE) %dopar% {
      if (max(sa[[i]])>1){warning(paste("Max of Aperture",i,"is",max(sa[[i]]))) }
      if (max(sfa[[i]])>1){warning(paste("Max of Convolved Aperture",i,"is",max(sfa[[i]]))) }
    }
  }#}}}

  #Do we want to plot a sample of the apertures? {{{
  if (plotsample) {
    #Set output name {{{
    pdf(paste(pathout,"PSF_Samples.pdf",sep=""))
    #}}}
    #Set Layout {{{
    par(mfrow=c(2,2))
    #}}}
    #Output a random 15 apertures to file {{{
    for (i in (runif(1:15)*npos)) {
      image(sa[[i]], main="Single Aperture (SA)", asp=1, col=heat.colors(1000), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
      image(sfa[[i]], main="Single Convolved Apertrure (SFA)", asp=1, col=heat.colors(1000), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
    }#}}}
    #Close the file {{{
    dev.off()
    #}}}
  }
  #}}}

  #Remove arrays that are no longer needed {{{
  rm(sa, envir=env)
  gc()
  #}}}

  #If wanted, output the Convolved & Weighted Aperture Mask {{{
  if (makefamask) {
    if (!quiet) { cat(paste('Outputting All Convolved Apertures Mask to',fafilename,"   ")) }
    timer=system.time(writefitsout(paste(pathout,fafilename,sep=""),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }#}}}

  #If wanted, create (& possibly output) the SourceMask {{{
  if (sourcemask) {
    #Details {{{
    #> In the case of flux measurements: we want SM to be 1 only where the sources are, and within the mask region
    #> In the case of creating a source mask: we want SM to be 1 within the mask region in between sources only
    #
    #To do this we:
    #         a) set the sm array = image mask (0s outside image, 1s inside)
    #         b) set all nonzero locations in aa to 0 in sm array }}}
    #Make Sourcemask {{{
    if (!quiet) { cat(paste("Creating Sourcemask    ")) }
    if (length(image.env$imm)!=1) { sm<-image.env$imm } else { sm<-array(1, dim=dim(image.env$im)) }
    if (nopsf) {
      sm[which(image.env$fa > 0)]<-0
    } else {
      sm[which(image.env$fa > get.confidence(psf,smConfidenceLim,value=TRUE))]<-0
    }
    if (!quiet) { cat(" - Done\n") }
    #}}}
    #-----Diagnoistic-----# {{{
    if (diagnostic) {
      message(paste("SourceMask Max/Min:",max(sm),min(sm)))
      message(paste("OLDMethod - SourceMask Max/Min:",max(1-image.env$fa),min(1-image.env$fa)))
    }#}}}
    #If wanted, output the SourceMask {{{
    if (!is.null(smfilename)){
      if (!quiet) { cat(paste('Outputting Source Mask to',smfilename,"   ")) }
      writefitsout(paste(pathout,smfilename, sep=""),sm,image.env$hdr_str,nochange=TRUE)
      if (!quiet) { cat(" - Done\n") }
    }#}}}

    #If we want the SourceMask only, then end here {{{
    if (sourcemaskonly) {
      if (!quiet) { cat("SourceMaskOnly Flag Set\n")  }
      return()
    } #}}}

    #Remove arrays that are no longer needed {{{
    rm(fa, envir=image.env)
    #}}}
  }#}}}
  #}}}

  #PART THREE:  DEBLENDING {{{

  #Perform Deblending of apertures {{{
  timer=system.time(dbw<-make_deblended_weightmap(environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Deblended Weightmap - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}

  #Remove arrays that are no longer needed {{{
  rm(wfa, envir=image.env)
  rm(wsfa, envir=env)
  gc()
  #}}}

  #Generate Deblended Flux Arrays {{{
  if (!quiet) { cat("Generating Deblended Flux Arrays") }
  timer<-proc.time()
  #Create the Deblended Flux Array {{{
  dfa<-foreach(i=1:npos) %dopar% { dbw[[i]]*sfa[[i]] }
  #}}}
  #Notify {{{
  if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Make DFA - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #}}}

  #If wanted, output the Deblended & Convolved Aperture Mask {{{
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
  }#}}}

  #Remove array that is no longer needed {{{
  rm(dbw, envir=env)
  rm(adfa, envir=image.env)
  gc()
  #}}}
  #}}}

  #PART FOUR: COMPUTE AND OUTPUT GALAXY-BY-GALAXY RESULTS {{{

  if (!quiet) { cat("Computing Galaxy-by-Galaxy Results {\n") }
  # Make Images Available {{{
  attach(image.env)
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic){
    if (length(which(is.na(im)))!=0) {        message(paste("Input Parameter 'im' contains NA elements")) }
    if (length(which(is.na(ime)))!=0) {        message(paste("Input Parameter 'ime' contains NA elements")) }
    if (length(which(is.na(sfa)))!=0) {        message(paste("Input Parameter 'sfa' contains NA elements")) }
    if (length(which(is.na(dfa)))!=0) {        message(paste("Input Parameter 'dfa' contains NA elements")) }
    if (length(which(is.na(fluxcorr)))!=0) {        message(paste("Input Parameter 'fluxcorr' contains NA elements")) }
  }#}}}

  #Perform Calculations {{{
  if (!quiet) { cat("   Calculating Fluxes {  \n") }
  timer<-proc.time()
#-----
  #Image Flux at central pixel {{{
  if (verbose) { cat("      Image Flux at central pixel ") }
  pixflux<-foreach(i=1:npos, .inorder=TRUE)%dopar% {
    if ((!(is.finite(x_p[i]) & is.finite(y_p[i]))) | (x_p[i] <= 0) | (y_p[i] <= 0) |
        (x_p[i] > (length(im[,1]))) | (y_p[i] > (length(im[1,])))) { -99 } else { im[x_p[i],y_p[i]]  }
  }
  pixflux<-array(unlist(pixflux),dim=c(dim(pixflux[[1]]),length(pixflux)))
  if (verbose) { cat(" - Done\n") }#}}}
#-----
  #Integral of the aperture {{{
  if (verbose) { cat("      Integral of the aperture") }
  ssa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(sa[[i]]) }
  ssa<-array(unlist(ssa),dim=c(dim(ssa[[1]]),length(ssa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the convolved aperture {{{
  if (verbose) { cat("      Integral of the convolved aperture") }
  ssfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(sfa[[i]]) }
  ssfa<-array(unlist(ssfa),dim=c(dim(ssfa[[1]]),length(ssfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture)^2] {{{
  if (verbose) { cat("      Integral of the [(convolved aperture)^2]") }
  ssfa2<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum((sfa[[i]])^2.)  }
  ssfa2<-array(unlist(ssfa2),dim=c(dim(ssfa2[[1]]),length(ssfa2)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * image) {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image)") }
  ssfad<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(sfa[[i]]*(im[xlo:xup,ylo:yup]))
  }
  ssfad<-array(unlist(ssfad),dim=c(dim(ssfad[[1]]),length(ssfad)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the deblended convolved aperture {{{
  if (verbose) { cat("      Integral of the deblended convolved aperture") }
  sdfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(dfa[[i]])  }
  sdfa<-array(unlist(sdfa),dim=c(dim(sdfa[[1]]),length(sdfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(deblended convolved aperture)^2] {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture)^2]") }
  sdfa2<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum((dfa[[i]])^2.)  }
  sdfa2<-array(unlist(sdfa2),dim=c(dim(sdfa2[[1]]),length(sdfa2)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (deblended convolved aperture * image) {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * image)") }
  sdfad<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(dfa[[i]]*(im[xlo:xup,ylo:yup]))
  }
  sdfad<-array(unlist(sdfad),dim=c(dim(sdfad[[1]]),length(sdfad)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (deblended convolved aperture * convolved aperture) {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * convolved aperture)") }
  sdfasfa<-foreach(i=1:npos, .inorder=TRUE)%dopar% { sum(dfa[[i]]*sfa[[i]]) }
  sdfasfa<-array(unlist(sdfasfa),dim=c(dim(sdfasfa[[1]]),length(sdfasfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * image error) {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime)!=1) {
    ssfae<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      sum(sfa[[i]]*ime[xlo:xup,ylo:yup])
    }
    ssfae<-array(unlist(ssfae),dim=c(dim(ssfae[[1]]),length(ssfae)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfae<-ssfa
    } else {
      ssfae<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]*ime)
      }
      ssfae<-array(unlist(ssfae),dim=c(dim(ssfae[[1]]),length(ssfae)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * image error) {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime)!=1) {
    ssfae2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]*(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfae2<-array(unlist(ssfae2),dim=c(dim(ssfae2[[1]]),length(ssfae2)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfae2<-ssfa
    } else {
      ssfae2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum(sfa[[i]]*(ime)^2.)
      }
      ssfae2<-array(unlist(ssfae2),dim=c(dim(ssfae2[[1]]),length(ssfae2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture * image error)^2] {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image error)^2]") }
  if (length(ime)!=1) {
    ssfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((sfa[[i]]*ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2e2<-array(unlist(ssfa2e2),dim=c(dim(ssfa2e2[[1]]),length(ssfa2e2)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfa2e2<-ssfa2
    } else {
      ssfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum((sfa[[i]]*ime)^2.)
      }
      ssfa2e2<-array(unlist(ssfa2e2),dim=c(dim(ssfa2e2[[1]]),length(ssfa2e2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(deblended convolved aperture * image error)^2] {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture * image error)^2]") }
  if (length(ime)!=1) {
    sdfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((dfa[[i]]*ime[xlo:xup,ylo:yup])^2.)
    }
    sdfa2e2<-array(unlist(sdfa2e2),dim=c(dim(sdfa2e2[[1]]),length(sdfa2e2)))
  } else if (length(ime)==1){
    if (ime==1) {
      sdfa2e2<-sdfa2
    } else {
      sdfa2e2<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum((dfa[[i]]*ime)^2.)
      }
      sdfa2e2<-array(unlist(sdfa2e2),dim=c(dim(sdfa2e2[[1]]),length(sdfa2e2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture * image) / {(image error)^2} ] {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image) / {(image error)^2} ]") }
  if (length(ime)!=1) {
    ssfadw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]*(im[xlo:xup,ylo:yup])/((ime[xlo:xup,ylo:yup])^2.))
    }
    ssfadw<-array(unlist(ssfadw),dim=c(dim(ssfadw[[1]]),length(ssfadw)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfadw<-ssfad
    } else {
      ssfadw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum(sfa[[i]]*(im[xlo:xup,ylo:yup])/((ime)^2.))
      }
      ssfadw<-array(unlist(ssfadw),dim=c(dim(ssfadw[[1]]),length(ssfadw)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [convolved aperture / (image error)^2 ] {{{
  if (verbose) { cat("      Integral of the [convolved aperture / (image error)^2 ]") }
  if (length(ime)!=1) {
    ssfaw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum(sfa[[i]]/(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfaw<-array(unlist(ssfaw),dim=c(dim(ssfaw[[1]]),length(ssfaw)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfaw<-ssfa
    } else {
      ssfaw<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum(sfa[[i]]/(ime)^2.)
      }
      ssfaw<-array(unlist(ssfaw),dim=c(dim(ssfaw[[1]]),length(ssfaw)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture / image error)^2] {{{
  if (verbose) { cat("      Integral of the [(convolved aperture / image error)^2]") }
  if (length(ime)!=1) {
    ssfa2w<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
        sum((sfa[[i]]/ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2w<-array(unlist(ssfa2w),dim=c(dim(ssfa2w[[1]]),length(ssfa2w)))
  } else if (length(ime)==1){
    if (ime==1) {
      ssfa2w<-ssfa2
    } else {
      ssfa2w<-foreach(i=1:npos,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
          sum((sfa[[i]]/ime)^2.)
      }
      ssfa2w<-array(unlist(ssfa2w),dim=c(dim(ssfa2w[[1]]),length(ssfa2w)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Final Flux & Error Calculations {{{
  if (verbose) { cat("      Final Fluxes and Error Calculations ") }
  #Convolved aperture flux = Int(fltAp*Im) * Int(fltAp) / Int(fltAp^2) {{{
  sfaflux<-ssfad*ssfa/ssfa2
  #}}}

  #Deblended convolved aperture flux = Int(DBfltAp*Im) * Int(fltAp) / Int(fltAp^2) {{{
  dfaflux<-sdfad*ssfa/ssfa2
  #}}}

  #Convolved aperture error {{{
  sfaerr<-sqrt((ssfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(ssfa)))
  #}}}

  #Deblended Convolved aperture error {{{
  dfaerr<-sqrt((sdfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(sdfa)))
  #}}}

  #Convolved Aperture Flux Weights & Error Weights {{{
  sfafluxw<-ssfadw*ssfa/ssfa2w
  sfaerrw<-sqrt(1/ssfaw*(ssfa/beamarea*ssfa/ssfa2)^2. + ((conf*beamarea)^2.*sqrt(ssfa)))
  #}}}
  if (verbose) { cat(" - Done\n") }
  #}}}
  if (showtime) { cat("   } - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Perform Calculations - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   } - Done\n") }
  #}}}

  #If wanted, perform sky estimation. Otherwise set to NA {{{
  if (doskyest||getskyrms) {
    #Get sky estimates {{{
    if (!quiet) { message("Perfoming Sky Estimation"); cat("   Performing Sky Estimation") }
    #Get PSF FWHM {{{
    if (nopsf) {
      #If no PSF, set FWHM to median aperture radius + buffer {{{
      psffwhm<-round(median(a_g[which(a_g>0)])*defbuff/asperpix)
      if (is.na(psffwhm)) {
        #All Apertures are point sources, force width = 5 pixels {{{
        #Details {{{
        #There are no apertures that are not point sources,
        #and no PSF convolution - everything is a single pixel flux
        #set to a default of 5pixels }}}
        psffwhm<-5
        #}}}
      }
      #}}}
    } else {
      #Get radius of FWHM using FWHM confidence value = erf(2*sqrt(2*log(2))/sqrt(2)) {{{
      psffwhm<-get.confidence(psf, confidence=(2*pnorm(2*sqrt(2*log(2)))-1))
      #}}}
    }#}}}
    #Perform Sky Estimation {{{
    timer<-system.time(skyest<-skyback(ra_g,dec_g,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,origim=list(dat=list(im)),maskim=list(dat=list(sm)),
                    astrom=astr_struc,clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm))
    #}}}
    #Get sky parameters {{{
    skylocal<-skyest[,'sky']
    skyerr<-skyest[,'skyerr']*correl.noise
    skyrms<-skyest[,'skyRMS']
    skypval<-skyest[,'skyRMSpval']
    skyflux<-skylocal*sdfa
    skyerr<-skyerr*sdfa
    #}}}
    #}}}
    #Notify {{{
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Sky Estimate - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Calculate Detection Thresholds {{{
    if (!quiet) { message("Calculating Detection Limits"); cat("   Calculating Detection Limits") }
    detecthres<-foreach(i=1:npos, .inorder=TRUE, .combine='c')%dopar% { 5*skyrms[i]*sqrt(length(which(sfa[[i]]>0))) }
    if (Magnitudes) {
      detecthres.mag<--2.5*(log10(detecthres)-log10(ABvegaflux))+magZP
    } else {
      detecthres.mag<-array(NA, dim=c(length(dfaflux)))
    }
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    #}}}
    #If wanted, do sky subtraction {{{
    if (doskyest) {
      if (!quiet) { message("Perfoming Sky Subtraction"); cat("   Performing Sky Subtraction") }
      #Subrtract Sky Flux
      dfaflux<-dfaflux-skyflux
      sfaflux<-sfaflux-skyflux
      dfaerr<-sqrt(dfaerr^2+skyerr^2)
      dfaerr<-sqrt(sfaerr^2+skyerr^2)
      if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    }#}}}
  } else {
    skylocal<-array(NA, dim=c(length(dfaflux)))
    skyflux<-array(NA, dim=c(length(dfaflux)))
    skyerr<-array(NA, dim=c(length(dfaflux)))
    skyrms<-array(NA, dim=c(length(dfaflux)))
    skypval<-array(NA, dim=c(length(dfaflux)))
    detecthres<-array(NA, dim=c(length(dfaflux)))
    detecthres.mag<-array(NA, dim=c(length(dfaflux)))
  }#}}}

  #Apply Flux Correction to the Finalised Values {{{
  if (!quiet) { cat("   Performing Final Calculations   ") }
  if (fluxcorr!=1) {
    sfaflux<-sfaflux*fluxcorr
    sfafluxw<-sfafluxw*fluxcorr
    dfaflux<-dfaflux*fluxcorr
    sfaerr<-sfaerr*fluxcorr
    sfaerrw<-sfaerrw*fluxcorr
    dfaerr<-dfaerr*fluxcorr
  }#}}}

  #Calculate Magnitudes {{{
  if (Magnitudes) {
    mags<--2.5*(log10(dfaflux)-log10(ABvegaflux))+magZP
  } else {
    mags<-array(NA, dim=c(length(dfaflux)))
  } #}}}

  #-----Diagnostic-----# {{{
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
  }#}}}

  #Check that final values of Deblended Convolved Apertures are not NA/NaN/Inf {{{
  if (length(which(!is.finite(dfaerr[which(sdfa > 0)]))) > 0) {
    message(paste(length(!is.finite(dfaerr[which(sdfa > 0)])), "elements of dfaerr are not finite"))
    sink(type="message")
    stop("NaN or Infs Produced in calculations")
  }#}}}

  if (!quiet) { cat(" - Done\n} Galaxy Results Complete\n") }
  #}}}

  #PART FIVE: OUTPUT {{{

  #If map was input in Jy/bm we need to convert it back before output in SourceSubtraction {{{
  if (Jybm) { ba=beamarea } else { ba=1. }
  #}}}

  #If wanted, make the Residual Map {{{
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
  }#}}}

  #If Subtracting Contaminants, then remove them before output {{{
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
    skylocal<-skylocal[which(contams==0)]
    skyerr <-skyerr[which(contams==0)]
    skyrms <-skyrms[which(contams==0)]
    skypval<-skypval[which(contams==0)]
    detecthres<-detecthres[which(contams==0)]
    detecthres.mag<-detecthres.mag[which(contams==0)]
    sfaerr <-sfaerr[which(contams==0)]
    sdfa   <-sdfa[which(contams==0)]
    sdfa2  <-sdfa2[which(contams==0)]
    sdfad  <-sdfad[which(contams==0)]
    sdfa2e2<-sdfa2e2[which(contams==0)]
    dfaflux<-dfaflux[which(contams==0)]
    dfaerr <-dfaerr[which(contams==0)]
    pixflux<-pixflux[which(contams==0)]
    mags   <-mags[which(contams==0)]
  }#}}}

  #If wanted, output the Results Table {{{
  if (writetab) {
    if (!quiet) { cat(paste('Writing Results Table to ',tableoutname,'   ')) }
    #Output the results table
    timer=system.time(writesfatableout(environment(), paste(pathout,tableoutname, sep="")) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Write Results Table - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }#}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) { message(paste('Sum of PSF = ',beamarea)) }
  #}}}

  #If wanted, open the browser for user inspection {{{
  if (interact) {
     cat(paste("Launching Interactive Mode: To end, type 'c'\n"))
     sink(type='message')
     browser()
     sink(sinkfile, type='message')
  }#}}}
  #}}}

  #Return {{{
  if (!quiet) { cat('\n') }
  detach(env)
  return=list(SFAflux=sfaflux,SFAerror=sfaerr, SFAfluxWeight=sfafluxw, SFAerrorWeight=sfaerrw, DFAflux=dfaflux,DFAerror=dfaerr)
  #}}}
}
