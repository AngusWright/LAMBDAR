fluxmeasurements <-
function(env=NULL) {
#Procedure measures the fluxes present in a supplied image
#inside catalogued apertures.

  #PART ZERO:  PSF DETERMINATION & GENERAL INITIALISATION /*fold*/ {{{
  # Load Parameter Space /*fold*/ {{{
  if (!is.null(env)) {
    attach(env)
  }
  # /*fend*/ }}}

  #Set function Environments /*fold*/ {{{
  environment(make_sa_mask)<-environment()
  environment(make_sfa_mask)<-environment()
  environment(make_a_mask)<-environment()
  environment(make_data_mask)<-environment()
  environment(readpsf)<-environment()
  environment(make_deblended_weightmap)<-environment()
  environment(writesfatableout)<-environment()
  environment(writedeblfractableout)<-environment()
  environment(writefitsout)<-environment()
  environment(sourcesubtraction)<-environment()
  # /*fend*/ }}}

  if (!quiet) { cat('Begining Flux Measurements\n') }
  message('} Initialisation Complete\nBegining Flux Measurements')
  dimim<-dim(image.env$im)

  #Get PSF details /*fold*/ {{{
  timer<-proc.time()
  if (!quiet) { cat('Getting PSF details') }
  message('Getting PSF details')
  if (!(nopsf)) {
    #There is a PSF specified /*fold*/ {{{
    #Details /*fold*/ {{{
    #Calculate PSF - if one has not be supplied,
    #then a gaussian PSF will be created. Stampsize of PSF
    #should be = maximum stampsize: max aperture * stamp mult /*fend*/ }}}
    #Get PSF /*fold*/ {{{
    psf<-readpsf(outenv=environment(),file.path(pathroot,pathwork,psfmap),asperpix,max(a_g),confidence,gauss_fwhm_as=gauss_fwhm_as)
    # /*fend*/ }}}
    #Notify /*fold*/ {{{
    if (verbose) { message(paste("Maxima of the PSF is at pixel", which(psf == max(psf)),"and has value",max(psf))) }
    # /*fend*/ }}}
    #Get radius of FWHM using FWHM confidence value = erf(2*sqrt(2*log(2))/sqrt(2)) /*fold*/ {{{
    psffwhm<-get.fwhm(psf)
    # /*fend*/ }}}
    #Normalise Beam Area /*fold*/ {{{
    beamarea_nn<-sumpsf
    beamarea_n<-as.single(sum(psf))
    # /*fend*/ }}}
    #-----Diagnostic-----## /*fold*/ {{{
    if (diagnostic) { message(paste('Beam area before/after norm: ',beamarea_nn,beamarea_n)) }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else {
    #There is no PSF specified /*fold*/ {{{
    #set relevant paramters manually
    beamarea_n<-1.0
    psfwidth<-0
    psf.clip<-0
    sumpsf<-0
    # /*fend*/ }}}
    #If no PSF, set FWHM to median aperture radius + buffer /*fold*/ {{{
    psffwhm<-round(min(a_g[which(a_g>0)])/asperpix)
    if ((is.na(psffwhm))|(!is.finite(psffwhm))) {
      #All Apertures are point sources, force width = 5 pixels /*fold*/ {{{
      #Details /*fold*/ {{{
      #There are no apertures that are not point sources,
      #and no PSF convolution - everything is a single pixel flux
      #set to a default of 5pixels /*fend*/ }}}
      psffwhm<-5
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
  }
  #If possible, Update Beamarea /*fold*/ {{{
  if (beamarea_pix == 0) {
    #Use the beamarea just determined /*fold*/ {{{
    beamarea<-beamarea_n
    # /*fend*/ }}}
  } else {
    # Otherwise just use the input beamarea /*fold*/ {{{
    beamarea<-beamarea_pix
    # /*fend*/ }}}
  }# /*fend*/ }}}
  # /*fend*/ }}}

  #-----Diagnostic-----# /*fold*/ {{{
  if (diagnostic) { message(paste('Beam area adopted: ',beamarea)) }
  # /*fend*/ }}}

  #Convert images from Jy/beam to Jy/pixel /*fold*/ {{{
  if (Jybm) {
    image.env$im<-image.env$im/beamarea
    if (length(image.env$ime)>1) { image.env$ime<-image.env$ime/beamarea }
    conf<-conf/beamarea
  } # /*fend*/ }}}

  #Notify /*fold*/ {{{
  if (showtime) { cat(paste(' - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )\n'))
    message(paste('Getting PSF details - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat(' - Done\n') }
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART ONE:  SINGLE APERTURE MASKS /*fold*/ {{{
  #Details /*fold*/ {{{
  #Create apertures for every object remaining in
  #the catalogue. /*fend*/ }}}
  #Get StampLengths for every galaxy /*fold*/ {{{
  if (psffilt) {
    stamplen<-(floor(ceiling(defbuff*a_g/asperpix)+(ceiling(psf.clip)/2))*2+5)
  } else {
    stamplen<-(floor(ceiling(defbuff*a_g/asperpix))*2+5)
  }
  # /*fend*/ }}}
  #Discard any contaminants that are beyond their nearest neighbours stamp /*fold*/ {{{
  #Number of Nearest Neighbors to search /*fold*/ {{{
  if (!exists("nNNs")) { nNNs<-10 }
  if (!exists("checkContam")) { checkContam<-TRUE }
  if (!exists("getLoners")) { getLoners<-FALSE }
  # /*fend*/ }}}
  if (getLoners) {
    timer<-proc.time()
    if (!quiet) { cat("Selecting out only objects that are completely alone ") }
    catlen<-length(x_g)
    nNNs<-max(nNNs,catlen-1)
    nearest<-nn2(data.frame(x_g,y_g),data.frame(x_g,y_g),k=nNNs+1)
    if (psffwhm!=0) {
      #If any of the nearest neighbors overlap with the object /*fold*/ {{{
      insidemask<-rowSums(nearest$nn.dist[,2:nNNs] < 4*psffwhm)==0
    } else {
      insidemask<-rep(NA,catlen)
      for(ind in 1:catlen) {
        #If any of the nearest neighbors overlap with the object /*fold*/ {{{
        insidemask[ind]<-any(nearest$nn.dist[ind,1:nNNs+1] < (sqrt(2)/2*(stamplen[nearest$nn.idx[ind,1:nNNs+1]]+stamplen[ind])))
        # /*fend*/ }}}
      }
      # /*fend*/ }}}
    }
    #Remove object catalogue entries /*fold*/ {{{
    x_g<-x_g[which(insidemask)]
    y_g<-y_g[which(insidemask)]
    id_g<-id_g[which(insidemask)]
    ra_g<-ra_g[which(insidemask)]
    dec_g<-dec_g[which(insidemask)]
    theta_g<-theta_g[which(insidemask)]
    a_g<-a_g[which(insidemask)]
    b_g<-b_g[which(insidemask)]
    if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
    if (exists("contams")) { contams<-contams[which(insidemask)] }
    chunkSize=length(id_g)/getDoParWorkers()
    mpiopts<-list(chunkSize=chunkSize)
    message("Number of objects per thread:",chunkSize)
    # /*fend*/ }}}
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(x_g),"supplied objects that are entirely unblended (",
                                  round(((catlen-length(x_g))/catlen)*100, digits=2),"% of objects were objects that were possibly blended)")) }
    # /*fend*/ }}}
    if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      message(paste('Remove irrelevant contaminants - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  } else if (checkContam && exists("contams") && length(which(contams==0))>1 && length(which(contams==1))>1) {
    timer<-proc.time()
    if (!quiet) { cat("Removing Contaminants that are irrelevant ") }
    catlen<-length(x_g)
    nNNs<-min(nNNs,length(which(contams==0))-1)
    nearest<-nn2(data.frame(x_g[which(contams==0)],y_g[which(contams==0)]),data.frame(x_g[which(contams==1)],y_g[which(contams==1)]),k=nNNs)
    contam.inside<-NULL
    for(ind in 1:length(which(contams==1))) {
      #If any of the nearest neighbors overlap with the contaminant /*fold*/ {{{
      contam.inside<-c(contam.inside,any(nearest$nn.dist[ind,] < (sqrt(2)/2*(stamplen[which(contams==0)][nearest$nn.idx[ind,]]+stamplen[which(contams==1)][ind]))))
    }
    insidemask<-rep(TRUE,catlen)
    insidemask[which(contams==1)]<-contam.inside
    # /*fend*/ }}}
    #Remove object catalogue entries /*fold*/ {{{
    x_g<-x_g[which(insidemask)]
    y_g<-y_g[which(insidemask)]
    id_g<-id_g[which(insidemask)]
    ra_g<-ra_g[which(insidemask)]
    dec_g<-dec_g[which(insidemask)]
    theta_g<-theta_g[which(insidemask)]
    a_g<-a_g[which(insidemask)]
    b_g<-b_g[which(insidemask)]
    if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
    contams<-contams[which(insidemask)]
    chunkSize=length(id_g)/getDoParWorkers()
    mpiopts<-list(chunkSize=chunkSize)
    message("Number of objects per thread:",chunkSize)
    # /*fend*/ }}}
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(x_g),"supplied objects & contaminants that intersect (",
                                  round(((catlen-length(x_g))/catlen)*100, digits=2),"% of objects were contaminants that didn't intersect galaxies)")) }
    # /*fend*/ }}}
    if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      message(paste('Remove irrelevant contaminants - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
  # /*fend*/ }}}
  #Convert decimal pixel values into actual pixel values /*fold*/ {{{
  x_p<-floor(x_g)
  y_p<-floor(y_g)
  # /*fend*/ }}}

  #Create an array of stamps containing the data image subsection for all objects /*fold*/ {{{
  timer=system.time(im_mask<-make_data_mask(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Data Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Remove Image arrays as they are no longer needed /*fold*/ {{{
  if (cutup) {
    imenvlist<-ls(envir=image.env)
    imenvlist<-imenvlist[which(imenvlist!="im"&imenvlist!="hdr_str")]
    rm(list=imenvlist, envir=image.env)
  }
  # /*fend*/ }}}
  #Create an array of stamps containing the apertures for all objects /*fold*/ {{{
  timer=system.time(sa<-make_sa_mask(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make SA Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Discard any apertures that were zero'd in the make process /*fold*/ {{{
  totsa<-foreach(sam=sa, i=1:length(sa), .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) } ; sum(sam) }
  insidemask<-(totsa > 0)
  if ((nloops<=1)&(length(which(insidemask==TRUE))==0)) { sink(type="message") ; stop("No Apertures are inside the Mask after Aperture Creation.") }  # Nothing inside the mask
  else if (length(which(insidemask==TRUE))==0) {
      warning("No Apertures are inside the Mask after Aperture Creation.")
      #Notify & Close Logfile /*fold*/ {{{
      if (!is.null(env)) {
        on.exit(detach(env), add=TRUE)
      }
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      if (!quiet) {
        cat(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      }
      return(NULL)
      # /*fend*/ }}}
  }
  # /*fend*/ }}}

  #Remove object that failed aperture creation /*fold*/ {{{
  sa<-sa[which(insidemask)]
  if (length(im_mask )>1) { im_mask<-im_mask[which(insidemask)] }
  if (length(imm_mask)>1) { imm_mask<-imm_mask[which(insidemask)] }
  if (length(ime_mask)>1) { ime_mask<-ime_mask[which(insidemask)] }
  stamplen<-stamplen[which(insidemask)]
  stamp_lims<-rbind(stamp_lims[which(insidemask),])
  mstamp_lims<-rbind(mstamp_lims[which(insidemask),])
  estamp_lims<-rbind(estamp_lims[which(insidemask),])
  im_stamp_lims<-rbind(im_stamp_lims[which(insidemask),])
  imm_stamp_lims<-rbind(imm_stamp_lims[which(insidemask),])
  ime_stamp_lims<-rbind(ime_stamp_lims[which(insidemask),])
  image_lims<-rbind(image_lims[which(insidemask),])
  mask_lims<-rbind(mask_lims[which(insidemask),])
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
  if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
  if (exists("contams")) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]
  chunkSize=ceiling(length(id_g)/getDoParWorkers())
  mpiopts<-list(chunkSize=chunkSize)
  message("Number of objects per thread:",chunkSize)
  # /*fend*/ }}}
  #Re-Initialise object count /*fold*/ {{{
  npos<-length(id_g)
  # /*fend*/ }}}
  #Create a full mask of all apertures in their correct image-space locations /*fold*/ {{{
  timer=system.time(image.env$aa<-make_a_mask(outenv=environment(), sa, dimim))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
   message(paste('Make A Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #If wanted, output the All Apertures Mask /*fold*/ {{{
  if (makeaamask) {
    if (!quiet) { cat(paste('Outputting All Apertures Mask to',aafilename,"   ")) }
    #Write All Apertures Mask to file /*fold*/ {{{
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,aafilename),image.env$aa,image.env$hdr_str,nochange=TRUE))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }# /*fend*/ }}}
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART TWO:  SINGLE FILTERED APERTURE MASKS /*fold*/ {{{
  #Details /*fold*/ {{{
  #If wanted, perform the convolution of the
  #apertures with the PSF. If not wanted,
  #duplicate the single apertures.
  #Also, perform flux weighting. /*fend*/ }}}
  #If wanted, use pixel fluxes as fluxweights /*fold*/ {{{
  if (usePixelFluxWeights) {
    #Image Flux at central pixel /*fold*/ {{{
    cat("Determine Image Flux at central pixel ")
    if (cutup) {
      pixflux<-foreach(xp=x_p-(im_stamp_lims[,1]-1),yp=y_p-(im_stamp_lims[,3]-1), im=im_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
            im[xp,yp]
      }
    } else {
      pixflux<-foreach(xp=x_p,yp=y_p, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment()), .export=c('image.env')) %dopar% {
            image.env$im[xp,yp]
      }
    }
    if (verbose) { cat(" - Done\n") }
    # /*fend*/ }}}
    #Determine Noise Characteristics /*fold*/ {{{
    if (verbose) { cat("Determine Image Noise Characteristics ") }
    if (any(pixflux==-99)) { warning("Some pixel flux determinations failed"); pixflux[which(pixflux==-99)]<-NA }
    quan<-quantile(image.env$im,c(0.1,0.9),na.rm=T)
    mode<-density(as.numeric(image.env$im),from=quan[1],to=quan[2],kernel='rect',bw=abs(quan[2]-quan[1])/100/sqrt(12),na.rm=T)
    mode<-mode$x[which.max(mode$y)]
    mad<-median(abs(image.env$im-mode),na.rm=T)
    if (length(which(pixflux>mode+mad))==0) {
      message("WARNING: No objects have pixel flux measurements that are > Pixel Mode+MAD. Pixel Flux weighting cannot be used\n")
        fluxweight<-1
    } else {
      fluxweight<-magmap(pixflux, lo=mode+mad, hi=max(pixflux), range=c(0.01,1), type="num", stretch='lin',bad=0.01)$map
    }
    cat(" - Done\n")
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #If needed, do Convolutions & Fluxweightings /*fold*/ {{{
  if ((!psffilt)&(length(which(fluxweight!=1))!=0)) {
    #No Convolution, Need Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If not convolving with psf, and if all the fluxweights are not unity,
    #then skip filtering, duplicate the arrays, and only weight the stamps/apertures
    # /*fend*/ }}}
    #No Convolution, duplicate apertures /*fold*/ {{{
    if (verbose) { message("No Convolution: Convolved Apertures are identical to Simple Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    # /*fend*/ }}}
    #Perform aperture weighting /*fold*/ {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(outenv=environment(), sa,fluxweightin=fluxweight,immask=imm_mask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if ((!psffilt)&(length(which(fluxweight!=1))==0)) {
    #No Convolution, No Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If not convolving with psf, and if all the fluxweights are unity,
    #then skip filtering and weighting, and simply duplicate the arrays /*fend*/ }}}
    #Duplicate arrays /*fold*/ {{{
    if (verbose) { message("No Convolution: Convolved Apertures are identical to Simple Apertures")
                   message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    wsfa<-sfa
    image.env$wfa<-image.env$fa
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if ((psffilt)&(length(which(fluxweight!=1))!=0)) {
    #Convolving PSF, Need Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If convolving with psf, and if all the fluxweights are not unity,
    #then colvolve and weight the stamps/apertures /*fend*/ }}}
    #Perform Convolution /*fold*/ {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(outenv=environment(), sa,immask=imm_mask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all convolved stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$fa<-make_a_mask(outenv=environment(), sfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Perform aperture weighting /*fold*/ {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(outenv=environment(), sa,fluxweightin=fluxweight,immask=imm_mask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all the convolved weighted stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if ((psffilt)&(length(which(fluxweight!=1))==0)) {
    #Convolve PSF, No Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If convolving with psf, and if all the fluxweights are unity,
    #colvolve, and then duplicate the stamps/apertures /*fend*/ }}}
    #Perform Convolution /*fold*/ {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(outenv=environment(), sa,immask=imm_mask))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all convolved stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$fa<-make_a_mask(outenv=environment(), sfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #No Fluxweights, duplicate apertures /*fold*/ {{{
    if (verbose) { message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    wsfa<-sfa
    image.env$wfa<-image.env$fa
    # /*fend*/ }}}
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #Discard any apertures that were zero'd in the make process /*fold*/ {{{
  totsfa<-foreach(sam=sfa, i=1:length(sfa), .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) } ; sum(sam) }
  insidemask<-(totsfa > 0)
  if ((nloops<=1)&(length(which(insidemask==TRUE))==0)) { sink(type="message") ; stop("No Apertures Remain within the mask after convolution.") }  # Nothing inside the mask
  else if (length(which(insidemask==TRUE))==0) {
      warning("No Apertures Remain within the mask after convolution.")
      #Notify & Close Logfile /*fold*/ {{{
      if (!is.null(env)) {
        on.exit(detach(env), add=TRUE)
      }
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      if (!quiet) {
        cat(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      }
      return(NULL)
      # /*fend*/ }}}
    }
  #Remove object that failed aperture creation /*fold*/ {{{
  sa<-sa[which(insidemask)]
  sfa<-sfa[which(insidemask)]
  wsfa<-wsfa[which(insidemask)]
  if (length(im_mask )>1) { im_mask<-im_mask[which(insidemask)] }
  if (length(imm_mask)>1) { imm_mask<-imm_mask[which(insidemask)] }
  if (length(ime_mask)>1) { ime_mask<-ime_mask[which(insidemask)] }
  stamplen<-stamplen[which(insidemask)]
  stamp_lims<-rbind(stamp_lims[which(insidemask),])
  mstamp_lims<-rbind(mstamp_lims[which(insidemask),])
  estamp_lims<-rbind(estamp_lims[which(insidemask),])
  im_stamp_lims<-rbind(im_stamp_lims[which(insidemask),])
  imm_stamp_lims<-rbind(imm_stamp_lims[which(insidemask),])
  ime_stamp_lims<-rbind(ime_stamp_lims[which(insidemask),])
  image_lims<-rbind(image_lims[which(insidemask),])
  mask_lims<-rbind(mask_lims[which(insidemask),])
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
  if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
  if (exists("contams")) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]
  chunkSize=ceiling(length(id_g)/getDoParWorkers())
  mpiopts<-list(chunkSize=chunkSize)
  message("Number of objects per thread:",chunkSize)
  # /*fend*/ }}}
  # /*fend*/ }}}
  #Re-Initialise object count /*fold*/ {{{
  npos<-length(id_g)
  # /*fend*/ }}}
  #Remove Uneeded Image /*fold*/ {{{
  rm(aa, envir=image.env)
  gc()
  # /*fend*/ }}}
  #-----Diagnostic-----# /*fold*/ {{{
  if (diagnostic) {
    if (verbose) { message("Checking Apertures for scaling errors") }
    foreach(sam=sa, sfam=sfa, i=1:length(sa), .combine=function(a,b){NULL},.inorder=FALSE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) }
      if (max(sfam)>1){warning(paste("Max of Convolved Aperture",i,"is",max(sfam))) }
    }
  }# /*fend*/ }}}
  #Do we want to plot a sample of the apertures? /*fold*/ {{{
  if (plotsample) {
    #Set output name /*fold*/ {{{
    pdf(file.path(pathroot,pathwork,pathout,"PSF_Samples.pdf"))
    # /*fend*/ }}}
    #Set Layout /*fold*/ {{{
    par(mfrow=c(2,2))
    # /*fend*/ }}}
    #Output a random 15 apertures to file /*fold*/ {{{
    ind1=which(a_g>0)
    ind1=ind1[order(runif(1:length(ind1)))][1:15]
    ind2=order(runif(1:npos))[1:15]
    ind<-c(ind1, ind2)
    rm(ind1)
    rm(ind2)
    ind<-ind[which(!is.na(ind))]
    for (i in ind) {
      image(sa[[i]], main="Single Aperture (SA)", asp=1, col=heat.colors(256), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
      image(sfa[[i]], main="Single Convolved Apertrure (SFA)", asp=1, col=heat.colors(256), useRaster=TRUE)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
    }# /*fend*/ }}}
    #Close the file /*fold*/ {{{
    dev.off()
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #Remove arrays that are no longer needed /*fold*/ {{{
  #rm(sa)
  #gc()
  # /*fend*/ }}}
  #If wanted, output the Convolved & Weighted Aperture Mask /*fold*/ {{{
  if (makefamask) {
    if (!quiet) { cat(paste('Outputting All Convolved Apertures Mask to',fafilename,"   ")) }
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,fafilename),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }# /*fend*/ }}}
  #If wanted, create (& possibly output) the SourceMask /*fold*/ {{{
  if (sourcemask) {
    #Details /*fold*/ {{{
    #> In the case of flux measurements: we want SM to be 1 only where the sources are, and within the mask region
    #> In the case of creating a source mask: we want SM to be 1 within the mask region in between sources only
    #To do this we:
    #         a) set the sm array = image mask (0s outside image, 1s inside)
    #         b) set all nonzero locations in fa to 0 in sm array
    # /*fend*/ }}}
    #Make Sourcemask /*fold*/ {{{
    if (!exists("TransmissionMap")) { TransmissionMap<-FALSE }
    if (!quiet) { cat(paste("Creating Sourcemask    ")) }
    if (length(image.env$imm)>1) { sm<-image.env$imm } else { sm<-array(1, dim=dimim) }
    if (!psffilt) {
      if (TransmissionMap) {
        sm[which(!image.env$fa > 0)]<-0
      } else {
        sm[which(image.env$fa > 0)]<-0
      }
    } else {
      if (TransmissionMap) {
        sm[which(!image.env$fa > max(psf, na.rm=TRUE)*(1-smConfidenceLim))]<-0
      } else {
        sm[which(image.env$fa > max(psf, na.rm=TRUE)*(1-smConfidenceLim))]<-0
      }
    }
    if (doskyest||getskyrms||BlankCor) {
      if (cutup) {
        imm_mask<-list(NULL)
        for (i in 1:npos) {
          imm_mask[[i]]<-sm[imm_stamp_lims[i,1]:imm_stamp_lims[i,2],imm_stamp_lims[i,3]:imm_stamp_lims[i,4]]
        }
      } else {
        image.env$imm<-sm
        imm_mask<-sm
      }
    }
    if (!quiet) { cat(" - Done\n") }
    # /*fend*/ }}}
    #-----Diagnoistic-----# /*fold*/ {{{
    if (diagnostic) {
      message(paste("SourceMask Max/Min:",max(sm),min(sm)))
      message(paste("OLDMethod - SourceMask Max/Min:",max(1-image.env$fa),min(1-image.env$fa)))
    }# /*fend*/ }}}
    #If wanted, output the SourceMask /*fold*/ {{{
    if (!exists("sourcemaskout")) { sourcemaskout<-FALSE }
    if (sourcemaskout){
      if (!quiet) { cat(paste('Outputting Source Mask to',smfilename,"   ")) }
      writefitsout(file.path(pathroot,pathwork,pathout,smfilename),sm,image.env$hdr_str,nochange=TRUE)
      if (!quiet) { cat(" - Done\n") }
    }# /*fend*/ }}}
    #If we want the SourceMask only, then end here /*fold*/ {{{
    if (sourcemaskonly) {
      if (!quiet) { cat("SourceMaskOnly Flag Set\n")  }
      return()
    } # /*fend*/ }}}
  } else if (plotsample) {
    #Set the mask to 1 everywhere /*fold*/ {{{
    sm<-1
    imm_mask<-1
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #Remove arrays that are no longer needed /*fold*/ {{{
  rm(fa, envir=image.env)
  gc()
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART THREE:  DEBLENDING /*fold*/ {{{
  #If wanted, perform sky estimation. Otherwise set to NA /*fold*/ {{{
  if (doskyest||getskyrms) {
    #Get sky estimates /*fold*/ {{{
    if (!quiet) { message("Perfoming Sky Estimation"); cat("Performing Sky Estimation") }
    #Perform Sky Estimation /*fold*/ {{{
    #timer<-system.time(skyest<-skyback(ra_g,dec_g,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,origim=list(dat=list(im)),maskim=list(dat=list(sm)),
    #                astrom=astr_struc,clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm))
    if (cutup) {
      timer<-system.time(skyest<-skyback.par(x_p-(im_stamp_lims[,1]-1),y_p-(im_stamp_lims[,3]-1),cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,im_mask=im_mask,imm_mask=imm_mask,
                      clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm, mpiopts=mpiopts))
    } else {
      timer<-system.time(skyest<-skyback.par(x_p,y_p,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,im_mask=image.env$im,imm_mask=image.env$imm,
                      clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm, mpiopts=mpiopts))
    }
    # /*fend*/ }}}
    #Get sky parameters /*fold*/ {{{
    skylocal<-skyest[,'sky']
    skylocal.mean<-skyest[,'sky.mean']
    if (correl.noise==0) {
      warning("Noise Correlation coefficient is set to 0; All flux errors will end up as 0")
      message("WARNING: Noise Correlation coefficient is set to 0; All flux errors will end up as 0")
    }
    skyerr<-skyest[,'skyerr']*correl.noise
    skyrms<-skyest[,'skyRMS']
    skypval<-skyest[,'skyRMSpval']
    skyNBinNear<-skyest[,'Nnearsky']
    skyerr.mean<-skyest[,'skyerr.mean']*correl.noise
    skyrms.mean<-skyest[,'skyRMS.mean']
    skypval.mean<-skyest[,'skyRMSpval.mean']
    skyNBinNear.mean<-skyest[,'Nnearsky.mean']
    if(!is.na(as.numeric(skydefault))) {
      skydefault<-as.numeric(skydefault)
    } else if (grepl("median",skydefault)) {
      skydefault<-median(skylocal,na.rm=TRUE)
    } else if (grepl("mean",skydefault)) {
      skydefault<-mean(skylocal, na.rm=TRUE)
    }
    if (is.na(skydefault)) {
      warning("All sky estimates have failed, and so the requested sky default is NA.\nTo stop bad behaviour, sky values are all being set to 0")
      skydefault<-0
    }
    skylocal[which(is.na(skylocal))]<-skydefault
    #skyflux<-skylocal*sdfa
    #skyerr<-skyerr*sdfa
    # /*fend*/ }}}
    #Notify /*fold*/ {{{
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Sky Estimate - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #If wanted, plot some of the Sky Estimates /*fold*/ {{{
    if (plotsample) {
      if (!quiet) { message("Plotting Sky Estimation"); cat("Plotting Sky Estimation") }
        if (cutup) {
          timer<-system.time(PlotSkyback(id_g,x_g,y_g,im_stamp_lims,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,im_mask=im_mask,imm_mask=imm_mask,
                          clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm,plotall=plotall,path=file.path(pathroot,pathwork,pathout),mpiopts=mpiopts))
        } else {
          timer<-system.time(PlotSkyback(id_g,x_g,y_g,im_stamp_lims,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,im_mask=image.env$im,imm_mask=image.env$imm,
                          clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm,plotall=plotall,path=file.path(pathroot,pathwork,pathout),mpiopts=mpiopts))
        }
      #Notify /*fold*/ {{{
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Plot Sky Estimate - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Calculate Detection Thresholds /*fold*/ {{{
    if (!quiet) { message("Calculating Detection Limits"); cat("Calculating Detection Limits") }
    detecthres<-foreach(sfam=sfa, srms=skyrms, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { 5*srms*sqrt(length(which(sfam>0))) }
    if (Magnitudes) {
      suppressWarnings(detecthres.mag<--2.5*(log10(detecthres)-log10(ABvegaflux))+magZP)
    } else {
      detecthres.mag<-array(NA, dim=c(length(sfa)))
    }
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    # /*fend*/ }}}
    #If wanted, do sky subtraction /*fold*/ {{{
    #if (doskyest) {
    #  if (!quiet) { message("Perfoming Sky Subtraction"); cat("   Performing Sky Subtraction") }
    #  #Subrtract Sky Flux
    #  dfaflux<-dfaflux-skyflux
    #  sfaflux<-sfaflux-skyflux
    #  dfaerr[which(!is.na(skyerr))]<-sqrt(dfaerr^2+skyerr^2)[which(!is.na(skyerr))]
    #  sfaerr[which(!is.na(skyerr))]<-sqrt(sfaerr^2+skyerr^2)[which(!is.na(skyerr))]
    #  if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    #}# /*fend*/ }}}
  } else {
    skylocal<-array(NA, dim=c(length(sfa)))
    skyflux<-array(NA, dim=c(length(sfa)))
    skyerr<-array(NA, dim=c(length(sfa)))
    skyrms<-array(NA, dim=c(length(sfa)))
    skypval<-array(NA, dim=c(length(sfa)))
    detecthres<-array(NA, dim=c(length(sfa)))
    detecthres.mag<-array(NA, dim=c(length(sfa)))
    skyNBinNear<-array(NA, dim=c(length(sfa)))
    skyNBinNear.mean<-array(NA, dim=c(length(sfa)))
    skyflux.mean<-array(NA, dim=c(length(sfa)))
    skyerr.mean   <- array(NA, dim=c(length(sfa)))
    skylocal.mean <- array(NA, dim=c(length(sfa)))
    skypval.mean  <- array(NA, dim=c(length(sfa)))
    skyrms.mean   <- array(NA, dim=c(length(sfa)))
  }# /*fend*/ }}}
  #Perform Deblending of apertures /*fold*/ {{{
  timer=system.time(dbw<-make_deblended_weightmap(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Deblended Weightmap - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Remove arrays that are no longer needed /*fold*/ {{{
  rm(wfa, envir=image.env)
  # /*fend*/ }}}
  #Finalise SFA Apertures /*fold*/ {{{
  sfabak<-sfa
  if (!exists("PSFWeighted")) { PSFWeighted<-FALSE }
  if (psffilt & !PSFWeighted) {
    # If we have convolved with a PSF, & not doing PSF Weighting, convert back to tophat apertures (now expanded by PSF convolution) /*fold*/ {{{
    # For each aperture, Binary filter at aplim*max(ap) /*fold*/ {{{
    sba<-foreach(sfam=sfa, .options.mpi=mpiopts, .noexport=ls(envir=environment()), .export="apLimit")%dopar%{
      apvals<-rev(sort(sfam))
      tempsum<-cumsum(apvals)
      tempfunc<-approxfun(tempsum,apvals)
      apLim<-tempfunc(apLimit*max(tempsum, na.rm=TRUE))
      #message(paste("ApLimit:",apLimit,"; PSFLimit:",psfLimit))
      sfam[which(sfam <  apLim)]<-0
      sfam[which(sfam >= apLim)]<-1
      return=sfam
    }
    # /*fend*/ }}}
    #Update the sfa /*fold*/ {{{
    sfa<-sba
    rm(sba)
    # /*fend*/ }}}
    #If wanted, output the Convolved & Weighted Aperture Mask /*fold*/ {{{
    if (makefamask) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), sfa, dimim))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      fafilename<-paste("final_",fafilename,sep="")
      if (!quiet) { cat(paste('Outputting Re-Boxcar-ed All Convolved Apertures Mask to',fafilename,"   ")) }
      timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,fafilename),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      rm(wfa, envir=image.env)
    }# /*fend*/ }}}
    #remove unneeded variables /*fold*/ {{{
    rm(psfvals)
    rm(tempsum)
    rm(tempfunc)
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if (psffilt & PSFWeighted) {
    #If we have convolved by the PSF & want PSF Weighting /*fold*/ {{{
    #If wanted, output the Convolved & Weighted Aperture Mask /*fold*/ {{{
    if (makefamask) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), sfa, dimim))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      fafilename<-paste("final_",fafilename,sep="")
      if (!quiet) { cat(paste('Outputting Final All Convolved Apertures Mask to',fafilename,"   ")) }
      timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,fafilename),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      rm(wfa, envir=image.env)
    }# /*fend*/ }}}
    #remove unneeded variables /*fold*/ {{{
    rm(psfvals)
    rm(tempsum)
    rm(tempfunc)
    # /*fend*/ }}}
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #If PSF Supplied & sample wanted, Plot PSF and Min Aperture Correction /*fold*/ {{{
  if (!(nopsf)&(plotsample)) {
    #PSF with Contours /*fold*/ {{{
    pdf(file.path(pathroot,pathwork,pathout,"PSF.pdf"))
    psfvals<-rev(sort(psf))
    tempsum<-cumsum(psfvals)
    tempfunc<-approxfun(tempsum,psfvals)
    psfLimit<-tempfunc(apLimit*max(tempsum, na.rm=TRUE))
    suppressWarnings(image(log10(psf),main="PSF & Binary Contour Levels", asp=1,col=heat.colors(256),useRaster=TRUE))
    contour(psf, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
    if (psffilt) { contour(psf, levels=c(psfLimit), labels=c(apLimit), col='green', add=TRUE) }
    dev.off()
    # /*fend*/ }}}
    #Plot Example of PS minimum aperture Correction; minimum source /*fold*/ {{{
    ind<-(which(a_g==0)[1])
    # /*fend*/ }}}
    #If no point sources, use the minimum aperture /*fold*/ {{{
    if (is.na(ind)) { ind<-which.min(a_g) }
    # /*fend*/ }}}
    #Make Sure PSF is centred on centre of stamp /*fold*/ {{{
    centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
    delta<-floor(stamplen[ind]/2)*c(-1,+1)
    lims<-rbind(centre[1]+delta,centre[2]+delta)
    dx<-lims[1,1]-1
    dy<-lims[2,1]-1
    dx<-ifelse(dx>0,dx+1,dx)
    dy<-ifelse(dy>0,dy+1,dy)
    #Reinterpolate the PSF at point source XcenYcen /*fold*/ {{{
    lenxpsf<-length(psf[,1])
    lenypsf<-length(psf[1,])
    lenx<-length(1:stamplen[ind])
    leny<-length(1:stamplen[ind])
    # /*fend*/ }}}
    #Make grid for psf at old pixel centres /*fold*/ {{{
    psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Make expanded grid of new pixel centres /*fold*/ {{{
    expanded<-expand.grid(seq(1,lenx),seq(1,leny))
    xnew<-expanded[,1]-x_g[ind]%%1
    ynew<-expanded[,2]-y_g[ind]%%1
    # /*fend*/ }}}
    #Interpolate /*fold*/ {{{
    ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
    # /*fend*/ }}}
    #Aperture Correction Plot /*fold*/ {{{
    pdf(file.path(pathroot,pathwork,pathout,"ApertureCorrection.pdf"))
    if (!PSFWeighted) {
      #Binary Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (smallest source)", asp=1,col=heat.colors(256),useRaster=TRUE))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      contour(ap, levels=c(psfLimit), labels=c(apLimit), col='green', add=TRUE)
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=TRUE))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    } else {
      #PSF Weighted Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (smallest source)", asp=1,col=heat.colors(256),useRaster=TRUE))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=TRUE))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
    #Plot Example of PS minimum aperture Correction; median source /*fold*/ {{{
    ind<-(which.min(abs(a_g-median(a_g)))[1])
    #Make Sure PSF is centred on centre of stamp /*fold*/ {{{
    centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
    delta<-floor(stamplen[ind]/2)*c(-1,+1)
    lims<-rbind(centre[1]+delta,centre[2]+delta)
    dx<-lims[1,1]-1
    dy<-lims[2,1]-1
    dx<-ifelse(dx>0,dx+1,dx)
    dy<-ifelse(dy>0,dy+1,dy)
    # /*fend*/ }}}
    #Reinterpolate the PSF at point source XcenYcen /*fold*/ {{{
    lenxpsf<-length(psf[,1])
    lenypsf<-length(psf[1,])
    lenx<-length(1:stamplen[ind])
    leny<-length(1:stamplen[ind])
    # /*fend*/ }}}
    #Make grid for psf at old pixel centres /*fold*/ {{{
    psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
    # /*fend*/ }}}
    #Make expanded grid of new pixel centres /*fold*/ {{{
    expanded<-expand.grid(seq(1,lenx),seq(1,leny))
    xnew<-expanded[,1]-x_g[ind]%%1
    ynew<-expanded[,2]-y_g[ind]%%1
    # /*fend*/ }}}
    #Interpolate /*fold*/ {{{
    ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Aperture Correction Plot /*fold*/ {{{
    if (!PSFWeighted) {
      #Binary Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (median source)", asp=1,col=heat.colors(256),useRaster=TRUE))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      if (psffilt) { contour(ap, levels=c(psfLimit), labels=c(apLimit), col='green', add=TRUE) }
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=TRUE))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    } else {
      #PSF Weighted Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap)-log10(sfa[[ind]]),main="Example: Minimum Aperture Correction (median source; shown as residual)", asp=1,col=heat.colors(256),useRaster=TRUE))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    }
    dev.off()
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #Generate Deblended Flux Arrays /*fold*/ {{{
  if (!quiet) { cat("Generating Deblended Flux Arrays") }
  timer<-proc.time()
  #Create the Deblended Flux Array /*fold*/ {{{
  dfa<-foreach(dbwm=dbw, sfam=sfa, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
  # /*fend*/ }}}
  #Notify /*fold*/ {{{
  if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Make DFA - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Check if we're just calculating the deblend fraction /*fold*/ {{{
  if (!exists("getDeblFrac")) { getDeblFrac<-FALSE }
  if (getDeblFrac) {
    if (!quiet) { cat("Getting Deblend Fraction Only flag is TRUE.\nCalculating Aperture Integrals & Outputting catalogue.\n") }
    #Integral of the aperture; ssa /*fold*/ {{{
    if (verbose) { cat("Integral of the aperture") }
    ssa<-foreach(sam=sa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the convolved aperture; ssfa /*fold*/ {{{
    if (verbose) { cat("Integral of the convolved aperture") }
    ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the deblended convolved aperture; sdfa /*fold*/ {{{
    if (verbose) { cat("Integral of the deblended convolved aperture") }
    sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    spsf<-rep(sumpsf, length(sdfa))
    #If wanted, output the Results Table /*fold*/ {{{
    if (writetab) {
      #Output the results table /*fold*/ {{{
      if ((nloops!=1)&&(length(param.env$tableoutname)!=nloops)) {
        if (!quiet) { cat(paste('Writing Deblend Fraction Table to ',file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,"_file",f,".csv")),'   ')) }
        timer=system.time(writedeblfractableout(filename=file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,"_file",f,".csv"))) )
      } else {
        if (!quiet) { cat(paste('Writing Deblend Fraction Table to ',file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,".csv")),'   ')) }
        timer=system.time(writedeblfractableout(filename=file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,".csv"))) )
      }
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Write Deblend Fraction Table - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }# /*fend*/ }}}
    if (!quiet) { cat("....Done\n") }
    return()
  }
  # /*fend*/ }}}
  # /*fend*/ }}}
  #If wanted; Iterate the fluxes to improve final measurements /*fold*/ {{{
  if (iterateFluxes) {
    #Notify /*fold*/ {{{
    message(paste("Iterating Flux Determination",nIterations,"times {\n"))
    if (!quiet) { cat("Iterating Flux Determination",nIterations,"times {\n") }
    # /*fend*/ }}}
    #For the number of desired iterations /*fold*/ {{{
    #Setup for iteration /*fold*/ {{{
    weightType='flux'
    quietbak<-quiet
    psffiltbak<-psffilt
    fluxiters<-matrix(as.numeric(NA),ncol=nIterations,nrow=length(id_g))
    erriters<-fluxiters
    sdfaiters<-fluxiters
    #Integral of the convolved aperture; ssfa /*fold*/ {{{
    if (verbose) { cat("  Integral of the convolved aperture") }
    ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the deblended convolved aperture; sdfa /*fold*/ {{{
    if (verbose) { cat("  Integral of the deblended convolved aperture") }
    sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    sdfad<-sdfa*NA
    sdfa2e2<-sdfa*NA
    attach(image.env)
    #/*fend*/ }}}
    for (iter in 1:nIterations) {
      #Determine objects to iterate over /*fold*/ {{{
      if (iter!=1) {
        xind<-which(sdfa!=ssfa & sdfa!=0)
        iterateLost<-sdfa==0
      } else {
        xind<-1:length(sdfa)
      }
      #Check we have something to do! /*fold*/ {{{
      if (length(xind)==0) {
        if (verbose) { cat(" - Breaking out; no objects overlap, so no point iterating!\n") }
        break
      }
      #/*fend*/ }}}
      #Update MPI options /*fold*/ {{{
      chunkSize=ceiling(length(xind)/getDoParWorkers())
      mpiopts<-list(chunkSize=chunkSize)
      #/*fend*/ }}}
      #/*fend*/ }}}
      #Calculate the flux per object /*fold*/ {{{
      #Notify /*fold*/ {{{
      message(paste("Calculating Flux (#",iter,")"))
      if (!quiet) { cat("  Calculating Flux (#",iter,")") }
      # /*fend*/ }}}
      timer<-proc.time()
      if (cutup) {
        sdfad[xind]<-foreach(dfam=dfa[xind],im=im_mask[xind], xlo=stamp_lims[xind,1],xup=stamp_lims[xind,2], ylo=stamp_lims[xind,3],yup=stamp_lims[xind,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
           sum(dfam*(im[xlo:xup,ylo:yup]))
        }
      } else {
        sdfad[xind]<-foreach(dfam=dfa[xind], xlo=stamp_lims[xind,1],xup=stamp_lims[xind,2], ylo=stamp_lims[xind,3],yup=stamp_lims[xind,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment()),.export='im') %dopar% {
           sum(dfam*(im[xlo:xup,ylo:yup]))
        }
      }
      #Notify /*fold*/ {{{
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #Calculate the error per object /*fold*/ {{{
      #Notify /*fold*/ {{{
      message(paste("Calculating Error (#",iter,")"))
      if (!quiet) { cat("  Calculating Error (#",iter,")") }
      timer<-proc.time()
      # /*fend*/ }}}
      if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
        sdfa2e2[xind]<-foreach(dfam=dfa[xind], xlo=estamp_lims[xind,1],xup=estamp_lims[xind,2],ylo=estamp_lims[xind,3],yup=estamp_lims[xind,4], ime=ime_mask[xind], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
            sum((dfam*ime[xlo:xup,ylo:yup])^2.)
        }
      } else if (length(ime_mask)==1){
        if (ime_mask==1) {
          sdfa2e2[xind]<-foreach(dfam=dfa[xind], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
            sum((dfam)^2)
          }
        } else if (ime_mask==0) {
          sdfa2e2<-rep(0,length(sdfad))
        } else {
          sdfa2e2[xind]<-foreach(dfam=dfa[xind], .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
              sum((dfam*ime_mask)^2.)
          }
        }
      } else {
        sdfa2e2[xind]<-foreach(dfam=dfa[xind], xlo=estamp_lims[xind,1],xup=estamp_lims[xind,2],ylo=estamp_lims[xind,3],yup=estamp_lims[xind,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
            sum((dfam*ime[xlo:xup,ylo:yup])^2.)
        }
      }
      #Notify /*fold*/ {{{
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #If wanted, remove sky estimates /*fold*/ {{{
      if (!quiet) { message(paste("Calculating Deblend Fraction (#",iter,")")); cat(paste("  Calculating Deblend Fraction (#",iter,")")) }
      timer<-proc.time()
      sdfa[xind]<-foreach(dfam=dfa[xind], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
      #Subtract Sky Flux /*fold*/ {{{
      if (doskyest) {
        sdfad[xind]<-sdfad[xind]-skylocal[xind]*sdfa[xind]
      }
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
        message(paste('Sky Estimate - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat(" - Done\n") }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #Save the values of the fluxes for output /*fold*/ {{{
      fluxiters[,iter]<-sdfad
      erriters[,iter]<-sqrt(sdfa2e2)
      sdfaiters[,iter]<-sdfa
      # /*fend*/ }}}
      #Re-calculate the weighted apertures /*fold*/ {{{
      #Notify /*fold*/ {{{
      message(paste("Calculating Weighting Apertures (#",iter,")"))
      if (!quiet) { cat("  Calculating Weighting Apertures (#",iter,")") }
      # /*fend*/ }}}
      quiet<-TRUE
      if (PSFWeighted) {
        #If using PSFWeighted, we may be able to skip an unnecessary convolution /*fold*/ {{{
        psffilt<-FALSE
        timer=system.time(wsfa[xind]<-make_sfa_mask(outenv=environment(), sfa,fluxweightin=sdfad,subs=xind))
        # /*fend*/ }}}
      } else {
        #If not using PSFWeighted, we may need to re-convolve /*fold*/ {{{
        psffilt<-psffiltbak
        timer=system.time(wsfa[xind]<-make_sfa_mask(outenv=environment(), sa,fluxweightin=sdfad,subs=xind))
        # /*fend*/ }}}
      }
      timer2=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dimim,subs=xind))
      psffilt<-psffiltbak
      quiet<-quietbak
      #Notify /*fold*/ {{{
      if (showtime) { cat(" - Done (",round(timer[3]+timer2[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #Re-calculate the deblending matricies /*fold*/ {{{
      #Notify /*fold*/ {{{
      message(paste("Calculating Deblend Matricies (#",iter,")"))
      if (!quiet) { cat("  Calculating Deblend Matricies (#",iter,")") }
      # /*fend*/ }}}
      quiet<-TRUE
      timer2=system.time(dbw[xind]<-make_deblended_weightmap(outenv=environment(),subs=xind))
      quiet<-quietbak
      timer<-proc.time()
      dfa[xind]<-foreach(dbwm=dbw[xind], sfam=sfa[xind], .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
      #Notify /*fold*/ {{{
      if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3]+timer2[3],digits=2),"sec )\n")
      } else if (!quiet) { cat("   - Done\n") }
      # /*fend*/ }}}
      # /*fend*/ }}}
    }
    #For objects completely removed, make final measured flux the one that is output /*fold*/ {{{
    fluxiters[which(fluxiters==0)]<-NA
    erriters[which(erriters==0)]<-NA
    sdfaiters[which(sdfaiters==0)]<-NA
    for (i in 1:length(fluxiters[,1])) {
      fluxiters[i,which(is.na(fluxiters[i,]))]<-rev(fluxiters[i,which(!is.na(fluxiters[i,]))])[1]
      erriters[i,which(is.na(erriters[i,]))]<-rev(erriters[i,which(!is.na(erriters[i,]))])[1]
      sdfaiters[i,which(is.na(sdfaiters[i,]))]<-rev(sdfaiters[i,which(!is.na(sdfaiters[i,]))])[1]
    }
    #Update MPI options /*fold*/ {{{
    chunkSize=ceiling(length(id_g)/getDoParWorkers())
    mpiopts<-list(chunkSize=chunkSize)
    #/*fend*/ }}}
    #/*fend*/ }}}
    detach(image.env)
    if (!quiet) { cat("} Done\n") }
  # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #If wanted, output the Deblended & Convolved Aperture Mask /*fold*/ {{{
  if (makedfamask) {
    if (!quiet) { cat(paste('Making All Deblended Convolved Apertures Mask - ')) }
    timer=system.time(image.env$adfa<-make_a_mask(outenv=environment(), dfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (!quiet) { cat(paste('Outputting All Deblended Convolved Apertures Mask to',dfafilename,"   ")) }
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,dfafilename),image.env$adfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }# /*fend*/ }}}
  #Remove array that is no longer needed /*fold*/ {{{
  if (!plotsample) { rm(dbw) }
  rm(adfa, envir=image.env)
  gc()
  # /*fend*/ }}}
  # /*fend*/ }}}
  #If wanted, Perform Randoms Correction /*fold*/ {{{
  if (!exists("RanCor")) { RanCor<-FALSE }
  if (RanCor=="execute") {
    .executeRanCor()
    RanCor<-TRUE
  }
  if (RanCor) {
    if (!quiet) { message("Perfoming Randoms Correction"); cat("Performing Randoms Correction") }
    #debug(rancor.par)
    if (cutup) {
      timer<-system.time(randoms<-rancor.par(im_mask=im_mask,imm_mask=imm_mask,ap_mask=sfa,stamplims=stamp_lims,masklims=mstamp_lims,numIters=nRandoms,mpiopts=mpiopts,remask=FALSE))
    } else {
      timer<-system.time(randoms<-rancor.par(im_mask=image.env$im,imm_mask=image.env$imm,ap_mask=sfa,stamplims=image_lims,masklims=mask_lims,numIters=nRandoms,mpiopts=mpiopts,remask=FALSE))
    }
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Randoms Correction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (plotsample) {
      if (!quiet) { message("Plotting Randoms Correction"); cat("Plotting Randoms Correction") }
      #debug(rancor.par)
      if (cutup) {
        timer<-system.time(plotRanCor(id_g,x_g,y_g,im_mask=im_mask,imm_mask=imm_mask,ap_mask=sfa,stamplims=stamp_lims,imstamplims=im_stamp_lims,masklims=mstamp_lims,numIters=nRandoms,remask=FALSE,path=file.path(pathroot,pathwork,pathout),plotall=plotall))
      } else {
        timer<-system.time(plotRanCor(id_g,x_g,y_g,im_mask=image.env$im,imm_mask=image.env$imm,ap_mask=sfa,stamplims=image_lims,imstamplims=im_stamp_lims,masklims=mask_lims,numIters=nRandoms,remask=FALSE,path=file.path(pathroot,pathwork,pathout),plotall=plotall))
      }
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Plotting Randoms Correction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
  }
  # /*fend*/ }}}
  #If wanted, Perform Blanks Correction /*fold*/ {{{
  if (!exists("BlankCor")) { BlankCor<-FALSE }
  if (BlankCor) {
    if (!quiet) { message("Perfoming Blanks Correction"); cat("Performing Blanks Correction") }
    #debug(rancor.par)
    if (cutup) {
      timer<-system.time(blanks<-rancor.par(im_mask=im_mask,imm_mask=imm_mask,ap_mask=sfa,stamplims=stamp_lims,masklims=mstamp_lims,numIters=nBlanks,mpiopts=mpiopts,remask=TRUE))
    } else {
      timer<-system.time(blanks<-rancor.par(im_mask=image.env$im,imm_mask=image.env$imm,ap_mask=sfa,stamplims=image_lims,masklims=mask_lims,numIters=nBlanks,mpiopts=mpiopts,remask=TRUE))
    }
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Blanks Correction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (plotsample) {
      if (!quiet) { message("Plotting Blanks Correction"); cat("Plotting Blanks Correction") }
      #debug(rancor.par)
    if (cutup) {
      timer<-system.time(plotRanCor(id_g,x_g,y_g,im_mask=im_mask,imm_mask=imm_mask,ap_mask=sfa,stamplims=stamp_lims,imstamplims=im_stamp_lims,masklims=mstamp_lims,numIters=nBlanks,remask=TRUE,path=file.path(pathroot,pathwork,pathout),plotall=plotall))
    } else {
      timer<-system.time(plotRanCor(id_g,x_g,y_g,im_mask=image.env$im,imm_mask=image.env$imm,ap_mask=sfa,stamplims=image_lims,imstamplims=im_stamp_lims,masklims=mask_lims,numIters=nBlanks,remask=TRUE,path=file.path(pathroot,pathwork,pathout),plotall=plotall))
    }
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Plotting Blanks Correction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
  }
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART FOUR: COMPUTE AND OUTPUT GALAXY-BY-GALAXY RESULTS /*fold*/ {{{
  if (!quiet) { cat("Computing Galaxy-by-Galaxy Results {\n") }
  # Make Images Available /*fold*/ {{{
  if (!cutup) { attach(image.env) }
  # /*fend*/ }}}
  #-----Diagnostic-----# /*fold*/ {{{
  if (diagnostic){
    if (length(which(is.na(image.env$im)))!=0) {        message(paste("Input Parameter 'im' contains NA elements")) }
    if (length(which(is.na(ime_mask)))!=0) {        message(paste("Input Parameter 'ime' contains NA elements")) }
    if (length(which(is.na(sfa)))!=0) {        message(paste("Input Parameter 'sfa' contains NA elements")) }
    if (length(which(is.na(dfa)))!=0) {        message(paste("Input Parameter 'dfa' contains NA elements")) }
    if (length(which(is.na(fluxcorr)))!=0) {        message(paste("Input Parameter 'fluxcorr' contains NA elements")) }
  }# /*fend*/ }}}
  #Perform Calculations /*fold*/ {{{
  if (!quiet) { cat("   Calculating Fluxes ") }
  if (verbose) { cat("{\n") }
  timer<-proc.time()
#-----
  #Check for Saturations; saturated /*fold*/ {{{
  if (is.finite(saturation)) {
    if (verbose) { cat("      Checking for Saturations ") }
    if (cutup) {
      saturated<-foreach(sfam=sfa, im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment()), .export='saturation') %dopar% {
          any(im[xlo:xup,ylo:yup][which(sfam!=0,arr.ind=T)]>=saturation)
      }
    } else {
      saturated<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export=c('im','saturation'), xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          any(im[xlo:xup,ylo:yup][which(sfam!=0,arr.ind=T)]>=saturation)
      }
    }
    if (verbose) { cat(" - Done\n") }
  } else {
    saturated<-rep(FALSE,length(id_g))
  }
  # /*fend*/ }}}
#-----
  #Image Flux at central pixel; pixflux /*fold*/ {{{
  if (!usePixelFluxWeights) {
    if (verbose) { cat("      Image Flux at central pixel ") }
    if (cutup) {
      pixflux<-foreach(xp=x_p-(im_stamp_lims[,1]-1),yp=y_p-(im_stamp_lims[,3]-1), im=im_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
            im[xp,yp]
      }
    } else {
      pixflux<-foreach(xp=x_p,yp=y_p, .noexport=ls(envir=environment()), .export=c('im'), .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
            im[xp,yp]
      }
    }
    if (verbose) { cat(" - Done\n") }
  }
  # /*fend*/ }}}
#-----
  #Integral of the aperture; ssa /*fold*/ {{{
  if (verbose) { cat("      Integral of the aperture") }
  ssa<-foreach(sam=sa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sam) }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the convolved aperture; ssfa /*fold*/ {{{
  if (verbose) { cat("      Integral of the convolved aperture") }
  ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(unmodified convolved aperture)*(modified convolved aperture)]; ssfau /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(unmodified convolved aperture)*(convolved aperture)]") }
  ssfau<-foreach(sfam=sfa, usfam=sfabak, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sfam*usfam)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(convolved aperture)^2]; ssfa2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(convolved aperture)^2]") }
  ssfa2<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum((sfam)^2.)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * image); ssfad /*fold*/ {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image)") }
  if (cutup) {
    ssfad<-foreach(sfam=sfa, im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup]))
    }
  } else {
    ssfad<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export='im', xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup]))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Quartered Photometry - Quartered Integral of the (convolved aperture * image); qssfad /*fold*/ {{{
  if (verbose) { cat("      Quartered Integral of the (convolved aperture * image)") }
  if (cutup) {
    qssfad<-foreach(sfam=sfa, im=im_mask, x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                          y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
        cbind(sum(sfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(sfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(sfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(sfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  } else {
    qssfad<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export='im', x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                          y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind") %dopar% {
        cbind(sum(sfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(sfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(sfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(sfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the reinterpolated psf; spsf /*fold*/ {{{
  if (!nopsf) {
    if (verbose) { cat("      Integral of the reinterpolated psf") }
    spsf<-rep(sum(psf),length(id_g))
    if (verbose) { cat(" - Done\n") }
  } else {
    spsf<-rep(NA, length(sfa))
  }
  # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * psf); ssfap /*fold*/ {{{
  if (!nopsf) {
    if (verbose) { cat("      Integral of the (convolved aperture * psf)") }
    ssfap<-foreach(sfam=sfa, slen=stamplen, xc=x_g, yc=y_g, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment()), .export="psf") %dopar% {
      #for (i in 1:npos) {
      #slen=stamplen[i]
      #Make Sure PSF is centred on centre of stamp /*fold*/ {{{
      centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
      delta<-floor(slen/2)*c(-1,+1)
      lims<-rbind(centre[1]+delta,centre[2]+delta)
      dx<-lims[1,1]-1
      dy<-lims[2,1]-1
      dx<-ifelse(dx>0,dx+1,dx)
      dy<-ifelse(dy>0,dy+1,dy)
      # /*fend*/ }}}
      #Reinterpolate the PSF at point source XcenYcen /*fold*/ {{{
      lenxpsf<-length(psf[,1])
      lenypsf<-length(psf[1,])
      lenx<-length(1:slen)
      leny<-length(1:slen)
      #Make grid for psf at old pixel centres /*fold*/ {{{
      psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
      # /*fend*/ }}}
      #Make expanded grid of new pixel centres /*fold*/ {{{
      expanded<-expand.grid(seq(1,lenx),seq(1,leny))
      xnew<-expanded[,1]-xc%%1
      ynew<-expanded[,2]-yc%%1
      # /*fend*/ }}}
      #Interpolate /*fold*/ {{{
      ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
      # /*fend*/ }}}
      # /*fend*/ }}}
      sum(sfam*ap)
    }
    #}
    if (verbose) { cat(" - Done\n") }
  } else {
    ssfap<-rep(NA, length(sfa))
  }
  # /*fend*/ }}}
#-----
  #Integral of the deblended convolved aperture; sdfa /*fold*/ {{{
  if (verbose) { cat("      Integral of the deblended convolved aperture") }
  sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(deblended convolved aperture)^2]; sdfa2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture)^2]") }
  sdfa2<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum((dfam)^2.)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (deblended convolved aperture * image); sdfad /*fold*/ {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * image)") }
  if (cutup) {
    sdfad<-foreach(dfam=dfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(dfam*(im[xlo:xup,ylo:yup]))
    }
  } else {
    sdfad<-foreach(dfam=dfa,.noexport=ls(envir=environment()), .export='im', xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(dfam*(im[xlo:xup,ylo:yup]))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Quartered Photometry - Quartered Integral of the (deblended convolved aperture * image); qsdfad /*fold*/ {{{
  if (verbose) { cat("      Quartered Integral of the (deblended convolved aperture * image)") }
  if (cutup) {
    qsdfad<-foreach(dfam=dfa, im=im_mask, x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                          y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
        cbind(sum(dfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(dfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(dfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(dfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  } else {
    qsdfad<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export='im', x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                          y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind") %dopar% {
        cbind(sum(dfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(dfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(dfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(dfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (deblended convolved aperture * convolved aperture); sdfasfa /*fold*/ {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * convolved aperture)") }
  sdfasfa<-foreach(dfam=dfa, sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam*sfam) }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * image error); ssfae /*fold*/ {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    ssfae<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      sum(sfam*ime[xlo:xup,ylo:yup])
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfae<-ssfa
    } else if (ime_mask==0) {
      ssfae<-rep(0,length(ssfa))
    } else {
      ssfae<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(sfam*ime_mask)
      }
    }
  } else {
    ssfae<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
      sum(sfam*ime[xlo:xup,ylo:yup])
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * (image error)^2); ssfae2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    ssfae2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*(ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfae2<-ssfa
    } else if (ime_mask==0) {
      ssfae2<-rep(0,length(ssfa))
    } else {
      ssfae2<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum(sfam*(ime_mask)^2.)
      }
    }
  } else {
    ssfae2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(sfam*(ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(convolved aperture * image error)^2]; ssfa2e2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image error)^2]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    ssfa2e2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((sfam*ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfa2e2<-ssfa2
    } else if (ime_mask==0) {
      ssfa2e2<-rep(0,length(ssfa))
    } else {
      ssfa2e2<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum((sfam*ime_mask)^2.)
      }
    }
  } else {
    ssfa2e2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum((sfam*ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (deblended convolved aperture * image error); sdfae /*fold*/ {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    sdfae<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      sum(dfam*ime[xlo:xup,ylo:yup])
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      sdfae<-sdfa
    } else if (ime_mask==0) {
      sdfae<-rep(0,length(sdfa))
    } else {
      sdfae<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(dfam*ime_mask)
      }
    }
  } else {
    sdfae<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
      sum(dfam*ime[xlo:xup,ylo:yup])
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [deblended convolved aperture * (image error)^2]; sdfae2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [deblended convolved aperture * (image error)^2]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    sdfae2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(dfam*(ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      sdfae2<-sdfa2
    } else if (ime_mask==0) {
      sdfae2<-rep(0,length(ssfa))
    } else {
      sdfae2<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum(dfam*(ime_mask)^2.)
      }
    }
  } else {
    sdfae2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(dfam*(ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(deblended convolved aperture * image error)^2]; sdfa2e2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture * image error)^2]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    sdfa2e2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      sdfa2e2<-sdfa2
    } else if (ime_mask==0) {
      sdfa2e2<-rep(0,length(ssfa))
    } else {
      sdfa2e2<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum((dfam*ime_mask)^2.)
      }
    }
  } else {
    sdfa2e2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(convolved aperture * image) / {(image error)^2} ]; ssfadw /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image) / {(image error)^2} ]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    if (cutup) {
      ssfadw<-foreach(sfam=sfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], xelo=estamp_lims[,1],xeup=estamp_lims[,2],yelo=estamp_lims[,3],yeup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum(sfam*(im[xlo:xup,ylo:yup])/((ime[xelo:xeup,yelo:yeup])^2.))
      }
    } else {
      ssfadw<-foreach(sfam=sfa,.noexport=ls(envir=environment()), .export=c('im','ime'), xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], xelo=estamp_lims[,1],xeup=estamp_lims[,2],yelo=estamp_lims[,3],yeup=estamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum(sfam*(im[xlo:xup,ylo:yup])/((ime[xelo:xeup,yelo:yeup])^2.))
      }
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfadw<-ssfad
    } else if (ime_mask==0) {
      ssfadw<-rep(Inf,length(ssfa))
    } else {
      if (cutup) {
        ssfadw<-foreach(sfam=sfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
            sum(sfam*(im[xlo:xup,ylo:yup])/((ime_mask)^2.))
        }
      } else {
        ssfadw<-foreach(sfam=sfa,.noexport=ls(envir=environment()), .export=c('im',"ime_mask"), xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
            sum(sfam*(im[xlo:xup,ylo:yup])/((ime_mask)^2.))
        }
      }
    }
  } else {
    ssfadw<-foreach(sfam=sfa,.noexport=ls(envir=environment()), .export=c('im','ime'), xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], xelo=estamp_lims[,1],xeup=estamp_lims[,2],yelo=estamp_lims[,3],yeup=estamp_lims[,4], .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup])/((ime[xelo:xeup,yelo:yeup])^2.))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [convolved aperture / (image error)^2 ]; ssfaw /*fold*/ {{{
  if (verbose) { cat("      Integral of the [convolved aperture / (image error)^2 ]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    ssfaw<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam/(ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfaw<-ssfa
    } else if (ime_mask==0) {
      ssfaw<-rep(Inf,length(ssfa))
    } else {
      ssfaw<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum(sfam/(ime_mask)^2.)
      }
    }
  } else {
    ssfaw<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum(sfam/(ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(convolved aperture / image error)^2]; ssfa2w /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(convolved aperture / image error)^2]") }
  if (length(ime_mask)>1|(length(ime_mask)==1 & is.list(ime_mask))) {
    ssfa2w<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((sfam/ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfa2w<-ssfa2
    } else if (ime_mask==0) {
      ssfa2w<-rep(Inf,length(ssfa))
    } else {
      ssfa2w<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="ime_mask", .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
          sum((sfam/ime_mask)^2.)
      }
    }
  } else {
    ssfa2w<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpiopts) %dopar% {
        sum((sfam/ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Final Flux & Error Calculations /*fold*/ {{{
  if (!cutup) { detach(image.env) }
  if (verbose) { cat("      Final Fluxes and Error Calculations ") }
  #Calculate Aperture-Weighting Correction /*fold*/ {{{
  WtCorr<-ssfa/ssfau
  # /*fend*/ }}}
  #Calculate Minimum Aperture Correction /*fold*/ {{{
  if (!nopsf) {
    ApCorr<-spsf/ssfap
  } else {
    ApCorr<-1
  }
  # /*fend*/ }}}

  #Convolved aperture flux = Int(fltAp*Im) /*fold*/ {{{
  sfaflux<-ssfad
  # /*fend*/ }}}

  #Deblended convolved aperture flux = Int(DBfltAp*Im) /*fold*/ {{{
  dfaflux<-sdfad
  # /*fend*/ }}}

  #Convolved aperture error /*fold*/ {{{
  sfaerr<-sqrt((ssfae2)*ApCorr^2 + ((conf*beamarea)^2.*sqrt(ssfa)))
  # /*fend*/ }}}

  #Deblend error /*fold*/ {{{
  deblerr<-((1-sdfa/ssfa)*(1/sqrt(12))*abs(sfaflux)*ApCorr)
  # /*fend*/ }}}

  #Deblended Convolved aperture error /*fold*/ {{{
  dfaerr<-sqrt((sdfa2e2)*ApCorr^2 + ((conf*beamarea)^2.*sqrt(sdfa)) + (deblerr)^2)
  # /*fend*/ }}}

  #Convolved Aperture Flux Weights & Error Weights /*fold*/ {{{
  sfafluxw<-ssfadw*ssfa/ssfa2w
  sfaerrw<-sqrt(1/ssfaw*(ssfa/beamarea*ssfa/ssfa2)^2. + ((conf*beamarea)^2.*sqrt(ssfa)))
  # /*fend*/ }}}
  if (verbose) { cat(" - Done\n   } ") }
  # /*fend*/ }}}
#-----
  #Notify /*fold*/ {{{
  if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Perform Calculations - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat(" - Done\n") }
  # /*fend*/ }}}
  #Final Calculations /*fold*/ {{{
  if (!quiet) { cat("   Performing Final Calculations {  ") }
  if (doskyest|getskyrms) {
    skyflux.mean<-skylocal.mean*sdfa
    skyflux<-skylocal*sdfa
    if (doskyest) {
      if (!quiet) { message("Perfoming Sky Subtraction"); cat("\n   Performing Sky Subtraction") }
      #Subrtract Sky Flux /*fold*/ {{{
      dfaflux<-dfaflux-skyflux
      sfaflux<-sfaflux-skyflux
      dfaerr[which(!is.na(skyerr))]<-sqrt(dfaerr^2+(skyerr*sdfa)^2)[which(!is.na(skyerr))]
      sfaerr[which(!is.na(skyerr))]<-sqrt(sfaerr^2+(skyerr*sdfa)^2)[which(!is.na(skyerr))]
      if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    }
  }
  # /*fend*/ }}}
  #Apply Aperture Correction to the Aperture Fluxess /*fold*/ {{{
  sfaflux<-sfaflux*ApCorr
  dfaflux<-dfaflux*ApCorr
  # /*fend*/ }}}
  #Apply Additional Flux Correction to the Finalised Values /*fold*/ {{{
  if (fluxcorr!=1) {
    sfaflux<-sfaflux*fluxcorr
    sfafluxw<-sfafluxw*fluxcorr
    dfaflux<-dfaflux*fluxcorr
    sfaerr<-sfaerr*fluxcorr
    sfaerrw<-sfaerrw*fluxcorr
    dfaerr<-dfaerr*fluxcorr
  }# /*fend*/ }}}
  #Calculate Magnitudes /*fold*/ {{{
  if (Magnitudes) {
    suppressWarnings(mags<--2.5*(log10(dfaflux)-log10(ABvegaflux))+magZP)
  } else {
    mags<-array(NA, dim=c(length(dfaflux)))
  } # /*fend*/ }}}
  #-----Diagnostic-----# /*fold*/ {{{
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
  }# /*fend*/ }}}
  #Check that final values of Deblended Convolved Apertures are not NA/NaN/Inf /*fold*/ {{{
  #if (length(which(!is.finite(dfaerr[which(sdfa > 0)]))) > 0) {
    #message(paste(length(!is.finite(dfaerr[which(sdfa > 0)])), "elements of dfaerr are not finite"))
    #sink(type="message")
    #stop("NaN or Infs Produced in calculations")
  #}# /*fend*/ }}}
  if (!quiet) { cat("} Galaxy Results Complete\n") }
  # /*fend*/ }}}
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART FIVE: OUTPUT /*fold*/ {{{
  #Do we want to plot a sample of the apertures? /*fold*/ {{{
  if (plotsample) {
    #Set output name /*fold*/ {{{
    dir.create(file.path(pathroot,pathwork,pathout,"COGs"),showWarnings=FALSE)
    # /*fend*/ }}}
    #Determine Sample to Plot /*fold*/ {{{
    if (!exists("plotall")) { plotall<-FALSE }
    if (!plotall) {
      if (!quiet) { message("Writing Sample of COGs to File"); cat("Writing Sample of COGs to File") }
      #Output a random 15 apertures to file /*fold*/ {{{
      ind1=which(a_g>0 & sdfa!=0)
      ind1=ind1[order(sdfad[ind1],decreasing=T)][1:15]
      ind2=which(sdfa!=0)[order(runif(1:length(which(sdfa!=0))))[1:15]]
      ind<-c(ind1, ind2)
      rm(ind1)
      rm(ind2)
      ind<-ind[which(!is.na(ind))]
      # /*fend*/ }}}
    } else {
      #All /*fold*/ {{{
      if (!quiet) { message("Writing All COGs to File"); cat("Writing All COGs to File") }
      ind=c(which(a_g>0),which(a_g==0))
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
    for (i in ind) {
      #Open Device /*fold*/ {{{
      #png(file.path(pathroot,pathwork,pathout,paste("COGs/",id_g[i],".png",sep="")),width=14*240,height=3.5*240,res=240)
      pdf(file.path(pathroot,pathwork,pathout,paste("COGs/",id_g[i],".pdf",sep="")),width=14,height=3.5)
      # /*fend*/ }}}
      #Set Layout /*fold*/ {{{
      layout(cbind(1,2,3,4))
      # /*fend*/ }}}
      #Axes Limits /*fold*/ {{{
      #xlims=round(max(abs(seq(1,(diff(range(stamp_lims[i,1]:stamp_lims[i,2]))+1))-(diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1)/2)*asperpix, na.rm=TRUE))*c(-1,1)
      #xlims=range((seq(1,(diff(range(stamp_lims[i,1]:stamp_lims[i,2]))+1))-(x_p[i]-im_stamp_lims[i,1]))*asperpix)
      #ylims=range((seq(1,(diff(range(stamp_lims[i,3]:stamp_lims[i,4]))+1))-(y_p[i]-im_stamp_lims[i,3]))*asperpix)
      xlims=stamplen[i]*c(-2/3,2/3)*asperpix
      ylims=stamplen[i]*c(-2/3,2/3)*asperpix
      # /*fend*/ }}}
      #Make Ap and Source Mask Block/Trans matricies /*fold*/ {{{
      apT<-sfa[[i]]
      apT[which(sfa[[i]]==0)]<-NA
      apT[which(sfa[[i]]!=0)]<-1
      apB<-sfa[[i]]
      apB[which(sfa[[i]]==0)]<-1
      apB[which(sfa[[i]]!=0)]<-NA
      if (sourcemask) {
        smB<-sm[imm_stamp_lims[i,1]:imm_stamp_lims[i,2],imm_stamp_lims[i,3]:imm_stamp_lims[i,4]]
        smB[which(smB==0)]<-NA
      } else {
        smB<-1
      }
      # /*fend*/ }}}
      #Plot Image in greyscale /*fold*/ {{{
      Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
      suppressWarnings(image(x=(seq(1,(diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1))-(x_p[i]-im_stamp_lims[i,1]))*asperpix,y=(seq(1,(diff(range(im_stamp_lims[i,3]:im_stamp_lims[i,4]))+1))-(y_p[i]-im_stamp_lims[i,3]))*asperpix, z=log10(image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast, xlab="", ylab="",xlim=xlims, ylim=ylims,axes=FALSE))
      # /*fend*/ }}}
      #Plot Aperture in Blue /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=(sfa[[i]]), main="Image & Aperture", asp=1, col=hsv(2/3,seq(0,1,length=256)), useRaster=Rast, axes=FALSE, xlab="", ylab="", add=TRUE))
      # /*fend*/ }}}
      #Plot Image in greyscale /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,(diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1))-(x_p[i]-im_stamp_lims[i,1]))*asperpix,y=(seq(1,(diff(range(im_stamp_lims[i,3]:im_stamp_lims[i,4]))+1))-(y_p[i]-im_stamp_lims[i,3]))*asperpix, z=log10(image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast,add=TRUE, xlab="", ylab=""))
      # /*fend*/ }}}
      #Overlay Sourcemask in Green /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,(diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1))-(x_p[i]-im_stamp_lims[i,1]))*asperpix,y=(seq(1,(diff(range(im_stamp_lims[i,3]:im_stamp_lims[i,4]))+1))-(y_p[i]-im_stamp_lims[i,3]))*asperpix, z=log10(smB*image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]]), main="", asp=1, useRaster=Rast,add=TRUE, xlab="", ylab="",col=cm.colors(256)))
      # /*fend*/ }}}
      #Plot +ve flux in aperture in Heat Colours /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(apT*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]), main="", asp=1, col=heat.colors(256), useRaster=Rast,add=TRUE, xlab="", ylab=""))
      # /*fend*/ }}}
      #Overlay Aperture in Black /*fold*/ {{{
      #suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(sfa[[i]]), main="Image & Aperture", asp=1, useRaster=FALSE, axes=FALSE, xlab="", ylab="", xlim=xlims, ylim=xlims,add=TRUE))
      # /*fend*/ }}}
      #Plot Sources /*fold*/ {{{
      points(x=(x_p-x_p[i]+1)*asperpix,y=(y_p-y_p[i]+1)*asperpix, pch=3)
      # /*fend*/ }}}
      #Label with ID /*fold*/ {{{
      label("topleft",lab=id_g[i],cex=1.5, col='red')
      # /*fend*/ }}}
      #Draw Axes /*fold*/ {{{
      magaxis(frame.plot=T,main="Image & Aperture",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      # /*fend*/ }}}
      #Generate COGs /*fold*/ {{{
      #Raw Image /*fold*/ {{{
      if (PSFWeighted) {
        cog<-get.cog(sfa[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        cog<-get.cog(image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]],centre=c((diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1)/2, (diff(range(im_stamp_lims[i,3]:im_stamp_lims[i,4]))+1)/2))
      }
      # /*fend*/ }}}
      #Sky Subtracted only /*fold*/ {{{
      if (PSFWeighted) {
        cog_nosky<-get.cog(sfa[[i]]*(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        cog_nosky<-get.cog(image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]]-skylocal[i],centre=c((diff(range(im_stamp_lims[i,1]:im_stamp_lims[i,2]))+1)/2, (diff(range(im_stamp_lims[i,3]:im_stamp_lims[i,4]))+1)/2))
      }
      # /*fend*/ }}}
      #Deblended only /*fold*/ {{{
      if (PSFWeighted) {
        debl.cog<-get.cog(sfa[[i]]*dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        debl.cog<-get.cog(dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      }
      # /*fend*/ }}}
      #Deblended and Sky Subtracted /*fold*/ {{{
      if (PSFWeighted) {
        debl.cog_nosky<-get.cog(sfa[[i]]*dbw[[i]]*(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        debl.cog_nosky<-get.cog(dbw[[i]]*(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2))
      }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #Plot COGs /*fold*/ {{{
      if (Magnitudes) {
        #Plot in Magnitude Space /*fold*/ {{{
        #Get y limits /*fold*/ {{{
        suppressWarnings(ylim<-c((-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP)+3, (-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP)-3))
        # /*fend*/ }}}
        #If a limit in ±Inf, plot around median /*fold*/ {{{
        if (!all(is.finite(ylim))) { ylim=median(-2.5*(log10(cog$y)-log10(ABvegaflux))+magZP, na.rm=TRUE)+c(-1,3) }
        if (!all(is.finite(ylim))) { warning("Cog flux is always < 0; Using arbitrary plot limits"); ylim=18+c(-1,3) }
        # /*fend*/ }}}
        #Plot Raw COG /*fold*/ {{{
        suppressWarnings(magplot(x=cog$x*asperpix, y=-2.5*(log10(cog$y)-log10(ABvegaflux))+magZP, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Magnitude",
                ylim=ylim,type='l',lty=2,main="Curve of Growth"))
        # /*fend*/ }}}
        #Add Lines showing fluxes /*fold*/ {{{
        #Undeblended /*fold*/ {{{
        suppressWarnings(abline(h=-2.5*(log10(sfaflux[i])-log10(ABvegaflux))+magZP, lwd=1, col='orange', lty=1))
        # /*fend*/ }}}
        #Deblended /*fold*/ {{{
        suppressWarnings(abline(h=-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP, lwd=1, col='green', lty=1))
        # /*fend*/ }}}
        #Redraw Raw Cog /*fold*/ {{{
        suppressWarnings(lines(x=cog$x*asperpix, y=-2.5*(log10(cog$y)-log10(ABvegaflux))+magZP, pch=20, col='grey',lty=2))
        # /*fend*/ }}}
        #Draw Deblended Cog /*fold*/ {{{
        suppressWarnings(lines(x=debl.cog$x*asperpix, y=-2.5*(log10(debl.cog$y)-log10(ABvegaflux))+magZP, pch=20, col='black',lty=2))
        # /*fend*/ }}}
        #Draw Sky Subtracted Cog /*fold*/ {{{
        suppressWarnings(lines(x=cog_nosky$x*asperpix, y=-2.5*(log10(cog_nosky$y)-log10(ABvegaflux))+magZP, pch=20, col='grey',lty=1))
        # /*fend*/ }}}
        #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
        suppressWarnings(lines(x=debl.cog_nosky$x*asperpix, y=-2.5*(log10(debl.cog_nosky$y)-log10(ABvegaflux))+magZP, pch=20, col='black',lty=1))
        # /*fend*/ }}}
        #Draw Legend /*fold*/ {{{
        legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended ApMag","Deblended ApMag"),lty=c(2,2,1,1,1,1),
               col=c('grey','black','grey','black','orange','green'),pch=-1, cex=0.5)
        # /*fend*/ }}}
        # /*fend*/ }}}
        # /*fend*/ }}}
      } else {
        #Plot in Flux Space /*fold*/ {{{
        #Plot Raw Cog /*fold*/ {{{
        ylim<-range(cog$y)
        #If a limit in ±Inf, plot around median /*fold*/ {{{
        if (!all(is.finite(ylim))) { ylim=median(cog$y, na.rm=TRUE)+c(-1E2,1E2) }
        magplot(x=cog$x*asperpix, y=cog$y, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Flux",ylim=ylim,main="Curve of Growth",type='l')
        # /*fend*/ }}}
        # /*fend*/ }}}
        #Draw Lines showing fluxes /*fold*/ {{{
        abline(h=dfaflux[i], lwd=1, col='green')
        abline(h=sfaflux[i], lwd=1, col='orange', lty=2)
        # /*fend*/ }}}
        #Redraw Raw Cog /*fold*/ {{{
        lines(x=cog$x*asperpix, y=cog$y, pch=20, col='grey',lty=2)
        # /*fend*/ }}}
        #Draw Deblended Cog /*fold*/ {{{
        lines(x=debl.cog$x*asperpix, y=debl.cog$y, pch=20, col='black',lty=2)
        # /*fend*/ }}}
        #Draw Sky Subtracted Cog /*fold*/ {{{
        lines(x=cog_nosky$x*asperpix, y=cog_nosky$y, pch=20, col='grey',lty=1)
        # /*fend*/ }}}
        #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
        lines(x=debl.cog_nosky$x*asperpix, y=debl.cog_nosky$y, pch=20, col='black',lty=1)
        # /*fend*/ }}}
        legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended Flux","Deblended Flux"),lty=c(2,2,1,1,1,1),
               col=c('grey','black','grey','black','orange','green'),pch=-1, cex=0.5)
        # /*fend*/ }}}
      }
      # /*fend*/ }}}
      #Plot the Deblended Image /*fold*/ {{{
      nc<-length(image_lims[i,1]:image_lims[i,2])
      nr<-length(image_lims[i,3]:image_lims[i,4])
      suppressWarnings(z<-matrix(magmap(dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],stretch='asinh')$map,ncol=nc,nrow=nr))
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apB, main="Image x Weight Matrix", asp=1, col=grey.colors(256), useRaster=FALSE, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims)
      # /*fend*/ }}}
      #Overlay the Aperture /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apT, main="Image x Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=FALSE, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims,add=TRUE))
      # /*fend*/ }}}
      #Draw the Axes and scalebar /*fold*/ {{{
      magaxis(frame.plot=T,main="Image x Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      #magbar("topright",col=rev(rainbow(256, start=0,end=2/3)), range=range(z),title='log(pixval)')
      # /*fend*/ }}}
      #Plot the Deblend Matrix /*fold*/ {{{
      z=dbw[[i]]
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apB, main="Weight Matrix", asp=1, col=grey.colors(256), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims)
      # /*fend*/ }}}
      #Overlay the Aperture /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apT, main="Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims,add=T))
      # /*fend*/ }}}
      #Draw the Axes and scalebar /*fold*/ {{{
      magaxis(frame.plot=T,main="Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      #magbar("topright",col=rev(rainbow(256, start=0,end=2/3)), range=c(0,1))
      # /*fend*/ }}}
      #Close the file /*fold*/ {{{
      dev.off()
      # /*fend*/ }}}
    }
    #Remove unneeded Arrays /*fold*/ {{{
    if (!makeresidmap) {
      sfabak<-NULL
    }
    rm(dbw)
    # /*fend*/ }}}
    #Notify /*fold*/ {{{
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #If map was input in Jy/bm we need to convert it back before output in SourceSubtraction /*fold*/ {{{
  if (Jybm) { ba=beamarea } else { ba=1. }
  # /*fend*/ }}}
  #If wanted, make the Residual Map /*fold*/ {{{
  if (makeresidmap) {
    if (!is.null(sfabak)) { sfa<-sfabak }
    if (filtcontam) {
      if (!quiet) { cat(paste("Writing Contaminant-subtracted Map to",nocontammap,"   ")) }
      #Perform Source Subtraction /*fold*/ {{{
      timer=system.time(sourcesubtraction(image.env$im,sfa,image_lims,dfaflux/ApCorr,file.path(pathroot,pathwork,pathout,nocontammap),image.env$hdr_str,ba,contams,diagnostic,verbose))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Contam Subtraction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
    if (!quiet) { cat(paste("Writing Source-subtracted Map to",residmap,"   ")) }
    # /*fend*/ }}}
    #Perform Source Subtraction /*fold*/ {{{
    timer=system.time(sourcesubtraction(image.env$im,sfa,image_lims,dfaflux/ApCorr,file.path(pathroot,pathwork,pathout,residmap),image.env$hdr_str,ba,insidemask,diagnostic,verbose))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Source Subtraction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #If Subtracting Contaminants, then remove them before output /*fold*/ {{{
  if ((writetab)&(filtcontam)) {
    id_g   <-id_g[which(contams==0)]
    ra_g   <-ra_g[which(contams==0)]
    dec_g  <-dec_g[which(contams==0)]
    theta_g<-theta_g[which(contams==0)]
    x_p    <-x_p[which(contams==0)]
    y_p    <-y_p[which(contams==0)]
    x_g    <-x_g[which(contams==0)]
    y_g    <-y_g[which(contams==0)]
    a_g    <-a_g[which(contams==0)]
    b_g    <-b_g[which(contams==0)]
    ssa    <-ssa[which(contams==0)]
    ssfa   <-ssfa[which(contams==0)]
    ssfa2  <-ssfa2[which(contams==0)]
    ssfad  <-ssfad[which(contams==0)]
    qssfad <-qssfad[which(contams==0),]
    spsf   <-spsf[which(contams==0)]
    ssfap  <-ssfap[which(contams==0)]
    ssfa2e2<-ssfa2e2[which(contams==0)]
    sfaflux<-sfaflux[which(contams==0)]
    skyflux<-skyflux[which(contams==0)]
    skylocal<-skylocal[which(contams==0)]
    skyerr <-skyerr[which(contams==0)]
    skyrms <-skyrms[which(contams==0)]
    skypval<-skypval[which(contams==0)]
    skyNBinNear <- skyNBinNear[which(contams==0)]
    skyNBinNear.mean  <- skyNBinNear.mean[which(contams==0)]
    skyflux.mean<-skyflux.mean[which(contams==0)]
    skyerr.mean   <- skyerr.mean[which(contams==0)]
    skylocal.mean <- skylocal.mean[which(contams==0)]
    skypval.mean  <- skypval.mean[which(contams==0)]
    skyrms.mean   <- skyrms.mean[which(contams==0)]
    ssfa2w <- ssfa2w[which(contams==0)]
    ssfadw <- ssfadw[which(contams==0)]
    ssfae  <- ssfae[which(contams==0)]
    ssfae2 <- ssfae2[which(contams==0)]
    ssfaw  <- ssfaw[which(contams==0)]
    detecthres<-detecthres[which(contams==0)]
    detecthres.mag<-detecthres.mag[which(contams==0)]
    sfaerr <-sfaerr[which(contams==0)]
    sdfa   <-sdfa[which(contams==0)]
    sdfa2  <-sdfa2[which(contams==0)]
    sdfad  <-sdfad[which(contams==0)]
    qsdfad <-qsdfad[which(contams==0),]
    if (iterateFluxes) {
      fluxiters<-fluxiters[which(contams==0),]
      erriters<-erriters[which(contams==0),]
      sdfaiters<-sdfaiters[which(contams==0),]
    }
    sdfae <-sdfae[which(contams==0)]
    sdfae2<-sdfae2[which(contams==0)]
    sdfa2e2<-sdfa2e2[which(contams==0)]
    dfaflux<-dfaflux[which(contams==0)]
    deblerr <-deblerr[which(contams==0)]
    dfaerr <-dfaerr[which(contams==0)]
    saturated<-saturated[which(contams==0)]
    pixflux<-pixflux[which(contams==0)]
    ApCorr <-ApCorr[which(contams==0)]
    WtCorr <-WtCorr[which(contams==0)]
    stamplen<-stamplen[which(contams==0)]
    mags   <-mags[which(contams==0)]
    if (length(fluxweight!=1)) { fluxweight<-fluxweight[which(contams==0)] }
    if (RanCor) { randoms<-randoms[which(contams==0),] }
    if (BlankCor) { blanks<-blanks[which(contams==0),] }
    contams <-contams[which(contams==0)]
  }# /*fend*/ }}}

  #Create Photometry Warning Flags /*fold*/ {{{
  photWarnFlag<-rep("",length(ra_g))
  #Bad Quartered Photometry /*fold*/ {{{
  Qbad<-(apply(qsdfad,1,'max',na.rm=T)/apply(qsdfad,1,'sum',na.rm=T))>=0.7
  photWarnFlag<-paste0(photWarnFlag,ifelse(Qbad,"Q",""))
  # /*fend*/ }}}
  #Saturation /*fold*/ {{{
  photWarnFlag<-paste0(photWarnFlag,ifelse(saturated,"X",""))
  # /*fend*/ }}}
  #Bad Sky Estimate /*fold*/ {{{
  if (doskyest|getskyrms) {
    photWarnFlag<-paste0(photWarnFlag,ifelse(skyNBinNear<=3,"S",""))
  }
  # /*fend*/ }}}
  #Iterative Deblend Warning /*fold*/ {{{
  if (iterateFluxes) {
    photWarnFlag<-paste0(photWarnFlag,ifelse(iterateLost,"I",""))
  }
  # /*fend*/ }}}
  # /*fend*/ }}}

  #If wanted, output the Results Table /*fold*/ {{{
  if (writetab) {
    #Output the results table /*fold*/ {{{
    if ((nloops!=1)&&(length(param.env$tableoutname)!=nloops)) {
      if (!quiet) { cat(paste('Writing Results Table to ',file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,"_file",f,".csv")),'   ')) }
      timer=system.time(writesfatableout(filename=file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,"_file",f,".csv"))) )
    } else {
      if (!quiet) { cat(paste('Writing Results Table to ',file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,".csv")),'   ')) }
      timer=system.time(writesfatableout(filename=file.path(pathroot,pathwork,pathout,paste(sep="",tableoutname,".csv"))) )
    }
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Write Results Table - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }# /*fend*/ }}}
  #-----Diagnostic-----# /*fold*/ {{{
  if (diagnostic) { message(paste('Sum of PSF = ',beamarea)) }
  # /*fend*/ }}}
  #If wanted, open the browser for user inspection /*fold*/ {{{
  if (interact) {
     cat(paste("Launching Interactive Mode: To end, type 'c'\n"))
     sink(type='message')
     browser()
     sink(sinkfile, type='message')
  }# /*fend*/ }}}
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART SIX: FINISH /*fold*/ {{{
  #Send Parameters to logfile /*fold*/ {{{
  sink(sinkfile, type="output")
  cat("Memory Hogs in this run:\n")
  print(lsos(envir=environment(), head=TRUE, n=10))
  cat("Images used in this run:\n")
  print(lsos(envir=image.env, head=FALSE))
  sink(type="output")
  # /*fend*/ }}}
  #Return /*fold*/ {{{
  if (!quiet) { cat('\n') }
  #on.exit(detach(image.env),add=TRUE)
  if (!is.null(env)) {
    on.exit(detach(env), add=TRUE)
  }
  if (!diagnostic) {
    return=list(SFAflux=sfaflux,SFAerror=sfaerr, DFAflux=dfaflux,DFAerror=dfaerr)
  } else {
    return=list(SFAflux=sfaflux,SFAerror=sfaerr, DFAflux=dfaflux,DFAerror=dfaerr, SA_Stamps=sa, SFA_Stamps=sfa, WSFA_Stamps=wsfa, DFA_Stamps=dfa)
  }
  # /*fend*/ }}}
  # /*fend*/ }}}

}
#Func /*fold*/ {{{
.executeRanCor<-function() {
cat("Executing RanCor...\n"); Sys.sleep(2); cat('\n    |   _______   _     _   _     _   _ __      _    |\n    |  |__   __| | |   | | | |   | | |  __ \\   | |   |\n    |     | |    | |___| | | |   | | | |  | |  | |   |\n    |     | |    |  ___  | | |   | | | |  | |  |_|   |\n    |     | |    | |   | | | |___| | | |__| |   _    |\n    |     |_|    |_|   |_| |_______| |_____/   |_|   |\n    |                                                |\n     \\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\n      \\/  \\/  \\/  \\/  ,------------,  \\/  \\/  \\/  \\/\n         ____        ({ XX      XX })        ____\n        /////|\\      _|   \\\\  //   |_      /|\\\\\\\\\\\n        VVVV | \\____/ { \\  \\\\//  / } \\____/ | VVVV\n        ////_|       ,|V=V=V=V=V=V=|,       |_\\\\\\\\\n     ___\\\\\\\\/_\\~~~~~/_{+^+^+^+^+^+^}_\\~~~~~/_\\////___\n\n'); Sys.sleep(2); cat("... Done, you heartless Jedi Scum\n\n"); Sys.sleep(2)
}
# /*fend*/ }}}
