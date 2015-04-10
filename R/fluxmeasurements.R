fluxmeasurements <-
function(env=NULL) {
#Procedure measures the fluxes present in a supplied image
#inside catalogued apertures.

  #PART ZERO:  PSF DETERMINATION & GENERAL INITIALISATION {{{
  # Load Parameter Space {{{
  if (!is.null(env)) {
    attach(env)
  }
  #}}}

  #Set function Environments {{{
  environment(make_sa_mask)<-environment()
  environment(make_sfa_mask)<-environment()
  environment(make_a_mask)<-environment()
  environment(make_data_mask)<-environment()
  environment(readpsf)<-environment()
  environment(make_deblended_weightmap)<-environment()
  environment(writesfatableout)<-environment()
  environment(writefitsout)<-environment()
  environment(sourcesubtraction)<-environment()
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
    psf<-readpsf(outenv=environment(),file.path(pathroot,pathwork,psfmap),asperpix,defbuff*max(a_g),confidence,gauss_fwhm_as=gauss_fwhm_as)
    #}}}
    #Notify {{{
    if (verbose) { message(paste("Maxima of the PSF is at pixel", which(psf == max(psf)),"and has value",max(psf))) }
    #}}}
    #Get radius of FWHM using FWHM confidence value = erf(2*sqrt(2*log(2))/sqrt(2)) {{{
    psffwhm<-get.fwhm(psf)
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
    #If no PSF, set FWHM to median aperture radius + buffer {{{
    psffwhm<-round(min(a_g[which(a_g>0)])/asperpix)
    if ((is.na(psffwhm))|(!is.finite(psffwhm))) {
      #All Apertures are point sources, force width = 5 pixels {{{
      #Details {{{
      #There are no apertures that are not point sources,
      #and no PSF convolution - everything is a single pixel flux
      #set to a default of 5pixels }}}
      psffwhm<-5
      #}}}
    }
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
    if (length(image.env$ime)>1) { image.env$ime<-image.env$ime/beamarea }
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
  #Get StampLengths for every galaxy {{{
  stamplen<-(floor((ceiling(defbuff*a_g*2/asperpix)+ceiling(psf.clip))/2)*2+5)
  #}}}
  #Discard any contaminants that are beyond their nearest neighbours stamp {{{
  #Number of Nearest Neighbors to search {{{
  if (!exists("nNNs")) { nNNs<-10 }
  if (!exists("checkContam")) { checkContam<-TRUE }
  #}}}
  if (filtcontam & checkContam) {
    timer<-proc.time()
    if (!quiet) { cat("Removing Contaminants that are irrelevant ") }
    catlen<-length(x_g)
    nearest<-nn2(data.frame(x_g[which(contams==0)],y_g[which(contams==0)]),data.frame(x_g[which(contams==1)],y_g[which(contams==1)]),k=nNNs)
    contam.inside<-NULL
    for(ind in 1:length(which(contams==1))) {
      #If any of the nearest neighbors overlap with the contaminant
      contam.inside<-c(contam.inside,any(nearest$nn.dist[ind,] < (sqrt(2)/2*(stamplen[which(contams==0)][nearest$nn.idx[ind,]]+stamplen[which(contams==1)][ind]))))
    }
    insidemask<-rep(TRUE,catlen)
    insidemask[which(contams==1)]<-contam.inside
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
    if (verbose) { message(paste("There are",length(x_g),"supplied objects & contaminants that intersect (",
                                  round(((catlen-length(x_g))/catlen)*100, digits=2),"% of objects were contaminants that didn't intersect galaxies)")) }
    #}}}
    if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      message(paste('Remove irrelevant contaminants - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
  #}}}
  #Convert decimal pixel values into actual pixel values {{{
  x_p<-floor(x_g)
  y_p<-floor(y_g)
  #}}}

  #Create an array of stamps containing the data image subsection for all objects {{{
  timer=system.time(im_mask<-make_data_mask(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Data Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #Remove Image arrays as they are no longer needed {{{
  imenvlist<-ls(envir=image.env)
  imenvlist<-imenvlist[which(imenvlist!="im"&imenvlist!="hdr_str")]
  rm(list=imenvlist, envir=image.env)
  #}}}
  #Create an array of stamps containing the apertures for all objects {{{
  timer=system.time(sa<-make_sa_mask(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make SA Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #Discard any apertures that were zero'd in the make process {{{
  totsa<-foreach(sam=sa, i=1:length(sa), .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) } ; sum(sam) }
  insidemask<-(totsa > 0)
  if ((nloops<=1)&(length(which(insidemask==TRUE))==0)) { sink(type="message") ; stop("No Single Apertures are inside the Mask.") }  # Nothing inside the mask
  else if (length(which(insidemask==TRUE))==0) {
      warning("No Single Apertures are inside the image.")
      #Notify & Close Logfile {{{
      if (!is.null(env)) {
        on.exit(detach(env), add=TRUE)
      }
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      if (!quiet) {
        cat(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures in the Mask\n'))
      }
      return(NULL)
      #}}}
    }
  sa<-sa[which(insidemask)]
  im_mask<-im_mask[which(insidemask)]
  if (length(imm_mask)>1) { imm_mask<-imm_mask[which(insidemask)] }
  if (length(ime_mask)>1) { ime_mask<-ime_mask[which(insidemask)] }
  stamplen<-stamplen[which(insidemask)]
  stamp_lims<-rbind(stamp_lims[which(insidemask),])
  mstamp_lims<-rbind(mstamp_lims[which(insidemask),])
  estamp_lims<-rbind(estamp_lims[which(insidemask),])
  #Remove object that failed aperture creation {{{
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
  if (filtcontam) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]
  chunkSize=ceiling(length(id_g)/getDoParWorkers())
  mpiopts<-list(chunkSize=chunkSize)
  #}}}
  #}}}
  #Re-Initialise object count {{{
  npos<-length(id_g)
  #}}}
  #Create a full mask of all apertures in their correct image-space locations {{{
  timer=system.time(image.env$aa<-make_a_mask(outenv=environment(), sa, dim(image.env$im)))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
   message(paste('Make A Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #If wanted, output the All Apertures Mask {{{
  if (makeaamask) {
    if (!quiet) { cat(paste('Outputting All Apertures Mask to',aafilename,"   ")) }
    #Write All Apertures Mask to file
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,aafilename),image.env$aa,image.env$hdr_str,nochange=TRUE))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }#}}}
  #}}}
  #PART TWO:  SINGLE FILTERED APERTURE MASKS {{{
  #Details {{{
  #If wanted, perform the convolution of the
  #apertures with the PSF. If not wanted,
  #duplicate the single apertures.
  #Also, perform flux weighting. }}}
  #If wanted, use pixel fluxes as fluxweights {{{
  if (usePixelFluxWeights) {
    #Image Flux at central pixel {{{
    cat("Determine Image Flux at central pixel ")
    pixflux<-foreach(xp=x_p-(im_stamp_lims[,1]-1),yp=y_p-(im_stamp_lims[,3]-1), im=im_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          im[xp,yp]
    }
    pixflux<-array(unlist(pixflux),dim=c(dim(pixflux[[1]]),length(pixflux)))
    if (verbose) { cat(" - Done\n") }
    #}}}
    #Determine Noise Characteristics {{{
    if (verbose) { cat("Determine Image Noise Characteristics ") }
    if (any(pixflux==-99)) { warning("Some pixel flux determinations failed"); pixflux[which(pixflux==-99)]<-NA }
    if (length(image.env$im) > 1E6) { index<-runif(1E6,min=1,max=length(image.env$im)) } else { index<-1:length(image.env$im) }
    xdat<-as.numeric(image.env$im[index])
    ind<-which(xdat < mean(xdat)*1.1)
    if (length(ind)>2) {
      x.dens=density(xdat[ind])
      x.mode=x.dens$x[which.max(x.dens$y)]
      # Select all data with value < mode - this will be exclusively noise
      xdat=xdat-x.mode
      noise=xdat[which(zapsmall(xdat) <= 0.0)]
      # Determine Gaussian Characteristics
      stdev=abs(quantile(noise,pnorm(-1)*2))
      message("Noise Profile: Mean = ", x.mode, " ; Stdv = ",abs(stdev)," ;\n")
      # Set weighting limits as all pixels beyond 3sigma of noise
      pixlimit<-x.mode+2*stdev
      if (pixlimit<0) { warning("Pixel Limit is < 0?! Using absolute value"); pixlimit<-abs(x.mode+2*stdev) }
      if (pixlimit==0) { warning("Pixel Limit is == 0?! Using next highest value as limit"); pixlimit<-min(pixflux[which(pixflux>pixlimit)]) }
      if (verbose) { cat(" - Done\n") }
      #}}}
      #Create Flux weights {{{
      if (verbose) { cat("Creating Pixel Flux Weightings ") }
      #hi=quantile(pixflux, pnorm(4), na.rm=TRUE)
      hi=max(pixflux, na.rm=TRUE)
      if (hi<pixlimit) {
        warning("Strangeness in pixel-flux weighting, using 1sigma lower limit")
        pixlimit<-x.mode+stdev
        if ((hi<pixlimit)&max(pixflux, na.rm=TRUE)>pixlimit) {
          warning("Pixel-flux weighting still bad. Using pixel mode as lower limit")
          pixlimit<-x.mode
        } else if ((hi<pixlimit)&max(pixflux, na.rm=TRUE)<pixlimit) {
          warning("Lower limit determination must be bad. Using pixel-flux minima as lower limit, and pixel-flux maxima as upper limit")
          pixlimit=min(pixflux, na.rm=TRUE)
          hi=max(pixflux, na.rm=TRUE)
        }
      }
      fluxweight<-magmap(pixflux, lo=pixlimit, hi=hi, range=c(0.0001, 1), bad=pixlimit, type="num", stretch='lin')$map
      rm(index)
      rm(xdat)
      rm(noise)
      rm(x.dens)
    } else {
      fluxweight<-magmap(pixflux, lo=min(pixflux), hi=max(pixflux), range=c(0.0001,1), type="num", stretch='lin')$map
    }
    cat(" - Done\n")
    #}}}
  }#}}}
  #If needed, do Convolutions & Fluxweightings {{{
  if ((!psffilt)&(length(which(fluxweight!=1))!=0)) {
    #No Convolution, Need Fluxweighting {{{
    #Details {{{
    #If not convolving with psf, and if all the fluxweights are not unity,
    #then skip filtering, duplicate the arrays, and only weight the stamps/apertures
    #}}}
    #No Convolution, duplicate apertures {{{
    if (verbose) { message("No Convolution: Convolved Apertures are identical to Simple Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    #}}}
    #Perform aperture weighting {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(outenv=environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of stamps/apertures {{{
    timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #}}}
  } else if ((!psffilt)&(length(which(fluxweight!=1))==0)) {
    #No Convolution, No Fluxweighting {{{
    #Details {{{
    #If not convolving with psf, and if all the fluxweights are unity,
    #then skip filtering and weighting, and simply duplicate the arrays }}}
    #Duplicate arrays {{{
    if (verbose) { message("No Convolution: Convolved Apertures are identical to Simple Apertures")
                   message("No Fluxweights: Weighted Convolved Apertures are identical to Convolved Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    wsfa<-sfa
    image.env$wfa<-image.env$fa
    #}}}
    #}}}
  } else if ((psffilt)&(length(which(fluxweight!=1))!=0)) {
    #Convolving PSF, Need Fluxweighting {{{
    #Details {{{
    #If convolving with psf, and if all the fluxweights are not unity,
    #then colvolve and weight the stamps/apertures }}}
    #Perform Convolution {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(outenv=environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all convolved stamps/apertures {{{
    timer=system.time(image.env$fa<-make_a_mask(outenv=environment(), sfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Perform aperture weighting {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make_sfa_mask(outenv=environment(), sa,fluxweightin=fluxweight))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all the convolved weighted stamps/apertures {{{
    timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #}}}
  } else if ((psffilt)&(length(which(fluxweight!=1))==0)) {
    #Convolve PSF, No Fluxweighting {{{
    #Details {{{
    #If convolving with psf, and if all the fluxweights are unity,
    #colvolve, and then duplicate the stamps/apertures }}}
    #Perform Convolution {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make_sfa_mask(outenv=environment(), sa))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Create a full mask of all convolved stamps/apertures {{{
    timer=system.time(image.env$fa<-make_a_mask(outenv=environment(), sfa, dim(image.env$im)))
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
  #Remove Uneeded Image {{{
  rm(aa, envir=image.env)
  gc()
  #}}}
  #-----Diagnostic-----# {{{
  if (diagnostic) {
    if (verbose) { message("Checking Apertures for scaling errors") }
    foreach(sam=sa, sfam=sfa, i=1:length(sa), .combine=function(a,b){NULL},.inorder=FALSE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) }
      if (max(sfam)>1){warning(paste("Max of Convolved Aperture",i,"is",max(sfam))) }
    }
  }#}}}
  #Do we want to plot a sample of the apertures? {{{
  if (plotsample) {
    #Set output name {{{
    pdf(file.path(pathroot,pathwork,pathout,"PSF_Samples.pdf"))
    #}}}
    #Set Layout {{{
    par(mfrow=c(2,2))
    #}}}
    #Output a random 15 apertures to file {{{
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
    }#}}}
    #Close the file {{{
    dev.off()
    #}}}
  }
  #}}}
  #Remove arrays that are no longer needed {{{
  #rm(sa)
  #gc()
  #}}}
  #If wanted, output the Convolved & Weighted Aperture Mask {{{
  if (makefamask) {
    if (!quiet) { cat(paste('Outputting All Convolved Apertures Mask to',fafilename,"   ")) }
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,fafilename),image.env$wfa,image.env$hdr_str,nochange=TRUE) )
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
    #         b) set all nonzero locations in fa to 0 in sm array }}}
    #Make Sourcemask {{{
    if (!exists("TransmissionMap")) { TransmissionMap<-FALSE }
    if (!quiet) { cat(paste("Creating Sourcemask    ")) }
    if (length(image.env$imm)>1) { sm<-image.env$imm } else { sm<-array(1, dim=dim(image.env$im)) }
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
    if (doskyest||getskyrms) {
      imm_mask<-list(NULL)
      for (i in 1:npos) {
        imm_mask[[i]]<-sm[imm_stamp_lims[i,1]:imm_stamp_lims[i,2],imm_stamp_lims[i,3]:imm_stamp_lims[i,4]]
      }
    }
    if (!quiet) { cat(" - Done\n") }
    #}}}
    #-----Diagnoistic-----# {{{
    if (diagnostic) {
      message(paste("SourceMask Max/Min:",max(sm),min(sm)))
      message(paste("OLDMethod - SourceMask Max/Min:",max(1-image.env$fa),min(1-image.env$fa)))
    }#}}}
    #If wanted, output the SourceMask {{{
    if (!exists("sourcemaskout")) { sourcemaskout<-FALSE }
    if (sourcemaskout){
      if (!quiet) { cat(paste('Outputting Source Mask to',smfilename,"   ")) }
      writefitsout(file.path(pathroot,pathwork,pathout,smfilename),sm,image.env$hdr_str,nochange=TRUE)
      if (!quiet) { cat(" - Done\n") }
    }#}}}
    #If we want the SourceMask only, then end here {{{
    if (sourcemaskonly) {
      if (!quiet) { cat("SourceMaskOnly Flag Set\n")  }
      return()
    } #}}}
  } else if (plotsample) {
    #Set the mask to 1 everywhere {{{
    sm<-1
    imm_mask<-1
    #}}}
  }#}}}
  #Remove arrays that are no longer needed {{{
  rm(fa, envir=image.env)
  gc()
  #}}}
  #}}}
  #PART THREE:  DEBLENDING {{{
  #If wanted, perform sky estimation. Otherwise set to NA {{{
  if (doskyest||getskyrms) {
    #Get sky estimates {{{
    if (!quiet) { message("Perfoming Sky Estimation"); cat("Performing Sky Estimation") }
    #Perform Sky Estimation {{{
    #timer<-system.time(skyest<-skyback(ra_g,dec_g,cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,origim=list(dat=list(im)),maskim=list(dat=list(sm)),
    #                astrom=astr_struc,clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm))
    timer<-system.time(skyest<-skyback.par(x_p-(im_stamp_lims[,1]-1),y_p-(im_stamp_lims[,3]-1),cutlo=(a_g/asperpix),cuthi=(a_g/asperpix)*5,im_mask=im_mask,imm_mask=imm_mask,
                    clipiters=skycutiters,probcut=skyprobcut,PSFFWHMinPIX=psffwhm, mpiopts=mpiopts))
    #}}}
    #Get sky parameters {{{
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
    skylocal[which(is.na(skylocal))]<-skydefault
    #skyflux<-skylocal*sdfa
    #skyerr<-skyerr*sdfa
    #}}}
    #}}}
    #Notify {{{
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Sky Estimate - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    #}}}
    #Calculate Detection Thresholds {{{
    if (!quiet) { message("Calculating Detection Limits"); cat("Calculating Detection Limits") }
    detecthres<-foreach(sfam=sfa, srms=skyrms, .inorder=TRUE, .combine='c', .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { 5*srms*sqrt(length(which(sfam>0))) }
    if (Magnitudes) {
      suppressWarnings(detecthres.mag<--2.5*(log10(detecthres)-log10(ABvegaflux))+magZP)
    } else {
      detecthres.mag<-array(NA, dim=c(length(sfa)))
    }
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    #}}}
    #If wanted, do sky subtraction {{{
    #if (doskyest) {
    #  if (!quiet) { message("Perfoming Sky Subtraction"); cat("   Performing Sky Subtraction") }
    #  #Subrtract Sky Flux
    #  dfaflux<-dfaflux-skyflux
    #  sfaflux<-sfaflux-skyflux
    #  dfaerr[which(!is.na(skyerr))]<-sqrt(dfaerr^2+skyerr^2)[which(!is.na(skyerr))]
    #  sfaerr[which(!is.na(skyerr))]<-sqrt(sfaerr^2+skyerr^2)[which(!is.na(skyerr))]
    #  if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    #}#}}}
  } else {
    skylocal<-array(NA, dim=c(length(sfa)))
    skyflux<-array(NA, dim=c(length(sfa)))
    skyerr<-array(NA, dim=c(length(sfa)))
    skyrms<-array(NA, dim=c(length(sfa)))
    skypval<-array(NA, dim=c(length(sfa)))
    detecthres<-array(NA, dim=c(length(sfa)))
    detecthres.mag<-array(NA, dim=c(length(sfa)))
  }#}}}
  #Perform Deblending of apertures {{{
  timer=system.time(dbw<-make_deblended_weightmap(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Deblended Weightmap - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #Remove arrays that are no longer needed {{{
  rm(wfa, envir=image.env)
  #}}}
  #Finalise SFA Apertures {{{
  sfabak<-NULL
  if (!exists("PSFWeighted")) { PSFWeighted<-FALSE }
  if (psffilt & !PSFWeighted) {
    # If we have convolved with a PSF, & not doing PSF Weighting, convert back to tophat apertures (now expanded by PSF convolution) {{{
    #If we want the residual map, save the model apertures {{{
    if ((plotsample)|(makeresidmap)) {
      sfabak<-sfa
    }
    #}}}
    # For each aperture, Binary filter at aplim*max(ap) {{{
    sba<-foreach(sfam=sfa, .options.mpi=mpiopts, .noexport=ls(envir=environment()))%dopar%{
      apvals<-rev(sort(sfam))
      tempsum<-cumsum(apvals)
      tempfunc<-approxfun(tempsum,apvals)
      apLim<-tempfunc(apLimit*max(tempsum, na.rm=TRUE))
      #message(paste("ApLimit:",apLimit,"; PSFLimit:",psfLimit))
      sfam[which(sfam <  apLim)]<-0
      sfam[which(sfam >= apLim)]<-1
      return=sfam
    }
    #}}}
    #Update the sfa {{{
    sfa<-sba
    rm(sba)
    #}}}
    #If wanted, output the Convolved & Weighted Aperture Mask {{{
    if (makefamask) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), sfa, dim(image.env$im)))
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
    }#}}}
    #remove unneeded variables {{{
    rm(psfvals)
    rm(tempsum)
    rm(tempfunc)
    #}}}
    #}}}
  } else if (psffilt & PSFWeighted) {
    #If we have convolved by the PSF & want PSF Weighting {{{
    #If we want the residual map, save the model apertures {{{
    if ((plotsample)|(makeresidmap)) {
      sfabak<-sfa
    }
    #}}}
    # For each resolved aperture, Binary filter the section out to the desired level, then fuzzy the edges {{{
    chunkSize=ceiling(length(which(a_g>0))/getDoParWorkers())
    tempopts<-list(chunkSize=chunkSize)
    sba<-foreach(sfam=sfa[which(a_g>0)], .options.mpi=tempopts, .noexport=ls(envir=environment()))%dopar%{
      apvals<-rev(sort(sfam))
      tempsum<-cumsum(apvals)
      tempfunc<-approxfun(tempsum,apvals)
      apLim<-tempfunc(apLimit*max(tempsum, na.rm=TRUE))
      #message(paste("ApLimit:",apLimit,"; PSFLimit:",psfLimit))
      sfam[which(sfam >= apLim)]<-1
      sfam[which(sfam <  apLim)]<-sfam[which(sfam <  apLim)]/apLim
      return=sfam
    }
    sfa[which(a_g>0)]<-sba
    rm(sba)
    #If wanted, output the Convolved & Weighted Aperture Mask {{{
    if (makefamask) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make_a_mask(outenv=environment(), sfa, dim(image.env$im)))
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
    }#}}}
    #remove unneeded variables {{{
    rm(psfvals)
    rm(tempsum)
    rm(tempfunc)
    #}}}
    #}}}
    #}}}
  }
  #}}}
  #If PSF Supplied & sample wanted, Plot PSF and Min Aperture Correction {{{
  if (!(nopsf)&(plotsample)) {
    #PSF with Contours {{{
    pdf(file.path(pathroot,pathwork,pathout,"PSF.pdf"))
    psfvals<-rev(sort(psf))
    tempsum<-cumsum(psfvals)
    tempfunc<-approxfun(tempsum,psfvals)
    psfLimit<-tempfunc(apLimit*max(tempsum, na.rm=TRUE))
    suppressWarnings(image(log10(psf),main="PSF & Binary Contour Levels", asp=1,col=heat.colors(256),useRaster=TRUE))
    contour(psf, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
    if (psffilt) { contour(psf, levels=c(psfLimit), labels=c(apLimit), col='green', add=TRUE) }
    dev.off()
    #}}}
    if (!PSFWeighted) {
      #Plot Example of PS minimum aperture Correction{{{
      ind<-(which(a_g==0)[1])
      #If no point sources, use the minimum aperture
      if (is.na(ind)) { ind<-which.min(a_g) }
      #Make Sure PSF is centred on centre of stamp
      centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
      delta<-floor(stamplen[ind]/2)*c(-1,+1)
      lims<-rbind(centre[1]+delta,centre[2]+delta)
      dx<-lims[1,1]-1
      dy<-lims[2,1]-1
      #Reinterpolate the PSF at point source XcenYcen {{{
      lenxpsf<-length(psf[,1])
      lenypsf<-length(psf[1,])
      lenx<-length(1:stamplen[ind])
      leny<-length(1:stamplen[ind])
      #Make grid for psf at old pixel centres {{{
      psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
      #}}}
      #Make expanded grid of new pixel centres {{{
      expanded<-expand.grid(seq(1,lenx),seq(1,leny))
      xnew<-expanded[,1]-x_g[ind]%%1
      ynew<-expanded[,2]-y_g[ind]%%1
      #}}}
      #Interpolate {{{
      ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
      #}}}
      #}}}
      pdf(file.path(pathroot,pathwork,pathout,"ApertureCorrection.pdf"))
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction", asp=1,col=heat.colors(256),useRaster=TRUE))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      if (psffilt) { contour(psf, levels=c(psfLimit), labels=c(apLimit), col='green', add=TRUE) }
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=TRUE))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round((1+((spsf-ssfap)/spsf)),digits=2),sep=""))
      dev.off()
      #}}}
    }
  }
  #}}}
  #Generate Deblended Flux Arrays {{{
  if (!quiet) { cat("Generating Deblended Flux Arrays") }
  timer<-proc.time()
  #Create the Deblended Flux Array {{{
  dfa<-foreach(dbwm=dbw, sfam=sfa, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
  #}}}
  #Notify {{{
  if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Make DFA - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   - Done\n") }
  #}}}
  #}}}
  #If wanted; Iterate the fluxes to improve final measurements {{{
  if (iterateFluxes) {
    #Notify {{{
    message(paste("Iterating Flux Determination",nIterations,"times {\n"))
    if (!quiet) { cat("Iterating Flux Determination",nIterations,"times {\n") }
    #}}}
    #For the number of desired iterations
    weightType='flux'
    quietbak<-quiet
    fluxiters<-matrix(as.numeric(NA),ncol=nIterations,nrow=length(im_mask))
    erriters<-fluxiters
    sdfaiters<-fluxiters
    for (iter in 1:nIterations) {
      #Calculate the flux per object {{{
      #Notify {{{
      message(paste("Calculating Flux (#",iter,")"))
      if (!quiet) { cat("  Calculating Flux (#",iter,")") }
      #}}}
      timer<-proc.time()
      sdfad<-foreach(dfam=dfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
         sum(dfam*(im[xlo:xup,ylo:yup]))
      }
      sdfad<-array(unlist(sdfad),dim=c(dim(sdfad[[1]]),length(sdfad)))
      #Notify {{{
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      #}}}
      #}}}
      #Calculate the error per object {{{
      #Notify {{{
      message(paste("Calculating Error (#",iter,")"))
      if (!quiet) { cat("  Calculating Error (#",iter,")") }
      timer<-proc.time()
      #}}}
      if (length(ime_mask)>1) {
        sdfae2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
            sum((dfam*ime[xlo:xup,ylo:yup])^2.)
        }
        sdfae2<-array(unlist(sdfae2),dim=c(dim(sdfae2[[1]]),length(sdfae2)))
      } else if (length(ime_mask)==1){
        if (ime_mask==1) {
          sdfae2<-sdfa2
        } else if (ime_mask==0) {
          sdfae2<-rep(0,length(sdfad))
        } else {
          sdfae2<-foreach(dfam=dfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
              sum((dfam*ime_mask)^2.)
          }
          sdfae2<-array(unlist(sdfae2),dim=c(dim(sdfae2[[1]]),length(sdfae2)))
        }
      } else {
        sink(type="message")
        stop("Image Error is NULL")
      }
      #Notify {{{
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      #}}}
      #}}}
      #If wanted, remove sky estimates {{{
      if (!quiet) { message(paste("Calculating Deblend Fraction (#",iter,")")); cat(paste("  Calculating Deblend Fraction (#",iter,")")) }
      timer<-proc.time()
      sdfa<-foreach(dfam=dfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
      sdfa<-array(unlist(sdfa),dim=c(dim(sdfa[[1]]),length(sdfa)))
      #Subrtract Sky Flux
      if (doskyest) {
        sdfad<-sdfad-skylocal*sdfa
      }
      if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
        message(paste('Sky Estimate - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat(" - Done\n") }
      #}}}
      #Save the values of the fluxes for output {{{
      fluxiters[,iter]<-sdfad
      erriters[,iter]<-sdfae2
      sdfaiters[,iter]<-sdfa
      #}}}
      #Re-calculate the weighted apertures {{{
      #Notify {{{
      message(paste("Calculating Weighting Apertures (#",iter,")"))
      if (!quiet) { cat("  Calculating Weighting Apertures (#",iter,")") }
      #}}}
      quiet<-TRUE
      timer=system.time(wsfa<-make_sfa_mask(outenv=environment(), sfa,fluxweightin=sdfad))
      timer2=system.time(image.env$wfa<-make_a_mask(outenv=environment(), wsfa, dim(image.env$im)))
      quiet<-quietbak
      #Notify {{{
      if (showtime) { cat(" - Done (",round(timer[3]+timer2[3],digits=2),"sec )\n")
      } else if (!quiet) { cat(" - Done\n") }
      #}}}
      #}}}
      #Re-calculate the deblending matricies {{{
      #Notify {{{
      message(paste("Calculating Deblend Matricies (#",iter,")"))
      if (!quiet) { cat("  Calculating Deblend Matricies (#",iter,")") }
      #}}}
      quiet<-TRUE
      timer2=system.time(dbw<-make_deblended_weightmap(outenv=environment()))
      quiet<-quietbak
      timer<-proc.time()
      dfa<-foreach(dbwm=dbw, sfam=sfa, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
      #Notify {{{
      if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3]+timer2[3],digits=2),"sec )\n")
      } else if (!quiet) { cat("   - Done\n") }
      #}}}
      #}}}
    }
    if (!quiet) { cat("} Done\n") }
  }
  #}}}
  #If wanted, output the Deblended & Convolved Aperture Mask {{{
  if (makedfamask) {
    if (!quiet) { cat(paste('Making All Deblended Convolved Apertures Mask - ')) }
    timer=system.time(image.env$adfa<-make_a_mask(outenv=environment(), dfa, dim(image.env$im)))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (!quiet) { cat(paste('Outputting All Deblended Convolved Apertures Mask to',dfafilename,"   ")) }
    timer=system.time(writefitsout(file.path(pathroot,pathwork,pathout,dfafilename),image.env$adfa,image.env$hdr_str,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }#}}}
  #Remove array that is no longer needed {{{
  if (!plotsample) { rm(dbw) }
  rm(adfa, envir=image.env)
  gc()
  #}}}
  #If wanted, Perform Randoms Correction {{{
  if (!exists("RanCor")) { RanCor<-FALSE }
  if (RanCor=="execute") {
    .executeRanCor()
    RanCor<-TRUE
  }
  if (RanCor) {
    if (!quiet) { message("Perfoming Randoms Correction"); cat("Performing Randoms Correction") }
    #debug(rancor.par)
    timer<-system.time(randoms<-rancor.par(im_mask=im_mask,imm_mask=imm_mask,ap_mask=sfa,stamplims=stamp_lims,masklims=mstamp_lims,numIters=nRandoms,mpiopts=mpiopts))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Randoms Correction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
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
    if (length(which(is.na(ime_mask)))!=0) {        message(paste("Input Parameter 'ime' contains NA elements")) }
    if (length(which(is.na(sfa)))!=0) {        message(paste("Input Parameter 'sfa' contains NA elements")) }
    if (length(which(is.na(dfa)))!=0) {        message(paste("Input Parameter 'dfa' contains NA elements")) }
    if (length(which(is.na(fluxcorr)))!=0) {        message(paste("Input Parameter 'fluxcorr' contains NA elements")) }
  }#}}}
  #Perform Calculations {{{
  if (!quiet) { cat("   Calculating Fluxes ") }
  if (verbose) { cat("{\n") }
  timer<-proc.time()
#-----
  #Image Flux at central pixel; pixflux {{{
  if (!usePixelFluxWeights) {
    if (verbose) { cat("      Image Flux at central pixel ") }
    pixflux<-foreach(xp=x_p-(im_stamp_lims[,1]-1),yp=y_p-(im_stamp_lims[,3]-1), im=im_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          im[xp,yp]
    }
    pixflux<-array(unlist(pixflux),dim=c(dim(pixflux[[1]]),length(pixflux)))
    if (verbose) { cat(" - Done\n") }
  }
  #}}}
#-----
  #Integral of the aperture; ssa {{{
  if (verbose) { cat("      Integral of the aperture") }
  ssa<-foreach(sam=sa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sam) }
  ssa<-array(unlist(ssa),dim=c(dim(ssa[[1]]),length(ssa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the convolved aperture; ssfa {{{
  if (verbose) { cat("      Integral of the convolved aperture") }
  ssfa<-foreach(sfam=sfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
  ssfa<-array(unlist(ssfa),dim=c(dim(ssfa[[1]]),length(ssfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture)^2]; ssfa2 {{{
  if (verbose) { cat("      Integral of the [(convolved aperture)^2]") }
  ssfa2<-foreach(sfam=sfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum((sfam)^2.)  }
  ssfa2<-array(unlist(ssfa2),dim=c(dim(ssfa2[[1]]),length(ssfa2)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * image); ssfad {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image)") }
  ssfad<-foreach(sfam=sfa, im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      sum(sfam*(im[xlo:xup,ylo:yup]))
  }
  ssfad<-array(unlist(ssfad),dim=c(dim(ssfad[[1]]),length(ssfad)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Quartered Photometry - Quartered Integral of the (convolved aperture * image); qssfad {{{
  if (verbose) { cat("      Quartered Integral of the (convolved aperture * image)") }
  qssfad<-foreach(sfam=sfa, im=im_mask, x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                        y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                       sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                       sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                        .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
      cbind(sum(sfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(sfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(sfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(sfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the reinterpolated psf; spsf {{{
  if (!nopsf) {
    if (verbose) { cat("      Integral of the reinterpolated psf") }
    spsf<-foreach(slen=stamplen, xc=x_g, yc=y_g, .combine='rbind', .options.mpi=mpiopts, .export="psf", .noexport=ls(envir=environment())) %dopar% {
      #for (i in 1:npos) {
      #slen=stamplen[i]
      #Make Sure PSF is centred on centre of stamp
      centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
      delta<-floor(slen/2)*c(-1,+1)
      lims<-rbind(centre[1]+delta,centre[2]+delta)
      dx<-lims[1,1]-1
      dy<-lims[2,1]-1
      #Reinterpolate the PSF at point source XcenYcen {{{
      lenxpsf<-length(psf[,1])
      lenypsf<-length(psf[1,])
      lenx<-length(1:slen)
      leny<-length(1:slen)
      #Make grid for psf at old pixel centres {{{
      psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
      #}}}
      #Make expanded grid of new pixel centres {{{
      expanded<-expand.grid(seq(1,lenx),seq(1,leny))
      xnew<-expanded[,1]-xc%%1
      ynew<-expanded[,2]-yc%%1
      #}}}
      #Interpolate {{{
      ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
      #}}}
      #}}}
      sum(ap)
    }
    #}
    spsf<-array(unlist(spsf),dim=c(dim(spsf[[1]]),length(spsf)))
    if (verbose) { cat(" - Done\n") }
  } else {
    spsf<-rep(NA, length(sfa))
  }
  #}}}
#-----
  #Integral of the (convolved aperture * psf); ssfap {{{
  if (!nopsf) {
    if (verbose) { cat("      Integral of the (convolved aperture * psf)") }
    ssfap<-foreach(sfam=sfa, slen=stamplen, xc=x_g, yc=y_g, .combine='rbind', .options.mpi=mpiopts, .export="psf", .noexport=ls(envir=environment())) %dopar% {
      #for (i in 1:npos) {
      #slen=stamplen[i]
      #Make Sure PSF is centred on centre of stamp
      centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
      delta<-floor(slen/2)*c(-1,+1)
      lims<-rbind(centre[1]+delta,centre[2]+delta)
      dx<-lims[1,1]-1
      dy<-lims[2,1]-1
      #Reinterpolate the PSF at point source XcenYcen {{{
      lenxpsf<-length(psf[,1])
      lenypsf<-length(psf[1,])
      lenx<-length(1:slen)
      leny<-length(1:slen)
      #Make grid for psf at old pixel centres {{{
      psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
      #}}}
      #Make expanded grid of new pixel centres {{{
      expanded<-expand.grid(seq(1,lenx),seq(1,leny))
      xnew<-expanded[,1]-xc%%1
      ynew<-expanded[,2]-yc%%1
      #}}}
      #Interpolate {{{
      ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=leny,nrow=lenx)
      #}}}
      #}}}
      sum(sfam*ap)
    }
    #}
    ssfap<-array(unlist(ssfap),dim=c(dim(ssfap[[1]]),length(ssfap)))
    if (verbose) { cat(" - Done\n") }
  } else {
    ssfap<-rep(NA, length(sfa))
  }
  #}}}
#-----
  #Integral of the deblended convolved aperture; sdfa {{{
  if (verbose) { cat("      Integral of the deblended convolved aperture") }
  sdfa<-foreach(dfam=dfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
  sdfa<-array(unlist(sdfa),dim=c(dim(sdfa[[1]]),length(sdfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(deblended convolved aperture)^2]; sdfa2 {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture)^2]") }
  sdfa2<-foreach(dfam=dfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum((dfam)^2.)  }
  sdfa2<-array(unlist(sdfa2),dim=c(dim(sdfa2[[1]]),length(sdfa2)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (deblended convolved aperture * image); sdfad {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * image)") }
  sdfad<-foreach(dfam=dfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      sum(dfam*(im[xlo:xup,ylo:yup]))
  }
  sdfad<-array(unlist(sdfad),dim=c(dim(sdfad[[1]]),length(sdfad)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Quartered Photometry - Quartered Integral of the (deblended convolved aperture * image); qsdfad {{{
  if (verbose) { cat("      Quartered Integral of the (convolved aperture * image)") }
  qsdfad<-foreach(dfam=dfa, im=im_mask, x1=stamp_lims[,1], x2=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2), x3=stamp_lims[,1]+floor((stamp_lims[,2]-stamp_lims[,1])/2)+1, x4=stamp_lims[,2],
                                        y1=stamp_lims[,3], y2=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2), y3=stamp_lims[,3]+floor((stamp_lims[,4]-stamp_lims[,3])/2)+1, y4=stamp_lims[,4],
                                       sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                       sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                        .inorder=TRUE, .options.mpi=mpiopts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
      cbind(sum(dfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(dfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(dfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(dfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (deblended convolved aperture * convolved aperture); sdfasfa {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * convolved aperture)") }
  sdfasfa<-foreach(dfam=dfa, sfam=sfa, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% { sum(dfam*sfam) }
  sdfasfa<-array(unlist(sdfasfa),dim=c(dim(sdfasfa[[1]]),length(sdfasfa)))
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * image error); ssfae {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime_mask)>1) {
    ssfae<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
      sum(sfam*ime[xlo:xup,ylo:yup])
    }
    ssfae<-array(unlist(ssfae),dim=c(dim(ssfae[[1]]),length(ssfae)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfae<-ssfa
    } else if (ime_mask==0) {
      ssfae<-rep(0,length(ssfa))
    } else {
      ssfae<-foreach(sfam=sfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*ime_mask)
      }
      ssfae<-array(unlist(ssfae),dim=c(dim(ssfae[[1]]),length(ssfae)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the (convolved aperture * (image error)^2); ssfae2 {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
  if (length(ime_mask)>1) {
    ssfae2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfae2<-array(unlist(ssfae2),dim=c(dim(ssfae2[[1]]),length(ssfae2)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfae2<-ssfa
    } else if (ime_mask==0) {
      ssfae2<-rep(0,length(ssfa))
    } else {
      ssfae2<-foreach(sfam=sfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum(sfam*(ime_mask)^2.)
      }
      ssfae2<-array(unlist(ssfae2),dim=c(dim(ssfae2[[1]]),length(ssfae2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture * image error)^2]; ssfa2e2 {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image error)^2]") }
  if (length(ime_mask)>1) {
    ssfa2e2<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((sfam*ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2e2<-array(unlist(ssfa2e2),dim=c(dim(ssfa2e2[[1]]),length(ssfa2e2)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfa2e2<-ssfa2
    } else if (ime_mask==0) {
      ssfa2e2<-rep(0,length(ssfa))
    } else {
      ssfa2e2<-foreach(sfam=sfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum((sfam*ime_mask)^2.)
      }
      ssfa2e2<-array(unlist(ssfa2e2),dim=c(dim(ssfa2e2[[1]]),length(ssfa2e2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [deblended convolved aperture * (image error)^2]; sdfae2 {{{
  if (verbose) { cat("      Integral of the [deblended convolved aperture * (image error)^2]") }
  if (length(ime_mask)>1) {
    sdfae2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
    sdfae2<-array(unlist(sdfae2),dim=c(dim(sdfae2[[1]]),length(sdfae2)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      sdfae2<-sdfa2
    } else if (ime_mask==0) {
      sdfae2<-rep(0,length(ssfa))
    } else {
      sdfae2<-foreach(dfam=dfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum((dfam*ime_mask)^2.)
      }
      sdfae2<-array(unlist(sdfae2),dim=c(dim(sdfae2[[1]]),length(sdfae2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(deblended convolved aperture * image error)^2]; sdfa2e2 {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture * image error)^2]") }
  if (length(ime_mask)>1) {
    sdfa2e2<-foreach(dfam=dfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
    sdfa2e2<-array(unlist(sdfa2e2),dim=c(dim(sdfa2e2[[1]]),length(sdfa2e2)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      sdfa2e2<-sdfa2
    } else if (ime_mask==0) {
      sdfa2e2<-rep(0,length(ssfa))
    } else {
      sdfa2e2<-foreach(dfam=dfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum((dfam*ime_mask)^2.)
      }
      sdfa2e2<-array(unlist(sdfa2e2),dim=c(dim(sdfa2e2[[1]]),length(sdfa2e2)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture * image) / {(image error)^2} ]; ssfadw {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image) / {(image error)^2} ]") }
  if (length(ime_mask)>1) {
    ssfadw<-foreach(sfam=sfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], xelo=estamp_lims[,1],xeup=estamp_lims[,2],yelo=estamp_lims[,3],yeup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup])/((ime[xelo:xeup,yelo:yeup])^2.))
    }
    ssfadw<-array(unlist(ssfadw),dim=c(dim(ssfadw[[1]]),length(ssfadw)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfadw<-ssfad
    } else if (ime_mask==0) {
      ssfadw<-rep(Inf,length(ssfa))
    } else {
      ssfadw<-foreach(sfam=sfa,im=im_mask, xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4], .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum(sfam*(im[xlo:xup,ylo:yup])/((ime_mask)^2.))
      }
      ssfadw<-array(unlist(ssfadw),dim=c(dim(ssfadw[[1]]),length(ssfadw)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [convolved aperture / (image error)^2 ]; ssfaw {{{
  if (verbose) { cat("      Integral of the [convolved aperture / (image error)^2 ]") }
  if (length(ime_mask)>1) {
    ssfaw<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam/(ime[xlo:xup,ylo:yup])^2.)
    }
    ssfaw<-array(unlist(ssfaw),dim=c(dim(ssfaw[[1]]),length(ssfaw)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfaw<-ssfa
    } else if (ime_mask==0) {
      ssfaw<-rep(Inf,length(ssfa))
    } else {
      ssfaw<-foreach(sfam=sfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum(sfam/(ime_mask)^2.)
      }
      ssfaw<-array(unlist(ssfaw),dim=c(dim(ssfaw[[1]]),length(ssfaw)))
    }
  } else {
    sink(type="message")
    stop("Image Error is NULL")
  }
  if (verbose) { cat(" - Done\n") } #}}}
#-----
  #Integral of the [(convolved aperture / image error)^2]; ssfa2w {{{
  if (verbose) { cat("      Integral of the [(convolved aperture / image error)^2]") }
  if (length(ime_mask)>1) {
    ssfa2w<-foreach(sfam=sfa, xlo=estamp_lims[,1],xup=estamp_lims[,2],ylo=estamp_lims[,3],yup=estamp_lims[,4], ime=ime_mask, .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
        sum((sfam/ime[xlo:xup,ylo:yup])^2.)
    }
    ssfa2w<-array(unlist(ssfa2w),dim=c(dim(ssfa2w[[1]]),length(ssfa2w)))
  } else if (length(ime_mask)==1){
    if (ime_mask==1) {
      ssfa2w<-ssfa2
    } else if (ime_mask==0) {
      ssfa2w<-rep(Inf,length(ssfa))
    } else {
      ssfa2w<-foreach(sfam=sfa, .export="ime_mask", .inorder=TRUE, .options.mpi=mpiopts, .noexport=ls(envir=environment())) %dopar% {
          sum((sfam/ime_mask)^2.)
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
  #Calculate Aperture Correction {{{
  if (!nopsf) {
    ApCorr<-spsf/ssfap
  } else {
    ApCorr<-1
  }
  #}}}

  #Convolved aperture flux = Int(fltAp*Im) {{{
  sfaflux<-ssfad
  #sfaflux2<-ssfad*ssfa/ssfa2
  #}}}

  #Deblended convolved aperture flux = Int(DBfltAp*Im) {{{
  dfaflux<-sdfad
  #dfaflux<-sdfad*ssfa/ssfa2
  #}}}

  #Convolved aperture error {{{
  sfaerr<-sqrt((ssfa2e2)*ApCorr^2 + ((conf*beamarea)^2.*sqrt(ssfa)))
  #sfaerr<-sqrt((ssfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(ssfa)))
  #}}}

  #Deblended Convolved aperture error {{{
  dfaerr<-sqrt((sdfa2e2)*ApCorr^2 + ((conf*beamarea)^2.*sqrt(sdfa)))
  #dfaerr<-sqrt((sdfa2e2 * (ssfa/ssfa2)^2.) + ((conf*beamarea)^2.*sqrt(sdfa)))
  #}}}

  #Convolved Aperture Flux Weights & Error Weights {{{
  sfafluxw<-ssfadw*ssfa/ssfa2w
  sfaerrw<-sqrt(1/ssfaw*(ssfa/beamarea*ssfa/ssfa2)^2. + ((conf*beamarea)^2.*sqrt(ssfa)))
  #}}}
  if (verbose) { cat(" - Done\n   } ") }
  #}}}
#-----
  #Notify {{{
  if (showtime) { cat(" - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Perform Calculations - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat(" - Done\n") }
  #}}}
  if (!quiet) { cat("   Performing Final Calculations {  ") }
  if (doskyest|getskyrms) {
    skyflux<-skylocal*sdfa
    skyerr<-skyerr
    if (doskyest) {
      if (!quiet) { message("Perfoming Sky Subtraction"); cat("\n   Performing Sky Subtraction") }
      #Subrtract Sky Flux
      dfaflux<-dfaflux-skyflux
      sfaflux<-sfaflux-skyflux
      dfaerr[which(!is.na(skyerr))]<-sqrt(dfaerr^2+(skyerr*sdfa)^2)[which(!is.na(skyerr))]
      sfaerr[which(!is.na(skyerr))]<-sqrt(sfaerr^2+(skyerr*sdfa)^2)[which(!is.na(skyerr))]
      if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    }
  }
  #}}}
  #Apply Aperture Correction to the Aperture Fluxess {{{
  sfaflux<-sfaflux*ApCorr
  dfaflux<-dfaflux*ApCorr
  #}}}
  #Apply Additional Flux Correction to the Finalised Values {{{
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
    suppressWarnings(mags<--2.5*(log10(dfaflux)-log10(ABvegaflux))+magZP)
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
  #if (length(which(!is.finite(dfaerr[which(sdfa > 0)]))) > 0) {
    #message(paste(length(!is.finite(dfaerr[which(sdfa > 0)])), "elements of dfaerr are not finite"))
    #sink(type="message")
    #stop("NaN or Infs Produced in calculations")
  #}#}}}
  if (!quiet) { cat("} Galaxy Results Complete\n") }
  #}}}
  #PART FIVE: OUTPUT {{{
  #Do we want to plot a sample of the apertures? {{{
  if (plotsample) {
    #Set output name {{{
    dir.create(file.path(pathroot,pathwork,pathout,"COGs"),showWarnings=FALSE)
    #}}}
    #Determine Sample to Plot {{{
    if (!exists("plotall")) { plotall<-FALSE }
    if (!plotall) {
      if (!quiet) { message("Writing Sample of COGs to File"); cat("Writing Sample of COGs to File") }
      #Output a random 15 apertures to file {{{
      ind1=which(a_g>0)
      ind1=ind1[order(sdfad[ind1],decreasing=T)][1:15]
      ind2=order(runif(1:npos))[1:15]
      ind<-c(ind1, ind2)
      rm(ind1)
      rm(ind2)
      ind<-ind[which(!is.na(ind))]
      #}}}
    } else {
      #All {{{
      if (!quiet) { message("Writing All COGs to File"); cat("Writing All COGs to File") }
      ind=c(which(a_g>0),which(a_g==0))
      #}}}
    }
    #}}}
    for (i in ind) {
      #Open Device {{{
      png(file.path(pathroot,pathwork,pathout,paste("COGs/",id_g[i],".png",sep="")),width=14*240,height=3.5*240,res=240)
      #}}}
      #Set Layout {{{
      layout(cbind(1,2,3,4))
      #}}}
      #Axes Limits {{{
      xlims=round(max(abs(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, na.rm=TRUE))*c(-1,1)
      #}}}
      #Make Ap and Source Mask Block/Trans matricies {{{
      apT<-sfa[[i]]
      apT[which(sfa[[i]]==0)]<-NA
      apT[which(sfa[[i]]!=0)]<-1
      apB<-sfa[[i]]
      apB[which(sfa[[i]]==0)]<-1
      apB[which(sfa[[i]]!=0)]<-NA
      smB<-sm[mask_lims[i,1]:mask_lims[i,2],mask_lims[i,3]:mask_lims[i,4]]
      smB[which(smB==0)]<-NA
      #}}}
      #Plot Aperture in Blue {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=(sfa[[i]]), main="Image & Aperture", asp=1, col=hsv(2/3,seq(0,1,length=256)), useRaster=FALSE, axes=FALSE, xlab="", ylab="", xlim=xlims, ylim=xlims))
      #}}}
      #Plot Image in greyscale {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=FALSE,add=TRUE, xlab="", ylab=""))
      #}}}
      #Overlay Sourcemask in Green {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(smB*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]), main="", asp=1, useRaster=FALSE,add=TRUE, xlab="", ylab="",col=cm.colors(256)))
      #}}}
      #Plot +ve flux in aperture in Heat Colours {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(apT*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]), main="", asp=1, col=heat.colors(256), useRaster=FALSE,add=TRUE, xlab="", ylab=""))
      #}}}
      #Overlay Aperture in Black {{{
      #suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=log10(sfa[[i]]), main="Image & Aperture", asp=1, useRaster=FALSE, axes=FALSE, xlab="", ylab="", xlim=xlims, ylim=xlims,add=TRUE))
      #}}}
      #Plot Sources {{{
      points(x=(x_p-x_p[i]+1)*asperpix,y=(y_p-y_p[i]+1)*asperpix, pch=3)
      #}}}
      #Label with ID {{{
      label("topleft",lab=id_g[i],cex=1.5, col='red')
      #}}}
      #Draw Axes {{{
      magaxis(frame.plot=T,main="Image & Aperture",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      #}}}
      #Generate COGs {{{
      #Sky Subtracted only {{{
      cog<-get.cog(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      #}}}
      #Raw Image {{{
      cog_nosky<-get.cog(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i],centre=c(stamplen[i]/2, stamplen[i]/2))
      #}}}
      #Deblended only {{{
      if (PSFWeighted) {
        debl.cog<-get.cog(sfa[[i]]*dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        debl.cog<-get.cog(dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2))
      }
      #}}}
      #Deblended and Sky Subtracted {{{
      if (PSFWeighted) {
        debl.cog_nosky<-get.cog(sfa[[i]]*dbw[[i]]*(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2))
      } else {
        debl.cog_nosky<-get.cog(dbw[[i]]*(image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2))
      }
      #}}}
      #}}}
      #Plot COGs {{{
      if (Magnitudes) {
        #Plot in Magnitude Space {{{
        #Get y limits {{{
        suppressWarnings(ylim<-c((-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP)+3, (-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP)-3))
        #If a limit in ±Inf, plot around median
        if (!all(is.finite(ylim))) { ylim=median(cog$y, na.rm=TRUE)+c(-1,3) }
        #}}}
        #Plot Raw COG {{{
        suppressWarnings(magplot(x=cog$x*asperpix, y=-2.5*(log10(cog$y)-log10(ABvegaflux))+magZP, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Magnitude",
                ylim=ylim,type='l',lty=2,main="Curve of Growth"))
        #}}}
        #Add Lines showing fluxes {{{
        #Undeblended
        suppressWarnings(abline(h=-2.5*(log10(sfaflux[i])-log10(ABvegaflux))+magZP, lwd=1, col='orange', lty=1))
        #Deblended
        suppressWarnings(abline(h=-2.5*(log10(dfaflux[i])-log10(ABvegaflux))+magZP, lwd=1, col='green', lty=1))
        #}}}
        #Redraw Raw Cog {{{
        suppressWarnings(lines(x=cog$x*asperpix, y=-2.5*(log10(cog$y)-log10(ABvegaflux))+magZP, pch=20, col='grey',lty=2))
        #}}}
        #Draw Deblended Cog {{{
        suppressWarnings(lines(x=debl.cog$x*asperpix, y=-2.5*(log10(debl.cog$y)-log10(ABvegaflux))+magZP, pch=20, col='black',lty=2))
        #}}}
        #Draw Sky Subtracted Cog {{{
        suppressWarnings(lines(x=cog_nosky$x*asperpix, y=-2.5*(log10(cog_nosky$y)-log10(ABvegaflux))+magZP, pch=20, col='grey',lty=1))
        #}}}
        #Draw Sky Subtracted & Deblended Cog {{{
        suppressWarnings(lines(x=debl.cog_nosky$x*asperpix, y=-2.5*(log10(debl.cog_nosky$y)-log10(ABvegaflux))+magZP, pch=20, col='black',lty=1))
        #}}}
        #Draw Legend {{{
        legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended ApMag","Deblended ApMag"),lty=c(2,2,1,1,1,1),
               col=c('grey','black','grey','black','orange','green'),pch=-1, cex=0.5)
        #}}}
        #}}}
      } else {
        #Plot in Flux Space {{{
        #Plot Raw Cog {{{
        magplot(x=cog$x*asperpix, y=cog$y, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Flux",ylim=rev(range(cog$y)),main="Curve of Growth")
        #}}}
        #Draw Lines showing fluxes {{{
        abline(h=dfaflux[i], lwd=1, col='green')
        abline(h=sfaflux[i], lwd=1, col='orange', lty=2)
        #}}}
        #Redraw Raw Cog {{{
        lines(x=cog$x*asperpix, y=cog$y, pch=20, col='grey',lty=2)
        #}}}
        #Draw Deblended Cog {{{
        lines(x=debl.cog$x*asperpix, y=debl.cog$y, pch=20, col='black',lty=2)
        #}}}
        #Draw Sky Subtracted Cog {{{
        lines(x=cog_nosky$x*asperpix, y=cog_nosky$y, pch=20, col='grey',lty=1)
        #}}}
        #Draw Sky Subtracted & Deblended Cog {{{
        lines(x=debl.cog_nosky$x*asperpix, y=debl.cog_nosky$y, pch=20, col='black',lty=1)
        #}}}
        legend('topleft',legend=c("Undeblended COG","Deblended COG","Undeblended ApFlux","Undeblended ApFlux - Sky","Deblended ApFlux","Deblended ApFlux - Sky"),lty=c(-1,-1,2,1,2,1),
               col=c('grey','black',col2alpha('orange',0.5),'orange',col2alpha('green',0.5),'green'),pch=c(20,20,-1,-1,-1,-1,-1))
        #}}}
      }
      #}}}
      #Plot the Deblended Image {{{
      suppressWarnings(z<-log10(dbw[[i]]*image.env$im[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]))
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apB, main="Image x Weight Matrix", asp=1, col=grey.colors(256), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(-4,4), xlim=xlims, ylim=xlims)
      #}}}
      #Overlay the Aperture {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apT, main="Image x Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(-4,4), xlim=xlims, ylim=xlims,add=TRUE))
      #}}}
      #Draw the Axes and scalebar {{{
      magaxis(frame.plot=T,main="Image x Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      magbar("topright",col=rev(rainbow(256, start=0,end=2/3)), range=c(-4,4),title='log(pixval)')
      #}}}
      #Plot the Deblend Matrix {{{
      z=dbw[[i]]
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apB, main="Weight Matrix", asp=1, col=grey.colors(256), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=xlims)
      #}}}
      #Overlay the Aperture {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*asperpix, z=z*apT, main="Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=FALSE, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=xlims,add=T))
      #}}}
      #Draw the Axes and scalebar {{{
      magaxis(frame.plot=T,main="Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      magbar("topright",col=rev(rainbow(256, start=0,end=2/3)), range=c(0,1))
      #}}}
      #Close the file {{{
      dev.off()
      #}}}
    }
    #Remove unneeded Arrays {{{
    if (!makeresidmap) {
      sfabak<-NULL
    }
    rm(dbw)
    #}}}
    #Notify {{{
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    #}}}
  }
  #}}}
  #If map was input in Jy/bm we need to convert it back before output in SourceSubtraction {{{
  if (Jybm) { ba=beamarea } else { ba=1. }
  #}}}
  #If wanted, make the Residual Map {{{
  if (makeresidmap) {
    if (!is.null(sfabak)) { sfa<-sfabak }
    if (filtcontam) {
      if (!quiet) { cat(paste("Writing Contaminant-subtracted Map to",nocontammap,"   ")) }
      #Perform Source Subtraction
      timer=system.time(sourcesubtraction(im,sfa,image_lims,dfaflux/ApCorr,file.path(pathroot,pathwork,pathout,nocontammap),hdr_str,ba,contams,diagnostic,verbose))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Contam Subtraction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
    if (!quiet) { cat(paste("Writing Source-subtracted Map to",residmap,"   ")) }
    #Perform Source Subtraction
    timer=system.time(sourcesubtraction(im,sfa,image_lims,dfaflux/ApCorr,file.path(pathroot,pathwork,pathout,residmap),hdr_str,ba,insidemask,diagnostic,verbose))
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
    sdfae2<-sdfae2[which(contams==0)]
    sdfa2e2<-sdfa2e2[which(contams==0)]
    dfaflux<-dfaflux[which(contams==0)]
    dfaerr <-dfaerr[which(contams==0)]
    pixflux<-pixflux[which(contams==0)]
    ApCorr <-ApCorr[which(contams==0)]
    stamplen<-stamplen[which(contams==0)]
    mags   <-mags[which(contams==0)]
    if (length(fluxweight!=1)) { fluxweight<-fluxweight[which(contams==0)] }
    if (RanCor) { randoms<-randoms[which(contams==0),] }
  }#}}}
  #If wanted, output the Results Table {{{
  if (writetab) {
    #Output the results table
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
  #PART SIX: FINISH {{{
  #Send Parameters to logfile {{{
  sink(sinkfile, type="output")
  cat("Memory Hogs in this run:\n")
  print(lsos(envir=environment(), head=TRUE, n=10))
  cat("Images used in this run:\n")
  print(lsos(envir=image.env, head=FALSE))
  sink(type="output")
  #}}}
  #Return {{{
  if (!quiet) { cat('\n') }
  on.exit(detach(image.env))
  if (!is.null(env)) {
    on.exit(detach(env), add=TRUE)
  }
  if (!diagnostic) {
    return=list(SFAflux=sfaflux,SFAerror=sfaerr, DFAflux=dfaflux,DFAerror=dfaerr)
  } else {
    return=list(SFAflux=sfaflux,SFAerror=sfaerr, DFAflux=dfaflux,DFAerror=dfaerr, SA_Stamps=sa, SFA_Stamps=sfa, WSFA_Stamps=wsfa, DFA_Stamps=dfa)
  }
  #}}}
  #}}}

}
#Func {{{
.executeRanCor<-function() {
cat("Executing RanCor...\n"); Sys.sleep(2); cat('\n    |   _______   _     _   _     _   _ __      _    |\n    |  |__   __| | |   | | | |   | | |  __ \\   | |   |\n    |     | |    | |___| | | |   | | | |  | |  | |   |\n    |     | |    |  ___  | | |   | | | |  | |  |_|   |\n    |     | |    | |   | | | |___| | | |__| |   _    |\n    |     |_|    |_|   |_| |_______| |_____/   |_|   |\n    |                                                |\n     \\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\n      \\/  \\/  \\/  \\/  ,------------,  \\/  \\/  \\/  \\/\n         ____        ({ XX      XX })        ____\n        /////|\\      _|   \\\\  //   |_      /|\\\\\\\\\\\n        VVVV | \\____/ { \\  \\\\//  / } \\____/ | VVVV\n        ////_|       ,|V=V=V=V=V=V=|,       |_\\\\\\\\\n     ___\\\\\\\\/_\\~~~~~/_{+^+^+^+^+^+^}_\\~~~~~/_\\////___\n\n'); Sys.sleep(2); cat("... Done, you heartless Jedi Scum\n\n"); Sys.sleep(2)
}
#}}}
