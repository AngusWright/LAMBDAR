flux.measurements <-
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
  environment(make.catalogue.apertures)<-environment()
  environment(make.convolved.apertures)<-environment()
  environment(make.aperture.map)<-environment()
  environment(make.data.array.maps)<-environment()
  environment(read.psf)<-environment()
  environment(make.deblend.weight.maps)<-environment()
  environment(write.flux.measurements.table)<-environment()
  environment(write.deblend.fraction.table)<-environment()
  environment(write.fits.image.file)<-environment()
  environment(source.subtraction)<-environment()
  # /*fend*/ }}}

  if (!quiet) { cat('Begining Flux Measurements\n') }
  message('} Initialisation Complete\nBegining Flux Measurements')
  dimim<-dim(image.env$im)

  #Get PSF details /*fold*/ {{{
  timer<-proc.time()
  if (!quiet) { cat('Getting PSF details') }
  message('Getting PSF details')
  if (!(no.psf)) {
    #There is a PSF specified /*fold*/ {{{
    #Details /*fold*/ {{{
    #Calculate PSF - if one has not be supplied,
    #then a gaussian PSF will be created. Stampsize of PSF
    #should be = maximum stampsize: max aperture * stamp mult /*fend*/ }}}
    #Get PSF /*fold*/ {{{
    psf<-read.psf(outenv=environment(),file.path(path.root,path.work,psf.map),arcsec.per.pix,max(cat.a),confidence,gauss.fwhm.arcsec=gauss.fwhm.arcsec)
    # /*fend*/ }}}
    #Notify /*fold*/ {{{
    if (verbose) { message(paste("Maxima of the PSF is at pixel", which(psf == max(psf)),"and has value",max(psf))) }
    # /*fend*/ }}}
    #Get radius of FWHM using FWHM confidence value = erf(2*sqrt(2*log(2))/sqrt(2)) /*fold*/ {{{
    psffwhm<-get.fwhm(psf)
    # /*fend*/ }}}
    #Normalise Beam Area /*fold*/ {{{
    beam.area.nn<-sumpsf
    beam.area.n<-as.single(sum(psf))
    # /*fend*/ }}}
    #-----Diagnostic-----## /*fold*/ {{{
    if (diagnostic) { message(paste('Beam area before/after norm: ',beam.area.nn,beam.area.n)) }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else {
    #There is no PSF specified /*fold*/ {{{
    #set relevant paramters manually
    beam.area.n<-1.0
    psfwidth<-0
    psf.clip<-0
    sumpsf<-0
    # /*fend*/ }}}
    #If no PSF, set FWHM to median aperture radius + buffer /*fold*/ {{{
    psffwhm<-round(min(cat.a[which(cat.a>0)])/arcsec.per.pix)
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
  if (beam.area.pix == 0) {
    #Use the beamarea just determined /*fold*/ {{{
    beamarea<-beam.area.n
    # /*fend*/ }}}
  } else {
    # Otherwise just use the input beamarea /*fold*/ {{{
    beamarea<-beam.area.pix
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
  if (psf.filt) {
    stamplen<-(floor(ceiling(def.buff*cat.a/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
  } else {
    stamplen<-(floor(ceiling(def.buff*cat.a/arcsec.per.pix))*2+5)
  }
  # /*fend*/ }}}
  #Discard any contaminants that are beyond their nearest neighbours stamp /*fold*/ {{{
  #Number of Nearest Neighbors to search /*fold*/ {{{
  if (!exists("num.nearest.neighbours")) { num.nearest.neighbours<-10 }
  if (!exists("check.contam")) { check.contam<-TRUE }
  if (!exists("getLoners")) { getLoners<-FALSE }
  # /*fend*/ }}}
  if (getLoners) {
    timer<-proc.time()
    if (!quiet) { cat("Selecting out only objects that are completely alone ") }
    cat.len<-length(cat.x)
    num.nearest.neighbours<-max(num.nearest.neighbours,cat.len-1)
    nearest<-nn2(data.frame(cat.x,cat.y),data.frame(cat.x,cat.y),k=num.nearest.neighbours+1)
    if (psffwhm!=0) {
      #If any of the nearest neighbors overlap with the object /*fold*/ {{{
      inside.mask<-rowSums(nearest$nn.dist[,2:num.nearest.neighbours] < 4*psffwhm)==0
    } else {
      inside.mask<-rep(NA,cat.len)
      for(ind in 1:cat.len) {
        #If any of the nearest neighbors overlap with the object /*fold*/ {{{
        inside.mask[ind]<-any(nearest$nn.dist[ind,1:num.nearest.neighbours+1] < (sqrt(2)/2*(stamplen[nearest$nn.idx[ind,1:num.nearest.neighbours+1]]+stamplen[ind])))
        # /*fend*/ }}}
      }
      # /*fend*/ }}}
    }
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
    if (exists("contams")) { contams<-contams[which(inside.mask)] }
    chunk.size=length(cat.id)/getDoParWorkers()
    mpi.opts<-list(chunkSize=chunk.size)
    message("Number of objects per thread:",chunk.size)
    # /*fend*/ }}}
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects that are entirely unblended (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of objects were objects that were possibly blended)")) }
    # /*fend*/ }}}
    if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      message(paste('Remove irrelevant contaminants - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  } else if (check.contam && exists("contams") && length(which(contams==0))>1 && length(which(contams==1))>1) {
    timer<-proc.time()
    if (!quiet) { cat("Removing Contaminants that are irrelevant ") }
    cat.len<-length(cat.x)
    num.nearest.neighbours<-min(num.nearest.neighbours,length(which(contams==0))-1)
    nearest<-nn2(data.frame(cat.x[which(contams==0)],cat.y[which(contams==0)]),data.frame(cat.x[which(contams==1)],cat.y[which(contams==1)]),k=num.nearest.neighbours)
    contam.inside<-NULL
    for(ind in 1:length(which(contams==1))) {
      #If any of the nearest neighbors overlap with the contaminant /*fold*/ {{{
      contam.inside<-c(contam.inside,any(nearest$nn.dist[ind,] < (sqrt(2)/2*(stamplen[which(contams==0)][nearest$nn.idx[ind,]]+stamplen[which(contams==1)][ind]))))
    }
    inside.mask<-rep(TRUE,cat.len)
    inside.mask[which(contams==1)]<-contam.inside
    # /*fend*/ }}}
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
    contams<-contams[which(inside.mask)]
    chunk.size=length(cat.id)/getDoParWorkers()
    mpi.opts<-list(chunkSize=chunk.size)
    message("Number of objects per thread:",chunk.size)
    # /*fend*/ }}}
    #Notify how many objects remain /*fold*/ {{{
    if (verbose) { message(paste("There are",length(cat.x),"supplied objects & contaminants that intersect (",
                                  round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of objects were contaminants that didn't intersect galaxies)")) }
    # /*fend*/ }}}
    if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
      message(paste('Remove irrelevant contaminants - Done (',round(proc.time()[3]-timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }
  # /*fend*/ }}}
  #Convert decimal pixel values into actual pixel values /*fold*/ {{{
  x.pix<-floor(cat.x)
  y.pix<-floor(cat.y)
  # /*fend*/ }}}

  #Create an array of stamps containing the data image subsection for all objects /*fold*/ {{{
  timer=system.time(data.stamp<-make.data.array.maps(outenv=environment()))
  if (((loop.total<=1)&(cutup)&(length(data.stamp))==0)) { sink(type="message") ; stop("No Aperture Stamps are aligned enough to be valid.") }
  else if ((cutup)&(length(data.stamp)==0)) {
      if (showtime) { cat("   - BREAK (",round(timer[3],digits=2),"sec )\n")
        message(paste('Make Data Mask - BREAK (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - BREAK\n") }
      warning("No Aperture Stamps are aligned enough to be valid.")
      #Notify & Close Logfile /*fold*/ {{{
      if (!is.null(env)) {
        on.exit(detach(env), add=TRUE)
      }
      message(paste('\n-----------------------------------------------------\nDatamap Skipped - No Aperture Stamps are aligned enough to be valid\n'))
      if (!quiet) {
        cat(paste('\n-----------------------------------------------------\nDatamap Skipped - No Apertures Stamps are aligned enough to be valid'))
      }
      return(NULL)
      # /*fend*/ }}}
  }
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Data Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Remove Image arrays as they are no longer needed /*fold*/ {{{
  #if (cutup) {
  #  imenvlist<-ls(envir=image.env)
  #  imenvlist<-imenvlist[which(imenvlist!="im"&imenvlist!="data.hdr")]
  #  rm(list=imenvlist, envir=image.env)
  #}
  # /*fend*/ }}}
  #Create an array of stamps containing the apertures for all objects /*fold*/ {{{
  timer=system.time(sa<-make.catalogue.apertures(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make SA Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Discard any apertures that were zero'd in the make process /*fold*/ {{{
  totsa<-foreach(sam=sa, i=1:length(sa), .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) } ; sum(sam) }
  inside.mask<-(totsa > 0)
  if ((loop.total<=1)&(length(which(inside.mask==TRUE))==0)) { sink(type="message") ; stop("No Apertures are inside the Mask after Aperture Creation.") }  # Nothing inside the mask
  else if (length(which(inside.mask==TRUE))==0) {
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
  sa<-sa[which(inside.mask)]
  if (length(data.stamp )>1) { data.stamp<-data.stamp[which(inside.mask)] }
  if (length(mask.stamp)>1) { mask.stamp<-mask.stamp[which(inside.mask)] }
  if (length(error.stamp)>1) { error.stamp<-error.stamp[which(inside.mask)] }
  stamplen<-stamplen[which(inside.mask)]
  ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
  ap.lims.mask.stamp<-rbind(ap.lims.mask.stamp[which(inside.mask),])
  ap.lims.error.stamp<-rbind(ap.lims.error.stamp[which(inside.mask),])
  data.stamp.lims<-rbind(data.stamp.lims[which(inside.mask),])
  mask.stamp.lims<-rbind(mask.stamp.lims[which(inside.mask),])
  error.stamp.lims<-rbind(error.stamp.lims[which(inside.mask),])
  ap.lims.data.map<-rbind(ap.lims.data.map[which(inside.mask),])
  ap.lims.mask.map<-rbind(ap.lims.mask.map[which(inside.mask),])
  x.pix<-x.pix[which(inside.mask)]
  y.pix<-y.pix[which(inside.mask)]
  cat.x<-cat.x[which(inside.mask)]
  cat.y<-cat.y[which(inside.mask)]
  cat.id<-cat.id[which(inside.mask)]
  cat.ra<-cat.ra[which(inside.mask)]
  cat.dec<-cat.dec[which(inside.mask)]
  cat.theta<-cat.theta[which(inside.mask)]
  theta.offset<-theta.offset[which(inside.mask)]
  cat.a<-cat.a[which(inside.mask)]
  cat.a.pix<-cat.a.pix[which(inside.mask)]
  cat.b<-cat.b[which(inside.mask)]
  if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
  if (exists("contams")) { contams<-contams[which(inside.mask)] }
  inside.mask<-inside.mask[which(inside.mask)]
  chunk.size=ceiling(length(cat.id)/getDoParWorkers())
  mpi.opts<-list(chunkSize=chunk.size)
  message("Number of objects per thread:",chunk.size)
  # /*fend*/ }}}
  #Re-Initialise object count /*fold*/ {{{
  npos<-length(cat.id)
  # /*fend*/ }}}
  #Create a full mask of all apertures in their correct image-space locations /*fold*/ {{{
  timer=system.time(image.env$aa<-make.aperture.map(outenv=environment(), sa, dimim))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
   message(paste('Make A Mask - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #If wanted, output the All Apertures Mask /*fold*/ {{{
  if (make.all.apertures.map) {
    if (!quiet) { cat(paste('Outputting All Apertures Mask to',all.apertures.map.filename,"   ")) }
    #Write All Apertures Mask to file /*fold*/ {{{
    timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,all.apertures.map.filename),image.env$aa,image.env$data.hdr,nochange=TRUE))
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
  #If wanted, use pixel fluxes as flux.weights /*fold*/ {{{
  if (use.pixel.fluxweight) {
    #Image Flux at central pixel /*fold*/ {{{
    cat("Determine Image Flux at central pixel ")
    if (cutup) {
      pixflux<-foreach(xp=x.pix-(data.stamp.lims[,1]-1),yp=y.pix-(data.stamp.lims[,3]-1), im=data.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
            im[xp,yp]
      }
    } else {
      pixflux<-foreach(xp=x.pix,yp=y.pix, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment()), .export=c('image.env')) %dopar% {
            image.env$im[xp,yp]
      }
    }
    if (verbose) { cat(" - Done\n") }
    # /*fend*/ }}}
    #Determine Noise Characteristics /*fold*/ {{{
    if (verbose) { cat("Determine Image Noise Characteristics ") }
    if (any(pixflux==-99)) { warning("Some pixel flux determinations failed"); pixflux[which(pixflux==-99)]<-NA }
    quan<-quantile(image.env$im,c(0.1,0.9),na.rm=TRUE)
    bw<-abs(quan[2]-quan[1])/100/sqrt(12)
    if (bw<=0) {
      bw<-diff(range(image.env$im))/100/sqrt(12)
    }
    mode<-density(as.numeric(image.env$im),from=quan[1],to=quan[2],kernel='rect',bw=bw,na.rm=TRUE)
    mode<-mode$x[which.max(mode$y)]
    mad<-median(abs(image.env$im-mode),na.rm=TRUE)
    if (length(which(pixflux>mode+mad))==0) {
      message("WARNING: No objects have pixel flux measurements that are > Pixel Mode+MAD. Pixel Flux weighting cannot be used\n")
        flux.weight<-1
    } else {
      flux.weight<-magmap(pixflux, lo=mode+mad, hi=max(pixflux), range=c(0.01,1), type="num", stretch='lin',bad=0.01)$map
    }
    cat(" - Done\n")
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #If needed, do Convolutions & Fluxweightings /*fold*/ {{{
  if ((!psf.filt)&(length(which(flux.weight!=1))!=0)) {
    #No Convolution, Need Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If not convolving with psf, and if all the flux.weights are not unity,
    #then skip filtering, duplicate the arrays, and only weight the stamps/apertures
    # /*fend*/ }}}
    #No Convolution, duplicate apertures /*fold*/ {{{
    if (verbose) { message("No Convolution: Convolved Apertures are identical to Simple Apertures") }
    sfa<-sa
    image.env$fa<-image.env$aa
    # /*fend*/ }}}
    #Perform aperture weighting /*fold*/ {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make.convolved.apertures(outenv=environment(), sa,flux.weightin=flux.weight,immask=mask.stamp))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$wfa<-make.aperture.map(outenv=environment(), wsfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if ((!psf.filt)&(length(which(flux.weight!=1))==0)) {
    #No Convolution, No Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If not convolving with psf, and if all the flux.weights are unity,
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
  } else if ((psf.filt)&(length(which(flux.weight!=1))!=0)) {
    #Convolving PSF, Need Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If convolving with psf, and if all the flux.weights are not unity,
    #then colvolve and weight the stamps/apertures /*fend*/ }}}
    #Perform Convolution /*fold*/ {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make.convolved.apertures(outenv=environment(), sa,immask=mask.stamp))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all convolved stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$fa<-make.aperture.map(outenv=environment(), sfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Perform aperture weighting /*fold*/ {{{
    if (verbose) { message("Fluxweights Present: Weighting Convolved Apertures") }
    timer=system.time(wsfa<-make.convolved.apertures(outenv=environment(), sa,flux.weightin=flux.weight,immask=mask.stamp))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WSFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all the convolved weighted stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$wfa<-make.aperture.map(outenv=environment(), wsfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make WFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    # /*fend*/ }}}
  } else if ((psf.filt)&(length(which(flux.weight!=1))==0)) {
    #Convolve PSF, No Fluxweighting /*fold*/ {{{
    #Details /*fold*/ {{{
    #If convolving with psf, and if all the flux.weights are unity,
    #colvolve, and then duplicate the stamps/apertures /*fend*/ }}}
    #Perform Convolution /*fold*/ {{{
    if (verbose) { message("PSF Present: Making Convolved Apertures") }
    timer=system.time(sfa<-make.convolved.apertures(outenv=environment(), sa,immask=mask.stamp))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make SFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #Create a full mask of all convolved stamps/apertures /*fold*/ {{{
    timer=system.time(image.env$fa<-make.aperture.map(outenv=environment(), sfa, dimim))
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
  totsfa<-foreach(sam=sfa, i=1:length(sfa), .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) } ; sum(sam) }
  inside.mask<-(totsfa > 0)
  if ((loop.total<=1)&(length(which(inside.mask==TRUE))==0)) { sink(type="message") ; stop("No Apertures Remain within the mask after convolution.") }  # Nothing inside the mask
  else if (length(which(inside.mask==TRUE))==0) {
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
  sa<-sa[which(inside.mask)]
  sfa<-sfa[which(inside.mask)]
  wsfa<-wsfa[which(inside.mask)]
  if (length(data.stamp )>1) { data.stamp<-data.stamp[which(inside.mask)] }
  if (length(mask.stamp)>1) { mask.stamp<-mask.stamp[which(inside.mask)] }
  if (length(error.stamp)>1) { error.stamp<-error.stamp[which(inside.mask)] }
  stamplen<-stamplen[which(inside.mask)]
  ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
  ap.lims.mask.stamp<-rbind(ap.lims.mask.stamp[which(inside.mask),])
  ap.lims.error.stamp<-rbind(ap.lims.error.stamp[which(inside.mask),])
  data.stamp.lims<-rbind(data.stamp.lims[which(inside.mask),])
  mask.stamp.lims<-rbind(mask.stamp.lims[which(inside.mask),])
  error.stamp.lims<-rbind(error.stamp.lims[which(inside.mask),])
  ap.lims.data.map<-rbind(ap.lims.data.map[which(inside.mask),])
  ap.lims.mask.map<-rbind(ap.lims.mask.map[which(inside.mask),])
  x.pix<-x.pix[which(inside.mask)]
  y.pix<-y.pix[which(inside.mask)]
  cat.x<-cat.x[which(inside.mask)]
  cat.y<-cat.y[which(inside.mask)]
  cat.id<-cat.id[which(inside.mask)]
  cat.ra<-cat.ra[which(inside.mask)]
  cat.dec<-cat.dec[which(inside.mask)]
  cat.theta<-cat.theta[which(inside.mask)]
  theta.offset<-theta.offset[which(inside.mask)]
  cat.a<-cat.a[which(inside.mask)]
  cat.a.pix<-cat.a.pix[which(inside.mask)]
  cat.b<-cat.b[which(inside.mask)]
  if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
  if (exists("contams")) { contams<-contams[which(inside.mask)] }
  inside.mask<-inside.mask[which(inside.mask)]
  chunk.size=ceiling(length(cat.id)/getDoParWorkers())
  mpi.opts<-list(chunkSize=chunk.size)
  message("Number of objects per thread:",chunk.size)
  # /*fend*/ }}}
  # /*fend*/ }}}
  #Re-Initialise object count /*fold*/ {{{
  npos<-length(cat.id)
  # /*fend*/ }}}
  #Remove Uneeded Image /*fold*/ {{{
  rm(aa, envir=image.env)
  gc()
  # /*fend*/ }}}
  #-----Diagnostic-----# /*fold*/ {{{
  if (diagnostic) {
    if (verbose) { message("Checking Apertures for scaling errors") }
    foreach(sam=sa, sfam=sfa, i=1:length(sa), .combine=function(a,b){NULL},.inorder=FALSE, .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
      if (max(sam)>1){warning(paste("Max of Aperture",i,"is",max(sam))) }
      if (max(sfam)>1){warning(paste("Max of Convolved Aperture",i,"is",max(sfam))) }
    }
  }# /*fend*/ }}}
  #Do we want to plot a sample of the apertures? /*fold*/ {{{
  if (plot.sample) {
    res<-120
    #Set output name /*fold*/ {{{
    CairoPNG(file.path(path.root,path.work,path.out,"PSF_Samples.png"))
    # /*fend*/ }}}
    #Set Layout /*fold*/ {{{
    par(mfrow=c(2,2))
    # /*fend*/ }}}
    #Output a random 15 apertures to file /*fold*/ {{{
    ind1=which(cat.a>0)
    ind1=ind1[order(runif(1:length(ind1)))][1:15]
    ind2=order(runif(1:npos))[1:15]
    ind<-c(ind1, ind2)
    rm(ind1)
    rm(ind2)
    ind<-ind[which(!is.na(ind))]
    for (i in ind) {
      Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
      image(sa[[i]], main="Single Aperture (SA)", asp=1, col=heat.colors(256), useRaster=Rast)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
      image(sfa[[i]], main="Single Convolved Apertrure (SFA)", asp=1, col=heat.colors(256), useRaster=Rast)
      points(ceiling(stamplen[i]/2)/stamplen[i],ceiling(stamplen[i]/2)/stamplen[i],pch="+",lw=2.0, col="red")
    }# /*fend*/ }}}
    #Close the file /*fold*/ {{{
    dev.off()
    # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #If wanted, output the Convolved & Weighted Aperture Mask /*fold*/ {{{
  if (make.convolved.apertures.map) {
    if (!quiet) { cat(paste('Outputting All Convolved Apertures Mask to',fa.filename,"   ")) }
    timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,fa.filename),image.env$wfa,image.env$data.hdr,nochange=TRUE) )
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
    if (!exists("transmission.map")) { transmission.map<-FALSE }
    if (!quiet) { cat(paste("Creating Sourcemask    ")) }
    if (length(image.env$imm)>1) { sm<-image.env$imm } else { sm<-array(1, dim=dimim) }
    #sm<-array(1, dim=dimim)
    #Get Mask as is {{{
    if (length(mask.stamp)>1) {
      for (i in 1:npos) {
        sm[mask.stamp.lims[i,1]:mask.stamp.lims[i,2],mask.stamp.lims[i,3]:mask.stamp.lims[i,4]]<-mask.stamp[[i]]
        #sm[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]<-mask.stamp[[i]]
      }
    } #}}}
    mask.xrange<-c(min(mask.stamp.lims[,1]),max(mask.stamp.lims[,2]))
    mask.yrange<-c(min(mask.stamp.lims[,3]),max(mask.stamp.lims[,4]))
    data.xrange<-c(min(data.stamp.lims[,1]),max(data.stamp.lims[,2]))
    data.yrange<-c(min(data.stamp.lims[,3]),max(data.stamp.lims[,4]))
    if (!psf.filt) {
      if (transmission.map) {
        sm[mask.xrange[1]:mask.xrange[2],mask.yrange[1]:mask.yrange[2]][which(!image.env$fa[data.xrange[1]:data.xrange[2],data.yrange[1]:data.yrange[2]] > 0)]<-0
        #sm[which(!image.env$fa > 0)]<-0
      } else {
        sm[mask.xrange[1]:mask.xrange[2],mask.yrange[1]:mask.yrange[2]][which(image.env$fa[data.xrange[1]:data.xrange[2],data.yrange[1]:data.yrange[2]] > 0)]<-0
        #sm[which(image.env$fa > 0)]<-0
      }
    } else {
      if (transmission.map) {
        sm[mask.xrange[1]:mask.xrange[2],mask.yrange[1]:mask.yrange[2]][which(!image.env$fa[data.xrange[1]:data.xrange[2],data.yrange[1]:data.yrange[2]] > max(psf, na.rm=TRUE)*(1-sourcemask.conf.lim))]<-0
        #sm[which(!image.env$fa > max(psf, na.rm=TRUE)*(1-sourcemask.conf.lim))]<-0
      } else {
        sm[mask.xrange[1]:mask.xrange[2],mask.yrange[1]:mask.yrange[2]][which(image.env$fa[data.xrange[1]:data.xrange[2],data.yrange[1]:data.yrange[2]] > max(psf, na.rm=TRUE)*(1-sourcemask.conf.lim))]<-0
        #sm[which(image.env$fa > max(psf, na.rm=TRUE)*(1-sourcemask.conf.lim))]<-0
      }
    }
    #Cutup again
    if (do.sky.est||get.sky.rms||blank.cor) {
      if (cutup) {
        mask.stamp<-list(NULL)
        for (i in 1:npos) {
          mask.stamp[[i]]<-sm[mask.stamp.lims[i,1]:mask.stamp.lims[i,2],mask.stamp.lims[i,3]:mask.stamp.lims[i,4]]
        }
      } else {
        image.env$imm.dimim<-array(1, dim=dimim)
        image.env$imm.dimim[data.xrange[1]:data.xrange[2],data.yrange[1]:data.yrange[2]]<-sm[mask.xrange[1]:mask.xrange[2],mask.yrange[1]:mask.yrange[2]]
        image.env$imm<-sm
        mask.stamp<-sm
      }
    }
    if (!quiet) { cat(" - Done\n") }
    # /*fend*/ }}}
    #-----Diagnostic-----# /*fold*/ {{{
    if (diagnostic) {
      message(paste("SourceMask Max/Min:",max(sm),min(sm)))
      message(paste("OLDMethod - SourceMask Max/Min:",max(1-image.env$fa),min(1-image.env$fa)))
    }# /*fend*/ }}}
    #If wanted, output the SourceMask /*fold*/ {{{
    if (!exists("sourcemask.out")) { sourcemask.out<-FALSE }
    if (sourcemask.out){
      if (!quiet) { cat(paste('Outputting Source Mask to',sourcemask.filename,"   ")) }
      write.fits.image.file(file.path(path.root,path.work,path.out,sourcemask.filename),sm,image.env$data.hdr,nochange=TRUE)
      if (!quiet) { cat(" - Done\n") }
    }# /*fend*/ }}}
    #If we want the SourceMask only, then end here /*fold*/ {{{
    if (sourcemask.only) {
      if (!quiet) { cat("SourceMaskOnly Flag Set\n")  }
      return()
    } # /*fend*/ }}}
  } else if (plot.sample) {
    #Set the mask to 1 everywhere /*fold*/ {{{
    sm<-1
    mask.stamp<-1
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #Remove arrays that are no longer needed /*fold*/ {{{
  rm(fa, envir=image.env)
  gc()
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART THREE:  DEBLENDING /*fold*/ {{{
  #If wanted, perform sky estimation. Otherwise set to NA /*fold*/ {{{
  if (do.sky.est||get.sky.rms) {
    #Get sky estimates /*fold*/ {{{
    if (!quiet) { message("Perfoming Sky Estimation"); cat("Performing Sky Estimation") }
    #Perform Sky Estimation /*fold*/ {{{
    if (cutup) {
      timer<-system.time(skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                      cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,data.stamp=data.stamp,mask.stamp=mask.stamp,
                      clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts))
    } else {
      timer<-system.time(skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                      cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                      data.stamp=image.env$im, mask.stamp=image.env$imm.dimim,
                      clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts))
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
    if(!is.na(as.numeric(sky.default))) {
      sky.default<-as.numeric(sky.default)
      rmsdefault<-0.0
      errdefault<-0.0
      message("Using a default value for Sky Estimaes that have failed. Error on this value is unknown, so is assumed 0. Similarly, RMS is assumed 0")
    } else if (grepl("median",sky.default)) {
      sky.default<-median(skylocal,na.rm=TRUE)
      rmsdefault<-median(skyrms,na.rm=TRUE)
      errdefault<-median(skyerr,na.rm=TRUE)
    } else if (grepl("mean",sky.default)) {
      sky.default<-mean(skylocal, na.rm=TRUE)
      rmsdefault<-mean(skyrms, na.rm=TRUE)
      errdefault<-mean(skyerr, na.rm=TRUE)
    }
    if (is.na(sky.default)) {
      warning("All sky estimates have failed, and so the requested sky default is NA.\nTo stop bad behaviour, sky values are all being set to 0")
      sky.default<-0
      rmsdefault<-0
      errdefault<-0
    }
    skylocal[which(is.na(skylocal))]<-sky.default
    skyrms[which(is.na(skyrms))]<-rmsdefault
    skyerr[which(is.na(skyerr))]<-errdefault
    # /*fend*/ }}}
    #Notify /*fold*/ {{{
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Sky Estimate - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
    #If wanted, plot some of the Sky Estimates /*fold*/ {{{
    if (plot.sample) {
      if (!quiet) { message("Plotting Sky Estimation"); cat("Plotting Sky Estimation") }
        if (cutup) {
          timer<-system.time(plot.sky.estimate(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,
                          data.stamp.lims=data.stamp.lims,
                          cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                          data.stamp=data.stamp,mask.stamp=mask.stamp,
                          clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm,plot.all=plot.all,
                          path=file.path(path.root,path.work,path.out),rem.mask=TRUE,toFile=TRUE))
        } else {
          timer<-system.time(plot.sky.estimate(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,
                          data.stamp.lims=data.stamp.lims,
                          cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                          data.stamp=image.env$im,mask.stamp=image.env$imm.dimim,
                          clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm,plot.all=plot.all,
                          path=file.path(path.root,path.work,path.out),rem.mask=TRUE,toFile=TRUE))
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
    detecthres<-foreach(sfam=sfa, srms=skyrms, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { 5*srms*sqrt(sum(sfam)) }
    if (magnitudes) {
      suppressWarnings(detecthres.mag<--2.5*(log10(detecthres)-log10(ab.vega.flux))+mag.zp)
    } else {
      detecthres.mag<-array(NA, dim=c(length(sfa)))
    }
    if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
    # /*fend*/ }}}
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
  timer=system.time(dbw<-make.deblend.weight.maps(outenv=environment()))
  if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
    message(paste('Make Deblended Weightmap - Done (',round(timer[3], digits=2),'sec )'))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Remove arrays that are no longer needed /*fold*/ {{{
  rm(wfa, envir=image.env)
  # /*fend*/ }}}
  #Finalise SFA Apertures /*fold*/ {{{
  sfabak<-sfa
  if (!exists("psf.weighted")) { psf.weighted<-FALSE }
  if (psf.filt & !psf.weighted) {
    # If we have convolved with a PSF, & not doing PSF Weighting, convert back to tophat apertures (now expanded by PSF convolution) /*fold*/ {{{
    # For each aperture, Binary filter at aplim*max(ap) /*fold*/ {{{
    sba<-foreach(sfam=sfa, .options.mpi=mpi.opts, .noexport=ls(envir=environment()), .export="ap.limit")%dopar%{
      apvals<-rev(sort(sfam))
      tempsum<-cumsum(apvals)
      tempfunc<-approxfun(tempsum,apvals)
      apLim<-tempfunc(ap.limit*max(tempsum, na.rm=TRUE))
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
    if (make.convolved.apertures.map) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make.aperture.map(outenv=environment(), sfa, dimim))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      fa.filename<-paste("final_",fa.filename,sep="")
      if (!quiet) { cat(paste('Outputting Re-Boxcar-ed All Convolved Apertures Mask to',fa.filename,"   ")) }
      timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,fa.filename),image.env$wfa,image.env$data.hdr,nochange=TRUE) )
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
  } else if (psf.filt & psf.weighted) {
    #If we have convolved by the PSF & want PSF Weighting /*fold*/ {{{
    #If wanted, output the Convolved & Weighted Aperture Mask /*fold*/ {{{
    if (make.convolved.apertures.map) {
      if (!quiet) { cat(paste('Updated Apertures; ')) }
      timer=system.time(image.env$wfa<-make.aperture.map(outenv=environment(), sfa, dimim))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Make FA Mask - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
      fa.filename<-paste("final_",fa.filename,sep="")
      if (!quiet) { cat(paste('Outputting Final All Convolved Apertures Mask to',fa.filename,"   ")) }
      timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,fa.filename),image.env$wfa,image.env$data.hdr,nochange=TRUE) )
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
  if (!(no.psf)&(plot.sample)) {
    #PSF with Contours /*fold*/ {{{
    CairoPNG(file.path(path.root,path.work,path.out,"PSF.png"))
    psfvals<-rev(sort(psf))
    tempsum<-cumsum(psfvals)
    tempfunc<-approxfun(tempsum,psfvals)
    psfLimit<-tempfunc(ap.limit*max(tempsum, na.rm=TRUE))
    suppressWarnings(image(log10(psf),main="PSF & Binary Contour Levels", asp=1,col=heat.colors(256),useRaster=ifelse(length(psf)>1E4,TRUE,FALSE)))
    contour(psf, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
    if (psf.filt) { contour(psf, levels=c(psfLimit), labels=c(ap.limit), col='green', add=TRUE) }
    dev.off()
    # /*fend*/ }}}
    #Plot Example of PS minimum aperture Correction; minimum source /*fold*/ {{{
    ind<-(which(cat.a==0)[1])
    # /*fend*/ }}}
    #If no point sources, use the minimum aperture /*fold*/ {{{
    if (is.na(ind)) { ind<-which.min(cat.a) }
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
    psf.obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Make expanded grid of new pixel centres /*fold*/ {{{
    expanded<-expand.grid(seq(1,lenx),seq(1,leny))
    xnew<-expanded[,1]-cat.x[ind]%%1
    ynew<-expanded[,2]-cat.y[ind]%%1
    # /*fend*/ }}}
    #Interpolate /*fold*/ {{{
    ap<-matrix(interp.2d(xnew, ynew, psf.obj)[,3], ncol=leny,nrow=lenx)
    # /*fend*/ }}}
    #Aperture Correction Plot /*fold*/ {{{
    CairoPNG(file.path(path.root,path.work,path.out,"ApertureCorrection.png"))
    Rast<-ifelse(length(ap)>1E4,TRUE,FALSE)
    if (!psf.weighted) {
      #Binary Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (smallest source)", asp=1,col=heat.colors(256),useRaster=Rast))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      contour(ap, levels=c(psfLimit), labels=c(ap.limit), col='green', add=TRUE)
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=Rast))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    } else {
      #PSF Weighted Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (smallest source)", asp=1,col=heat.colors(256),useRaster=Rast))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=Rast))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
    #Plot Example of PS minimum aperture Correction; median source /*fold*/ {{{
    ind<-(which.min(abs(cat.a-median(cat.a)))[1])
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
    psf.obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
    # /*fend*/ }}}
    #Make expanded grid of new pixel centres /*fold*/ {{{
    expanded<-expand.grid(seq(1,lenx),seq(1,leny))
    xnew<-expanded[,1]-cat.x[ind]%%1
    ynew<-expanded[,2]-cat.y[ind]%%1
    # /*fend*/ }}}
    #Interpolate /*fold*/ {{{
    ap<-matrix(interp.2d(xnew, ynew, psf.obj)[,3], ncol=leny,nrow=lenx)
    # /*fend*/ }}}
    # /*fend*/ }}}
    #Aperture Correction Plot /*fold*/ {{{
    Rast<-ifelse(length(ap)>1E4,TRUE,FALSE)
    if (!psf.weighted) {
      #Binary Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap),main="Example: Minimum Aperture Correction (median source)", asp=1,col=heat.colors(256),useRaster=Rast))
      contour(ap, levels=(tempfunc(c(0.5,0.9,0.95,0.99,0.999,0.9999)*max(tempsum))), labels=c(0.5,0.9,0.95,0.99,0.999,0.9999), col='blue', add=TRUE)
      if (psf.filt) { contour(ap, levels=c(psfLimit), labels=c(ap.limit), col='green', add=TRUE) }
      suppressWarnings(image(log10(sfa[[ind]]),col=col2alpha('blue',0.3),add=TRUE,useRaster=Rast))
      spsf<-sum(ap)
      ssfap<-sum(sfa[[ind]]*ap)
      label("topright",lab=paste("SumPSF=",round(spsf,digits=2),"\nSum(PSF*Ap)=",round(ssfap,digits=2),"\nApCorr=",round(spsf/ssfap,digits=2),sep=""))
      # /*fend*/ }}}
    } else {
      #PSF Weighted Aperture /*fold*/ {{{
      suppressWarnings(image(log10(ap)-log10(sfa[[ind]]),main="Example: Minimum Aperture Correction (median source; shown as residual)", asp=1,col=heat.colors(256),useRaster=Rast))
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
  dfa<-foreach(dbwm=dbw, sfam=sfa, .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
  # /*fend*/ }}}
  #Notify /*fold*/ {{{
  if (showtime) { cat("   - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )\n")
    message(paste("Make DFA - Done (",round(proc.time()[3]-timer[3],digits=2),"sec )"))
  } else if (!quiet) { cat("   - Done\n") }
  # /*fend*/ }}}
  #Check if we're just calculating the deblend fraction /*fold*/ {{{
  if (!exists("get.debl.frac")) { get.debl.frac<-FALSE }
  if (get.debl.frac) {
    if (!quiet) { cat("Getting Deblend Fraction Only flag is TRUE.\nCalculating Aperture Integrals & Outputting catalogue.\n") }
    #Integral of the aperture; ssa /*fold*/ {{{
    if (verbose) { cat("Integral of the aperture") }
    ssa<-foreach(sam=sa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the convolved aperture; ssfa /*fold*/ {{{
    if (verbose) { cat("Integral of the convolved aperture") }
    ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the deblended convolved aperture; sdfa /*fold*/ {{{
    if (verbose) { cat("Integral of the deblended convolved aperture") }
    sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(dfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    spsf<-rep(sumpsf, length(sdfa))
    #If wanted, output the Results Table /*fold*/ {{{
    if (write.tab) {
      #Output the results table /*fold*/ {{{
      if ((loop.total!=1)&&(length(param.env$tableout.name)!=loop.total)) {
        if (!quiet) { cat(paste('Writing Deblend Fraction Table to ',file.path(path.root,path.work,path.out,paste(sep="",tableout.name,"_file",f,".csv")),'   ')) }
        timer=system.time(write.deblend.fraction.table(filename=file.path(path.root,path.work,path.out,paste(sep="",tableout.name,"_file",f,".csv"))) )
      } else {
        if (!quiet) { cat(paste('Writing Deblend Fraction Table to ',file.path(path.root,path.work,path.out,paste(sep="",tableout.name,".csv")),'   ')) }
        timer=system.time(write.deblend.fraction.table(filename=file.path(path.root,path.work,path.out,paste(sep="",tableout.name,".csv"))) )
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
  if (iterate.fluxes) {
    #Notify /*fold*/ {{{
    message(paste("Iterating Flux Determination",num.iterations,"times {\n"))
    if (!quiet) { cat("Iterating Flux Determination",num.iterations,"times {\n") }
    # /*fend*/ }}}
    #For the number of desired iterations /*fold*/ {{{
    #Setup for iteration /*fold*/ {{{
    weight.type='flux'
    quietbak<-quiet
    psf.filtbak<-psf.filt
    fluxiters<-matrix(as.numeric(NA),ncol=num.iterations,nrow=length(cat.id))
    erriters<-fluxiters
    sdfaiters<-fluxiters
    #Integral of the convolved aperture; ssfa /*fold*/ {{{
    if (verbose) { cat("  Integral of the convolved aperture") }
    ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    #Integral of the deblended convolved aperture; sdfa /*fold*/ {{{
    if (verbose) { cat("  Integral of the deblended convolved aperture") }
    sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(dfam) }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
    sdfad<-sdfa*NA
    sdfa2e2<-sdfa*NA
    attach(image.env)
    #/*fend*/ }}}
    for (iter in 1:num.iterations) {
      #Determine objects to iterate over /*fold*/ {{{
      if (iter!=1) {
        xind<-which(sdfa!=ssfa & sdfa!=0)
        iterateLost<-sdfa==0
      } else {
        xind<-1:length(sdfa)
      }
      #Check we have something to do! /*fold*/ {{{
      if (length(xind)==0) {
        if (verbose) { cat(" - Breaking out; no objects remaining overlap, so no point iterating!\n") }
        break
      }
      #/*fend*/ }}}
      #Update MPI options /*fold*/ {{{
      chunk.size=ceiling(length(xind)/getDoParWorkers())
      mpi.opts<-list(chunkSize=chunk.size)
      #/*fend*/ }}}
      #/*fend*/ }}}
      #Calculate the flux per object /*fold*/ {{{
      #Notify /*fold*/ {{{
      message(paste("Calculating Flux (#",iter,")"))
      if (!quiet) { cat("  Calculating Flux (#",iter,")") }
      # /*fend*/ }}}
      timer<-proc.time()
      if (cutup) {
        sdfad[xind]<-foreach(dfam=dfa[xind],im=data.stamp[xind], xlo=ap.lims.data.stamp[xind,1],xup=ap.lims.data.stamp[xind,2], ylo=ap.lims.data.stamp[xind,3],yup=ap.lims.data.stamp[xind,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
           sum(dfam*(im[xlo:xup,ylo:yup]))
        }
      } else {
        sdfad[xind]<-foreach(dfam=dfa[xind], xlo=ap.lims.data.map[xind,1],xup=ap.lims.data.map[xind,2], ylo=ap.lims.data.map[xind,3],yup=ap.lims.data.map[xind,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment()),.export='im') %dopar% {
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
      if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
        sdfa2e2[xind]<-foreach(dfam=dfa[xind], xlo=ap.lims.error.stamp[xind,1],xup=ap.lims.error.stamp[xind,2],ylo=ap.lims.error.stamp[xind,3],yup=ap.lims.error.stamp[xind,4], ime=error.stamp[xind], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
            sum((dfam*ime[xlo:xup,ylo:yup])^2.)
        }
      } else if (length(error.stamp)==1){
        if (error.stamp==1) {
          sdfa2e2[xind]<-foreach(dfam=dfa[xind], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
            sum((dfam)^2)
          }
        } else if (error.stamp==0) {
          sdfa2e2<-rep(0,length(sdfad))
        } else {
          sdfa2e2[xind]<-foreach(dfam=dfa[xind], .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
              sum((dfam*error.stamp)^2.)
          }
        }
      } else {
        sdfa2e2[xind]<-foreach(dfam=dfa[xind], xlo=ap.lims.error.map[xind,1],xup=ap.lims.error.map[xind,2],ylo=ap.lims.error.map[xind,3],yup=ap.lims.error.map[xind,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
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
      sdfa[xind]<-foreach(dfam=dfa[xind], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
      #Subtract Sky Flux /*fold*/ {{{
      if (do.sky.est) {
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
      if (psf.weighted) {
        #If using psf.weighted, we may be able to skip an unnecessary convolution /*fold*/ {{{
        psf.filt<-FALSE
        timer=system.time(wsfa[xind]<-make.convolved.apertures(outenv=environment(), sfa,flux.weightin=sdfad,subs=xind))
        # /*fend*/ }}}
      } else {
        #If not using psf.weighted, we may need to re-convolve /*fold*/ {{{
        psf.filt<-psf.filtbak
        timer=system.time(wsfa[xind]<-make.convolved.apertures(outenv=environment(), sa,flux.weightin=sdfad,subs=xind))
        # /*fend*/ }}}
      }
      timer2=system.time(image.env$wfa<-make.aperture.map(outenv=environment(), wsfa, dimim,subs=xind))
      psf.filt<-psf.filtbak
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
      timer2=system.time(dbw[xind]<-make.deblend.weight.maps(outenv=environment(),subs=xind))
      quiet<-quietbak
      timer<-proc.time()
      dfa[xind]<-foreach(dbwm=dbw[xind], sfam=sfa[xind], .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { dbwm*sfam }
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
    chunk.size=ceiling(length(cat.id)/getDoParWorkers())
    mpi.opts<-list(chunkSize=chunk.size)
    #/*fend*/ }}}
    #/*fend*/ }}}
    detach(image.env)
    if (!quiet) { cat("} Done\n") }
  # /*fend*/ }}}
  }
  # /*fend*/ }}}
  #If wanted, output the Deblended & Convolved Aperture Mask /*fold*/ {{{
  if (make.debelended.apertures.map) {
    if (!quiet) { cat(paste('Making All Deblended Convolved Apertures Mask - ')) }
    timer=system.time(image.env$adfa<-make.aperture.map(outenv=environment(), dfa, dimim))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Make ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (!quiet) { cat(paste('Outputting All Deblended Convolved Apertures Mask to',dfa.filename,"   ")) }
    timer=system.time(write.fits.image.file(file.path(path.root,path.work,path.out,dfa.filename),image.env$adfa,image.env$data.hdr,nochange=TRUE) )
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Output ADFA Mask - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
  }# /*fend*/ }}}
  #Remove array that is no longer needed /*fold*/ {{{
  if (!plot.sample) { rm(dbw) }
  rm(adfa, envir=image.env)
  gc()
  # /*fend*/ }}}
  # /*fend*/ }}}
  #If wanted, Perform Randoms Correction /*fold*/ {{{
  if (!exists("ran.cor")) { ran.cor<-FALSE }
  if (ran.cor=="execute") {
    .executeran.cor()
    ran.cor<-TRUE
  }
  if (ran.cor) {
    if (!quiet) { message("Perfoming Randoms Correction"); cat("Performing Randoms Correction") }
    if (cutup) {
      timer<-system.time(randoms<-ran.cor(data.stamp=data.stamp,mask.stamp=mask.stamp,ap.stamp=sfa,ap.stamp.lims=ap.lims.data.stamp,numIters=num.randoms,mpi.opts=mpi.opts,rem.mask=FALSE))
    } else {
      timer<-system.time(randoms<-ran.cor(data.stamp=image.env$im[data.stamp.lims[1,1]:data.stamp.lims[1,2],data.stamp.lims[1,3]:data.stamp.lims[1,4]],
                      mask.stamp=image.env$imm[mask.stamp.lims[1,1]:mask.stamp.lims[1,2],mask.stamp.lims[1,3]:mask.stamp.lims[1,4]],
      ap.stamp=sfa,ap.stamp.lims=ap.lims.data.map,numIters=num.randoms,mpi.opts=mpi.opts,rem.mask=FALSE))
    }
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Randoms Correction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (plot.sample) {
      if (!quiet) { message("Plotting Randoms Correction"); cat("Plotting Randoms Correction") }
      if (cutup) {
        timer<-system.time(plot.ran.cor(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,data.stamp=data.stamp,mask.stamp=mask.stamp,ap.stamp=sfa,ap.stamp.lims=ap.lims.data.stamp,data.stamp.lims=data.stamp.lims,numIters=num.randoms,rem.mask=FALSE,path=file.path(path.root,path.work,path.out),plot.all=plot.all,toFile=TRUE))
      } else {
        timer<-system.time(plot.ran.cor(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,data.stamp=image.env$im[data.stamp.lims[1,1]:data.stamp.lims[1,2],data.stamp.lims[1,3]:data.stamp.lims[1,4]],
                      mask.stamp=image.env$imm[mask.stamp.lims[1,1]:mask.stamp.lims[1,2],mask.stamp.lims[1,3]:mask.stamp.lims[1,4]],ap.stamp=sfa,ap.stamp.lims=ap.lims.data.map,data.stamp.lims=data.stamp.lims,numIters=num.randoms,rem.mask=FALSE,path=file.path(path.root,path.work,path.out),plot.all=plot.all,toFile=TRUE))
      }
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Plotting Randoms Correction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
  }
  # /*fend*/ }}}
  #If wanted, Perform Blanks Correction /*fold*/ {{{
  if (!exists("blank.cor")) { blank.cor<-FALSE }
  if (blank.cor) {
    if (!quiet) { message("Perfoming Blanks Correction"); cat("Performing Blanks Correction") }
    if (cutup) {
      timer<-system.time(blanks<-ran.cor(data.stamp=data.stamp,mask.stamp=mask.stamp,ap.stamp=sfa,ap.stamp.lims=ap.lims.data.stamp,numIters=num.blanks,mpi.opts=mpi.opts,rem.mask=TRUE))
    } else {
      timer<-system.time(blanks<-ran.cor(data.stamp=image.env$im[data.stamp.lims[1,1]:data.stamp.lims[1,2],data.stamp.lims[1,3]:data.stamp.lims[1,4]],
                      mask.stamp=image.env$imm[mask.stamp.lims[1,1]:mask.stamp.lims[1,2],mask.stamp.lims[1,3]:mask.stamp.lims[1,4]],ap.stamp=sfa,ap.stamp.lims=ap.lims.data.map,numIters=num.blanks,mpi.opts=mpi.opts,rem.mask=TRUE))
    }
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Blanks Correction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    if (plot.sample) {
      if (!quiet) { message("Plotting Blanks Correction"); cat("Plotting Blanks Correction") }
    if (cutup) {
      timer<-system.time(plot.ran.cor(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,data.stamp=data.stamp,mask.stamp=mask.stamp,ap.stamp=sfa,ap.stamp.lims=ap.lims.data.stamp,data.stamp.lims=data.stamp.lims,numIters=num.blanks,rem.mask=TRUE,path=file.path(path.root,path.work,path.out),plot.all=plot.all,toFile=TRUE))
    } else {
      timer<-system.time(plot.ran.cor(cat.id=cat.id,cat.x=cat.x,cat.y=cat.y,data.stamp=image.env$im[data.stamp.lims[1,1]:data.stamp.lims[1,2],data.stamp.lims[1,3]:data.stamp.lims[1,4]],
                      mask.stamp=image.env$imm[mask.stamp.lims[1,1]:mask.stamp.lims[1,2],mask.stamp.lims[1,3]:mask.stamp.lims[1,4]],ap.stamp=sfa,ap.stamp.lims=ap.lims.data.map,data.stamp.lims=data.stamp.lims,numIters=num.blanks,rem.mask=TRUE,path=file.path(path.root,path.work,path.out),plot.all=plot.all,toFile=TRUE))
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
    if (length(which(is.na(error.stamp)))!=0) {        message(paste("Input Parameter 'ime' contains NA elements")) }
    if (length(which(is.na(sfa)))!=0) {        message(paste("Input Parameter 'sfa' contains NA elements")) }
    if (length(which(is.na(dfa)))!=0) {        message(paste("Input Parameter 'dfa' contains NA elements")) }
    if (length(which(is.na(flux.corr)))!=0) {        message(paste("Input Parameter 'flux.corr' contains NA elements")) }
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
      saturated<-foreach(sfam=sfa, im=data.stamp, xlo=ap.lims.data.stamp[,1],xup=ap.lims.data.stamp[,2], ylo=ap.lims.data.stamp[,3],yup=ap.lims.data.stamp[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment()), .export='saturation') %dopar% {
          any(im[xlo:xup,ylo:yup][which(sfam!=0,arr.ind=TRUE)]>=saturation)
      }
    } else {
      saturated<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export=c('im','saturation'), xlo=ap.lims.data.map[,1],xup=ap.lims.data.map[,2], ylo=ap.lims.data.map[,3],yup=ap.lims.data.map[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          any(im[xlo:xup,ylo:yup][which(sfam!=0,arr.ind=TRUE)]>=saturation)
      }
    }
    if (verbose) { cat(" - Done\n") }
  } else {
    saturated<-rep(FALSE,length(cat.id))
  }
  # /*fend*/ }}}
#-----
  #Image Flux at central pixel; pixflux /*fold*/ {{{
  if (!use.pixel.fluxweight) {
    if (verbose) { cat("      Image Flux at central pixel ") }
    if (cutup) {
      pixflux<-foreach(xp=x.pix-(data.stamp.lims[,1]-1),yp=y.pix-(data.stamp.lims[,3]-1), im=data.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
            im[xp,yp]
      }
    } else {
      pixflux<-foreach(xp=x.pix,yp=y.pix, .noexport=ls(envir=environment()), .export=c('im'), .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
            im[xp,yp]
      }
    }
    if (verbose) { cat(" - Done\n") }
  }
  # /*fend*/ }}}
#-----
  #Integral of the aperture; ssa /*fold*/ {{{
  if (verbose) { cat("      Integral of the aperture") }
  ssa<-foreach(sam=sa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sam) }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the convolved aperture; ssfa /*fold*/ {{{
  if (verbose) { cat("      Integral of the convolved aperture") }
  ssfa<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sfam) }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(unmodified convolved aperture)*(modified convolved aperture)]; ssfau /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(unmodified convolved aperture)*(convolved aperture)]") }
  ssfau<-foreach(sfam=sfa, usfam=sfabak, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(sfam*usfam)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * image); ssfad /*fold*/ {{{
  if (verbose) { cat("      Integral of the (convolved aperture * image)") }
  if (cutup) {
    ssfad<-foreach(sfam=sfa, im=data.stamp, xlo=ap.lims.data.stamp[,1],xup=ap.lims.data.stamp[,2], ylo=ap.lims.data.stamp[,3],yup=ap.lims.data.stamp[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup]))
    }
  } else {
    ssfad<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export='im', xlo=ap.lims.data.map[,1],xup=ap.lims.data.map[,2], ylo=ap.lims.data.map[,3],yup=ap.lims.data.map[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum(sfam*(im[xlo:xup,ylo:yup]))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Quartered Photometry - Quartered Integral of the (convolved aperture * image); qssfad /*fold*/ {{{
  if (verbose) { cat("      Quartered Integral of the (convolved aperture * image)") }
  if (cutup) {
    qssfad<-foreach(sfam=sfa, im=data.stamp, x1=ap.lims.data.stamp[,1], x2=ap.lims.data.stamp[,1]+floor((ap.lims.data.stamp[,2]-ap.lims.data.stamp[,1])/2), x3=ap.lims.data.stamp[,1]+floor((ap.lims.data.stamp[,2]-ap.lims.data.stamp[,1])/2)+1, x4=ap.lims.data.stamp[,2],
                                          y1=ap.lims.data.stamp[,3], y2=ap.lims.data.stamp[,3]+floor((ap.lims.data.stamp[,4]-ap.lims.data.stamp[,3])/2), y3=ap.lims.data.stamp[,3]+floor((ap.lims.data.stamp[,4]-ap.lims.data.stamp[,3])/2)+1, y4=ap.lims.data.stamp[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpi.opts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
        cbind(sum(sfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(sfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(sfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(sfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  } else {
    qssfad<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export='im', x1=ap.lims.data.map[,1], x2=ap.lims.data.map[,1]+floor((ap.lims.data.map[,2]-ap.lims.data.map[,1])/2), x3=ap.lims.data.map[,1]+floor((ap.lims.data.map[,2]-ap.lims.data.map[,1])/2)+1, x4=ap.lims.data.map[,2],
                                          y1=ap.lims.data.map[,3], y2=ap.lims.data.map[,3]+floor((ap.lims.data.map[,4]-ap.lims.data.map[,3])/2), y3=ap.lims.data.map[,3]+floor((ap.lims.data.map[,4]-ap.lims.data.map[,3])/2)+1, y4=ap.lims.data.map[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpi.opts, .combine="rbind") %dopar% {
        cbind(sum(sfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(sfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(sfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(sfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the reinterpolated psf; spsf /*fold*/ {{{
  if (!no.psf) {
    if (verbose) { cat("      Integral of the reinterpolated psf") }
    spsf<-rep(sum(psf),length(cat.id))
    if (verbose) { cat(" - Done\n") }
  } else {
    spsf<-rep(NA, length(sfa))
  }
  # /*fend*/ }}}
#-----
  #Integral of the (convolved aperture * psf); ssfap /*fold*/ {{{
  if (!no.psf) {
    if (verbose) { cat("      Integral of the (convolved aperture * psf)") }
    ssfap<-foreach(sfam=sfa, slen=stamplen, xc=cat.x, yc=cat.y, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment()), .export="psf") %dopar% {
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
      psf.obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[(1:(lenxpsf+1)+dx-1)%%(lenxpsf+1),(1:(lenypsf+1)+dy-1)%%(lenypsf+1)][1:lenx,1:leny])
      # /*fend*/ }}}
      #Make expanded grid of new pixel centres /*fold*/ {{{
      expanded<-expand.grid(seq(1,lenx),seq(1,leny))
      xnew<-expanded[,1]-xc%%1
      ynew<-expanded[,2]-yc%%1
      # /*fend*/ }}}
      #Interpolate /*fold*/ {{{
      ap<-matrix(interp.2d(xnew, ynew, psf.obj)[,3], ncol=leny,nrow=lenx)
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
  sdfa<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum(dfam)  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the (deblended convolved aperture * image); sdfad /*fold*/ {{{
  if (verbose) { cat("      Integral of the (deblended convolved aperture * image)") }
  if (cutup) {
    sdfad<-foreach(dfam=dfa,im=data.stamp, xlo=ap.lims.data.stamp[,1],xup=ap.lims.data.stamp[,2], ylo=ap.lims.data.stamp[,3],yup=ap.lims.data.stamp[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum(dfam*(im[xlo:xup,ylo:yup]))
    }
  } else {
    sdfad<-foreach(dfam=dfa,.noexport=ls(envir=environment()), .export='im', xlo=ap.lims.data.map[,1],xup=ap.lims.data.map[,2], ylo=ap.lims.data.map[,3],yup=ap.lims.data.map[,4], .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum(dfam*(im[xlo:xup,ylo:yup]))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Quartered Photometry - Quartered Integral of the (deblended convolved aperture * image); qsdfad /*fold*/ {{{
  if (verbose) { cat("      Quartered Integral of the (deblended convolved aperture * image)") }
  if (cutup) {
    qsdfad<-foreach(dfam=dfa, im=data.stamp, x1=ap.lims.data.stamp[,1], x2=ap.lims.data.stamp[,1]+floor((ap.lims.data.stamp[,2]-ap.lims.data.stamp[,1])/2), x3=ap.lims.data.stamp[,1]+floor((ap.lims.data.stamp[,2]-ap.lims.data.stamp[,1])/2)+1, x4=ap.lims.data.stamp[,2],
                                          y1=ap.lims.data.stamp[,3], y2=ap.lims.data.stamp[,3]+floor((ap.lims.data.stamp[,4]-ap.lims.data.stamp[,3])/2), y3=ap.lims.data.stamp[,3]+floor((ap.lims.data.stamp[,4]-ap.lims.data.stamp[,3])/2)+1, y4=ap.lims.data.stamp[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpi.opts, .combine="rbind", .noexport=ls(envir=environment())) %dopar% {
        cbind(sum(dfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(dfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(dfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(dfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  } else {
    qsdfad<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export='im', x1=ap.lims.data.map[,1], x2=ap.lims.data.map[,1]+floor((ap.lims.data.map[,2]-ap.lims.data.map[,1])/2), x3=ap.lims.data.map[,1]+floor((ap.lims.data.map[,2]-ap.lims.data.map[,1])/2)+1, x4=ap.lims.data.map[,2],
                                          y1=ap.lims.data.map[,3], y2=ap.lims.data.map[,3]+floor((ap.lims.data.map[,4]-ap.lims.data.map[,3])/2), y3=ap.lims.data.map[,3]+floor((ap.lims.data.map[,4]-ap.lims.data.map[,3])/2)+1, y4=ap.lims.data.map[,4],
                                         sx1=rep(1,length(stamplen)),sx2=1+floor(stamplen-1)/2,sx3=1+floor(stamplen-1)/2+1,sx4=stamplen,
                                         sy1=rep(1,length(stamplen)),sy2=1+floor(stamplen-1)/2,sy3=1+floor(stamplen-1)/2+1,sy4=stamplen,
                                          .inorder=TRUE, .options.mpi=mpi.opts, .combine="rbind") %dopar% {
        cbind(sum(dfam[sx1:sx2,sy1:sy2]*(im[x1:x2,y1:y2])),sum(dfam[sx3:sx4,sy1:sy2]*(im[x3:x4,y1:y2])),sum(dfam[sx1:sx2,sy3:sy4]*(im[x1:x2,y3:y4])),sum(dfam[sx3:sx4,sy3:sy4]*(im[x3:x4,y3:y4])))
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(convolved aperture * image error)^2]; ssfa2e2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(convolved aperture * image error)^2]") }
  if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
    ssfa2e2<-foreach(sfam=sfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum((sfam*ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(error.stamp)==1){
    if (error.stamp==1) {
      ssfa2e2<-ssfa2
    } else if (error.stamp==0) {
      ssfa2e2<-rep(0,length(ssfa))
    } else {
      ssfa2e2<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum((sfam*error.stamp)^2.)
      }
    }
  } else {
    ssfa2e2<-foreach(sfam=sfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum((sfam*ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
  #Integral of the [(deblended convolved aperture * image error)^2]; sdfa2e2 /*fold*/ {{{
  if (verbose) { cat("      Integral of the [(deblended convolved aperture * image error)^2]") }
  if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
    sdfa2e2<-foreach(dfam=dfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
  } else if (length(error.stamp)==1){
    if (error.stamp==1) {
      sdfa2e2<-sdfa2
    } else if (error.stamp==0) {
      sdfa2e2<-rep(0,length(ssfa))
    } else {
      sdfa2e2<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum((dfam*error.stamp)^2.)
      }
    }
  } else {
    sdfa2e2<-foreach(dfam=dfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum((dfam*ime[xlo:xup,ylo:yup])^2.)
    }
  }
  if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
#-----
#----- DIAGNOSTIC ITEMS -----# /*fold*/ {{{
  if (diagnostic) {
  #-----
    #Integral of the (convolved aperture * (image error)^2); ssfae2 /*fold*/ {{{
    if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
    if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
      ssfae2<-foreach(sfam=sfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
          sum(sfam*(ime[xlo:xup,ylo:yup])^2.)
      }
    } else if (length(error.stamp)==1){
      if (error.stamp==1) {
        ssfae2<-ssfa
      } else if (error.stamp==0) {
        ssfae2<-rep(0,length(ssfa))
      } else {
        ssfae2<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
            sum(sfam*(error.stamp)^2.)
        }
      }
    } else {
      ssfae2<-foreach(sfam=sfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum(sfam*(ime[xlo:xup,ylo:yup])^2.)
      }
    }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
    #Integral of the [deblended convolved aperture * (image error)^2]; sdfae2 /*fold*/ {{{
    if (verbose) { cat("      Integral of the [deblended convolved aperture * (image error)^2]") }
    if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
      sdfae2<-foreach(dfam=dfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
          sum(dfam*(ime[xlo:xup,ylo:yup])^2.)
      }
    } else if (length(error.stamp)==1){
      if (error.stamp==1) {
        sdfae2<-sdfa2
      } else if (error.stamp==0) {
        sdfae2<-rep(0,length(ssfa))
      } else {
        sdfae2<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
            sum(dfam*(error.stamp)^2.)
        }
      }
    } else {
      sdfae2<-foreach(dfam=dfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum(dfam*(ime[xlo:xup,ylo:yup])^2.)
      }
    }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
    #Integral of the (convolved aperture * image error); ssfae /*fold*/ {{{
    if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
    if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
      ssfae<-foreach(sfam=sfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum(sfam*ime[xlo:xup,ylo:yup])
      }
    } else if (length(error.stamp)==1){
      if (error.stamp==1) {
        ssfae<-ssfa
      } else if (error.stamp==0) {
        ssfae<-rep(0,length(ssfa))
      } else {
        ssfae<-foreach(sfam=sfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum(sfam*error.stamp)
        }
      }
    } else {
      ssfae<-foreach(sfam=sfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum(sfam*ime[xlo:xup,ylo:yup])
      }
    }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
    #Integral of the (deblended convolved aperture * image error); sdfae /*fold*/ {{{
    if (verbose) { cat("      Integral of the (convolved aperture * image error)") }
    if (length(error.stamp)>1|(length(error.stamp)==1 & is.list(error.stamp))) {
      sdfae<-foreach(dfam=dfa, xlo=ap.lims.error.stamp[,1],xup=ap.lims.error.stamp[,2],ylo=ap.lims.error.stamp[,3],yup=ap.lims.error.stamp[,4], ime=error.stamp, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% {
        sum(dfam*ime[xlo:xup,ylo:yup])
      }
    } else if (length(error.stamp)==1){
      if (error.stamp==1) {
        sdfae<-sdfa
      } else if (error.stamp==0) {
        sdfae<-rep(0,length(sdfa))
      } else {
        sdfae<-foreach(dfam=dfa, .noexport=ls(envir=environment()), .export="error.stamp", .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
          sum(dfam*error.stamp)
        }
      }
    } else {
      sdfae<-foreach(dfam=dfa, xlo=ap.lims.error.map[,1],xup=ap.lims.error.map[,2],ylo=ap.lims.error.map[,3],yup=ap.lims.error.map[,4], .noexport=ls(envir=environment()), .export='ime', .inorder=TRUE, .combine='c', .options.mpi=mpi.opts) %dopar% {
        sum(dfam*ime[xlo:xup,ylo:yup])
      }
    }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
    #Integral of the [(convolved aperture)^2]; ssfa2 /*fold*/ {{{
    if (verbose) { cat("      Integral of the [(convolved aperture)^2]") }
    ssfa2<-foreach(sfam=sfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum((sfam)^2.)  }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
    #Integral of the [(deblended convolved aperture)^2]; sdfa2 /*fold*/ {{{
    if (verbose) { cat("      Integral of the [(deblended convolved aperture)^2]") }
    sdfa2<-foreach(dfam=dfa, .inorder=TRUE, .combine='c', .options.mpi=mpi.opts, .noexport=ls(envir=environment())) %dopar% { sum((dfam)^2.)  }
    if (verbose) { cat(" - Done\n") } # /*fend*/ }}}
  #-----
  }
#/*fend*/ }}}
#-----
  #Final Flux & Error Calculations /*fold*/ {{{
  if (!cutup) { detach(image.env) }
  if (verbose) { cat("      Final Fluxes and Error Calculations ") }
  #Calculate Aperture-Weighting Correction /*fold*/ {{{
  WtCorr<-ssfa/ssfau
  # /*fend*/ }}}
  #Calculate Minimum Aperture Correction /*fold*/ {{{
  if (!no.psf) {
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
  sfaerr<-ssfa2e2*ApCorr^2 + ((conf*beamarea)^2.*sqrt(ssfa)*ApCorr)^2
  if (blank.cor) { sfaerr<-sfaerr+(blanks$randMean.MAD*ApCorr)^2 }
  # /*fend*/ }}}

  #Deblended Convolved aperture error /*fold*/ {{{
  dfaerr<-sdfa2e2*ApCorr^2 + ((conf*beamarea)^2.*sqrt(sdfa)*ApCorr)^2
  if (blank.cor) { dfaerr<-dfaerr+(blanks$randMean.MAD*ApCorr)^2 }
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
  if (!quiet) { cat("   Performing Final Calculations   ") }
  if (do.sky.est|get.sky.rms) {
    skyflux.mean<-skylocal.mean*sdfa
    skyflux<-skylocal*sdfa
    if (!blank.cor) {
      sfaerr[which(!is.na(skyrms))]<-(sfaerr+(skyrms*sqrt(ssfa)*ApCorr)^2)[which(!is.na(skyrms))]
      dfaerr[which(!is.na(skyrms))]<-(dfaerr+(skyrms*sqrt(sdfa)*ApCorr)^2)[which(!is.na(skyrms))]
    }
    if (do.sky.est) {
      if (!quiet) { message("Perfoming Sky Subtraction"); cat("{\n   Performing Sky Subtraction") }
      #Subtract Sky Flux /*fold*/ {{{
      dfaflux<-dfaflux-skyflux
      sfaflux<-sfaflux-skylocal*ssfa
      dfaerr[which(!is.na(skyerr))]<-(dfaerr+(skyerr*sdfa*ApCorr)^2)[which(!is.na(skyerr))]
      sfaerr[which(!is.na(skyerr))]<-(sfaerr+(skyerr*sdfa*ApCorr)^2)[which(!is.na(skyerr))]
      if (!quiet) { message(paste("   - Done\n")); cat("   - Done\n")}
      # /*fend*/ }}}
    }
  }

  #Deblend error /*fold*/ {{{
  deblerr<-((1-sdfa/ssfa)*(1/sqrt(12))*abs(sfaflux)*ApCorr)
  dfaerr<-dfaerr + (deblerr)^2
  # /*fend*/ }}}

  #Finalise Errors /*fold*/ {{{
  sfaerr<-sqrt(sfaerr)
  dfaerr<-sqrt(dfaerr)
  # /*fend*/ }}}
  # /*fend*/ }}}
  #Apply Aperture Correction to the Aperture Fluxess /*fold*/ {{{
  sfaflux<-sfaflux*ApCorr
  dfaflux<-dfaflux*ApCorr
  # /*fend*/ }}}
  #Apply Additional Flux Correction to the Finalised Values /*fold*/ {{{
  if (flux.corr!=1) {
    sfaflux<-sfaflux*flux.corr
    dfaflux<-dfaflux*flux.corr
    sfaerr<-sfaerr*flux.corr
    dfaerr<-dfaerr*flux.corr
  }# /*fend*/ }}}
  #Calculate magnitudes /*fold*/ {{{
  if (magnitudes) {
    suppressWarnings(mags<--2.5*(log10(dfaflux)-log10(ab.vega.flux))+mag.zp)
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
    message(paste("After assignment",round(length(which(is.na(pixflux)))/length(pixflux)*100,digits=2),"% of the pixflux matrix are NA"))
  }# /*fend*/ }}}
  if (!quiet) { cat("\n} Galaxy Results Complete\n") }
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART FIVE: OUTPUT /*fold*/ {{{
  #Do we want to plot a sample of the apertures? /*fold*/ {{{
  if (plot.sample) {
    #Set output name /*fold*/ {{{
    dir.create(file.path(path.root,path.work,path.out,"COGs"),showWarnings=FALSE)
    # /*fend*/ }}}
    #Determine Sample to Plot /*fold*/ {{{
    if (!exists("plot.all")) { plot.all<-FALSE }
    if (!plot.all) {
      if (!quiet) { message("Writing Sample of COGs to File"); cat("Writing Sample of COGs to File") }
      #Output a random 15 apertures to file /*fold*/ {{{
      ind1=which(cat.a>0 & sdfa!=0)
      ind1=ind1[order(sdfad[ind1],decreasing=TRUE)][1:15]
      ind2=which(sdfa!=0)[order(runif(1:length(which(sdfa!=0))))[1:15]]
      ind<-c(ind1, ind2)
      rm(ind1)
      rm(ind2)
      ind<-ind[which(!is.na(ind))]
      # /*fend*/ }}}
    } else {
      #All /*fold*/ {{{
      if (!quiet) { message("Writing All COGs to File"); cat("Writing All COGs to File") }
      ind=c(which(cat.a>0),which(cat.a==0))
      # /*fend*/ }}}
    }
    # /*fend*/ }}}
    for (i in ind) {
      #Open Device /*fold*/ {{{
      CairoPNG(file.path(path.root,path.work,path.out,paste("COGs/",cat.id[i],".png",sep="")),width=8*res,height=8*res,res=res)
      #pdf(file.path(path.root,path.work,path.out,paste("COGs/",cat.id[i],".pdf",sep="")),width=14,height=3.5)
      # /*fend*/ }}}
      #Set Layout /*fold*/ {{{
      mar<-par("mar")
      par(mar=mar*c(0.6,0.8,0.5,1))
      layout(rbind(c(1,2),c(3,4),c(0,0)),heights=c(1,1,0.05))
      # /*fend*/ }}}
      #Axes Limits /*fold*/ {{{
      xlims=stamplen[i]*c(-2/3,2/3)*arcsec.per.pix
      ylims=stamplen[i]*c(-2/3,2/3)*arcsec.per.pix
      # /*fend*/ }}}
      #Make Ap and Source Mask Block/Trans matricies /*fold*/ {{{
      apT<-sfa[[i]]
      apT[which(sfa[[i]]==0)]<-NA
      apT[which(sfa[[i]]!=0)]<-1
      apB<-sfa[[i]]
      apB[which(sfa[[i]]==0)]<-1
      apB[which(sfa[[i]]!=0)]<-NA
      if (sourcemask) {
        smB<-sm[mask.stamp.lims[i,1]:mask.stamp.lims[i,2],mask.stamp.lims[i,3]:mask.stamp.lims[i,4]]
        smB[which(smB==0)]<-NA
      } else {
        smB<-1
      }
      # /*fend*/ }}}
      #Plot Image in greyscale /*fold*/ {{{
      Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
      suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast, xlab="", ylab="",xlim=xlims, ylim=ylims,axes=FALSE))
      # /*fend*/ }}}
      #Plot Aperture in Blue /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=(sfa[[i]]), main="Image & Aperture", asp=1, col=hsv(2/3,seq(0,1,length=256)), useRaster=Rast, axes=FALSE, xlab="", ylab="", add=TRUE))
      # /*fend*/ }}}
      #Plot Image in greyscale /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, col=grey.colors(1000), useRaster=Rast,add=TRUE, xlab="", ylab=""))
      # /*fend*/ }}}
      #Overlay Sourcemask in Green /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,(diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1))-(x.pix[i]-data.stamp.lims[i,1]))*arcsec.per.pix,y=(seq(1,(diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1))-(y.pix[i]-data.stamp.lims[i,3]))*arcsec.per.pix, z=log10(smB*image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]), main="", asp=1, useRaster=Rast,add=TRUE, xlab="", ylab="",col=cm.colors(256)))
      # /*fend*/ }}}
      #Plot +ve flux in aperture in Heat Colours /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=log10(apT*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]), main="", asp=1, col=heat.colors(256), useRaster=Rast,add=TRUE, xlab="", ylab=""))
      # /*fend*/ }}}
      #Plot Sources /*fold*/ {{{
      points(x=(x.pix-x.pix[i]+1)*arcsec.per.pix,y=(y.pix-y.pix[i]+1)*arcsec.per.pix, pch=3)
      # /*fend*/ }}}
      #Label with ID /*fold*/ {{{
      label("topleft",lab=cat.id[i],cex=1.5, col='red')
      # /*fend*/ }}}
      #Label Panel /*fold*/ {{{
      label("topleft",lab="(a)",cex=2.5,inset=c(0.1,0.23))
      # /*fend*/ }}}
      #Draw Axes /*fold*/ {{{
      magaxis(frame.plot=TRUE,main="Image & Aperture",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      # /*fend*/ }}}
      #Generate COGs /*fold*/ {{{
      #Raw Image /*fold*/ {{{
      if (psf.weighted) {
        cog<-get.cog(sfa[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      } else {
        cog<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]],centre=c((diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1)/2, (diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1)/2),sample=1E3)
      }
      # /*fend*/ }}}
      #Sky Subtracted only /*fold*/ {{{
      if (psf.weighted) {
        cog.nosky<-get.cog(sfa[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      } else {
        cog.nosky<-get.cog(image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]-skylocal[i],centre=c((diff(range(data.stamp.lims[i,1]:data.stamp.lims[i,2]))+1)/2, (diff(range(data.stamp.lims[i,3]:data.stamp.lims[i,4]))+1)/2),sample=1E3)
      }
      # /*fend*/ }}}
      #Deblended only /*fold*/ {{{
      if (psf.weighted) {
        debl.cog<-get.cog(sfa[[i]]*dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      } else {
        debl.cog<-get.cog(dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      }
      # /*fend*/ }}}
      #Deblended and Sky Subtracted /*fold*/ {{{
      if (psf.weighted) {
        debl.cog.nosky<-get.cog(sfa[[i]]*dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      } else {
        debl.cog.nosky<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3)
      }
      # /*fend*/ }}}
      # /*fend*/ }}}
      #Plot COGs /*fold*/ {{{
      if (magnitudes) {
        #Plot in Magnitude Space /*fold*/ {{{
        #Get y limits /*fold*/ {{{
        suppressWarnings(ylim<-c((-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp)+3, (-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp)-3))
        # /*fend*/ }}}
        #If a limit in ±Inf, plot around median /*fold*/ {{{
        if (!all(is.finite(ylim))) { ylim=median(-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, na.rm=TRUE)+c(-1,3) }
        if (!all(is.finite(ylim))) { warning("Cog flux is always < 0; Using arbitrary plot limits"); ylim=18+c(-1,3) }
        # /*fend*/ }}}
        #Plot Raw COG /*fold*/ {{{
        suppressWarnings(magplot(x=cog$x*arcsec.per.pix, y=-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Magnitude",
                ylim=ylim,type='l',lty=2,main="Curve of Growth"))
        # /*fend*/ }}}
        #Add Lines showing fluxes /*fold*/ {{{
        #Undeblended /*fold*/ {{{
        suppressWarnings(abline(h=-2.5*(log10(sfaflux[i])-log10(ab.vega.flux))+mag.zp, lwd=1, col='orange', lty=1))
        # /*fend*/ }}}
        #Deblended /*fold*/ {{{
        suppressWarnings(abline(h=-2.5*(log10(dfaflux[i])-log10(ab.vega.flux))+mag.zp, lwd=1, col='green', lty=1))
        # /*fend*/ }}}
        #Redraw Raw Cog /*fold*/ {{{
        suppressWarnings(lines(x=cog$x*arcsec.per.pix, y=-2.5*(log10(cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey',lty=2))
        # /*fend*/ }}}
        #Draw Deblended Cog /*fold*/ {{{
        suppressWarnings(lines(x=debl.cog$x*arcsec.per.pix, y=-2.5*(log10(debl.cog$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='black',lty=2))
        # /*fend*/ }}}
        #Draw Sky Subtracted Cog /*fold*/ {{{
        suppressWarnings(lines(x=cog.nosky$x*arcsec.per.pix, y=-2.5*(log10(cog.nosky$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='grey',lty=1))
        # /*fend*/ }}}
        #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
        suppressWarnings(lines(x=debl.cog.nosky$x*arcsec.per.pix, y=-2.5*(log10(debl.cog.nosky$y)-log10(ab.vega.flux))+mag.zp, pch=20, col='black',lty=1))
        # /*fend*/ }}}
        #Draw Legend /*fold*/ {{{
        legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended ApMag","Deblended ApMag"),lty=c(2,2,1,1,1,1),
               col=c('grey','black','grey','black','orange','green'),pch=-1, cex=1)
        # /*fend*/ }}}
        #Note Half Light Radius /*fold*/ {{{
        if (do.sky.est) {
          deproj.debl.cog.nosky<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=c(cat.b[i]/cat.a[i],theta.offset[i]))
          hlr<-which(debl.cog.nosky$y>=max(debl.cog.nosky$y,na.rm=TRUE)/2)
          dhlr<-which(deproj.debl.cog.nosky$y>max(debl.cog.nosky$y,na.rm=TRUE)/2)
          label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2),"\nDeprojected:",round(min(debl.cog.nosky$x[dhlr]),digits=2)))
        } else {
          deproj.debl.cog<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=c(cat.b[i]/cat.a[i],theta.offset[i]))
          hlr<-which(debl.cog$y>=max(debl.cog$y,na.rm=TRUE)/2)
          dhlr<-which(deproj.debl.cog$y>max(debl.cog$y,na.rm=TRUE)/2)
          label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2),"\nDeprojected:",round(min(debl.cog$x[dhlr]),digits=2)))
        }
        # /*fend*/ }}}
        # /*fend*/ }}}
        # /*fend*/ }}}
      } else {
        #Plot in Flux Space /*fold*/ {{{
        #Plot Raw Cog /*fold*/ {{{
        ylim<-range(cog$y)
        #If a limit in ±Inf, plot around median /*fold*/ {{{
        if (!all(is.finite(ylim))) { ylim=median(cog$y, na.rm=TRUE)+c(-1E2,1E2) }
        magplot(x=cog$x*arcsec.per.pix, y=cog$y, pch=20, col='grey', xlab="Radius (arcsec)", ylab="Enclosed Flux",ylim=ylim,main="Curve of Growth",type='l')
        # /*fend*/ }}}
        # /*fend*/ }}}
        #Draw Lines showing fluxes /*fold*/ {{{
        abline(h=dfaflux[i], lwd=1, col='green')
        abline(h=sfaflux[i], lwd=1, col='orange', lty=2)
        # /*fend*/ }}}
        #Redraw Raw Cog /*fold*/ {{{
        lines(x=cog$x*arcsec.per.pix, y=cog$y, pch=20, col='grey',lty=2)
        # /*fend*/ }}}
        #Draw Deblended Cog /*fold*/ {{{
        lines(x=debl.cog$x*arcsec.per.pix, y=debl.cog$y, pch=20, col='black',lty=2)
        # /*fend*/ }}}
        #Draw Sky Subtracted Cog /*fold*/ {{{
        lines(x=cog.nosky$x*arcsec.per.pix, y=cog.nosky$y, pch=20, col='grey',lty=1)
        # /*fend*/ }}}
        #Draw Sky Subtracted & Deblended Cog /*fold*/ {{{
        lines(x=debl.cog.nosky$x*arcsec.per.pix, y=debl.cog.nosky$y, pch=20, col='black',lty=1)
        # /*fend*/ }}}
        legend('bottomright',legend=c("Image COG","Deblended COG","Sky removed COG","Deblended & Sky Rem. COG","Undeblended Flux","Deblended Flux"),lty=c(2,2,1,1,1,1),
               col=c('grey','black','grey','black','orange','green'),pch=-1, cex=1)
        #Note Half Light Radius /*fold*/ {{{
        if (do.sky.est) {
          deproj.debl.cog.nosky<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]-skylocal[i]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=c(cat.b[i]/cat.a[i],theta.offset[i]))
          hlr<-which(debl.cog.nosky$y>=max(cog$y,na.rm=TRUE)/2)
          dhlr<-which(deproj.debl.cog.nosky$y>max(debl.cog$y,na.rm=TRUE)/2)
          label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog.nosky$x[hlr]),digits=2),"\nDeprojected:",round(min(debl.cog.nosky$x[dhlr]),digits=2)))
        } else {
          deproj.debl.cog<-get.cog(dbw[[i]]*(image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]),centre=c(stamplen[i]/2, stamplen[i]/2),sample=1E3,proj=c(cat.b[i]/cat.a[i],theta.offset[i]))
          hlr<-which(debl.cog$y>=max(cog$y,na.rm=TRUE)/2)
          dhlr<-which(deproj.debl.cog$y>max(debl.cog$y,na.rm=TRUE)/2)
          label('topright',lab=paste0("Deblended Half-Light Radius:\nImage:",round(min(debl.cog$x[hlr]),digits=2),"\nDeprojected:",round(min(debl.cog$x[dhlr]),digits=2)))
        }
        # /*fend*/ }}}
        # /*fend*/ }}}
      }
      # /*fend*/ }}}
      #Label Panel /*fold*/ {{{
      label("topleft",lab="(b)",cex=2.5,inset=c(0.1,0.23))
      # /*fend*/ }}}
      #Plot the Deblended Image /*fold*/ {{{
      nc<-length(ap.lims.data.map[i,1]:ap.lims.data.map[i,2])
      nr<-length(ap.lims.data.map[i,3]:ap.lims.data.map[i,4])
      Rast<-ifelse(stamplen[i]>100,TRUE,FALSE)
      suppressWarnings(z<-matrix(magmap(dbw[[i]]*image.env$im[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]],stretch='asinh')$map,ncol=nc,nrow=nr))
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apB, main="Image x Weight Matrix", asp=1, col=grey.colors(256), useRaster=Rast, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims)
      # /*fend*/ }}}
      #Overlay the Aperture /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apT, main="Image x Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=Rast, xlab="", ylab="", axes=FALSE, xlim=xlims, ylim=ylims,add=TRUE))
      # /*fend*/ }}}
      #Draw the projected half-light ellipse {{{
      lines(ellipse(a=min(debl.cog$x[dhlr]),e=1-cat.b[i]/cat.a[i],pa=90-theta.offset[i],xcen=0,ycen=0),col='black',lty=3,lwd=2)
      #}}}
      #Draw the Axes and scalebar /*fold*/ {{{
      magaxis(frame.plot=TRUE,main="Image x Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      # /*fend*/ }}}
      #Label Panel /*fold*/ {{{
      label("topleft",lab="(c)",cex=2.5,inset=c(0.1,0.23))
      # /*fend*/ }}}
      #Plot the Deblend Matrix /*fold*/ {{{
      z=dbw[[i]]
      image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apB, main="Weight Matrix", asp=1, col=grey.colors(256), useRaster=Rast, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims)
      # /*fend*/ }}}
      #Overlay the Aperture /*fold*/ {{{
      suppressWarnings(image(x=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix,y=(seq(1,length(sfa[[i]][,1]))-length(sfa[[i]][,1])/2)*arcsec.per.pix, z=z*apT, main="Weight Matrix", asp=1, col=rev(rainbow(256, start=0,end=2/3)), useRaster=Rast, xlab="", ylab="", axes=FALSE, zlim=c(0,1), xlim=xlims, ylim=ylims,add=TRUE))
      # /*fend*/ }}}
      #Draw the Axes and scalebar /*fold*/ {{{
      magaxis(frame.plot=TRUE,main="Weight Matrix",xlab="Delta RA (arcsec)",ylab="Delta Dec (arcsec)")
      # /*fend*/ }}}
      #Label Panel /*fold*/ {{{
      label("topleft",lab="(d)",cex=2.5,inset=c(0.1,0.23))
      # /*fend*/ }}}
      #Close the file /*fold*/ {{{
      dev.off()
      # /*fend*/ }}}
    }
    #Remove unneeded Arrays /*fold*/ {{{
    if (!make.resid.map) {
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
  if (make.resid.map) {
    if (!is.null(sfabak)) { sfa<-sfabak }
    if (filt.contam) {
      if (!quiet) { cat(paste("Writing Contaminant-subtracted Map to",no.contam.map,"   ")) }
      #Perform Source Subtraction /*fold*/ {{{
      timer=system.time(source.subtraction(image.env$im,sfa,ap.lims.data.map,dfaflux/ApCorr,file.path(path.root,path.work,path.out,no.contam.map),image.env$data.hdr,ba,contams,diagnostic,verbose))
      if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
        message(paste('Contam Subtraction - Done (',round(timer[3], digits=2),'sec )'))
      } else if (!quiet) { cat("   - Done\n") }
    }
    if (!quiet) { cat(paste("Writing Source-subtracted Map to",residual.map,"   ")) }
    # /*fend*/ }}}
    #Perform Source Subtraction /*fold*/ {{{
    timer=system.time(source.subtraction(image.env$im,sfa,ap.lims.data.map,dfaflux/ApCorr,file.path(path.root,path.work,path.out,residual.map),image.env$data.hdr,ba,inside.mask,diagnostic,verbose))
    if (showtime) { cat("   - Done (",round(timer[3],digits=2),"sec )\n")
      message(paste('Source Subtraction - Done (',round(timer[3], digits=2),'sec )'))
    } else if (!quiet) { cat("   - Done\n") }
    # /*fend*/ }}}
  }# /*fend*/ }}}
  #If Subtracting Contaminants, then remove them before output /*fold*/ {{{
  if ((write.tab)&(filt.contam)) {
    cat.id   <-cat.id[which(contams==0)]
    cat.ra   <-cat.ra[which(contams==0)]
    cat.dec  <-cat.dec[which(contams==0)]
    cat.theta<-cat.theta[which(contams==0)]
    theta.offset<-theta.offset[which(contams==0)]
    x.pix    <-x.pix[which(contams==0)]
    y.pix    <-y.pix[which(contams==0)]
    cat.x    <-cat.x[which(contams==0)]
    cat.y    <-cat.y[which(contams==0)]
    cat.a    <-cat.a[which(contams==0)]
    cat.b    <-cat.b[which(contams==0)]
    ssa    <-ssa[which(contams==0)]
    ssfa   <-ssfa[which(contams==0)]
    ssfad  <-ssfad[which(contams==0)]
    qssfad <-qssfad[which(contams==0),]
    spsf   <-spsf[which(contams==0)]
    ssfap  <-ssfap[which(contams==0)]
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
    ssfa2e2<-ssfa2e2[which(contams==0)]
    sdfa2e2<-sdfa2e2[which(contams==0)]
    detecthres<-detecthres[which(contams==0)]
    detecthres.mag<-detecthres.mag[which(contams==0)]
    sfaerr <-sfaerr[which(contams==0)]
    sdfa   <-sdfa[which(contams==0)]
    sdfad  <-sdfad[which(contams==0)]
    qsdfad <-qsdfad[which(contams==0),]
    if (iterate.fluxes) {
      fluxiters<-fluxiters[which(contams==0),]
      erriters<-erriters[which(contams==0),]
      sdfaiters<-sdfaiters[which(contams==0),]
      iterateLost<-iterateLost[which(contams==0)]
    }
    dfaflux<-dfaflux[which(contams==0)]
    deblerr <-deblerr[which(contams==0)]
    dfaerr <-dfaerr[which(contams==0)]
    saturated<-saturated[which(contams==0)]
    pixflux<-pixflux[which(contams==0)]
    ApCorr <-ApCorr[which(contams==0)]
    WtCorr <-WtCorr[which(contams==0)]
    stamplen<-stamplen[which(contams==0)]
    mags   <-mags[which(contams==0)]
    if (length(flux.weight!=1)) { flux.weight<-flux.weight[which(contams==0)] }
    if (ran.cor) { randoms<-randoms[which(contams==0),] }
    if (blank.cor) { blanks<-blanks[which(contams==0),] }
    if (diagnostic) {
      ssfa2 <-ssfa2[which(contams==0)]
      sdfa2 <-sdfa2[which(contams==0)]
      ssfae <-ssfae[which(contams==0)]
      sdfae <-sdfae[which(contams==0)]
      ssfae2<-ssfae2[which(contams==0)]
      sdfae2<-sdfae2[which(contams==0)]
    }
    contams <-contams[which(contams==0)]
  }# /*fend*/ }}}

  #Create Photometry Warning Flags /*fold*/ {{{
  photWarnFlag<-rep("",length(cat.ra))
  #Bad Quartered Photometry /*fold*/ {{{
  Qbad<-(apply(qsdfad,1,'max',na.rm=TRUE)/apply(qsdfad,1,'sum',na.rm=TRUE))>=0.7
  Qbad[which(!is.finite(Qbad))]<-1
  photWarnFlag<-paste0(photWarnFlag,ifelse(Qbad,"Q",""))
  # /*fend*/ }}}
  #Saturation /*fold*/ {{{
  photWarnFlag<-paste0(photWarnFlag,ifelse(saturated,"X",""))
  # /*fend*/ }}}
  #Bad Sky Estimate /*fold*/ {{{
  if (do.sky.est|get.sky.rms) {
    temp<-skyNBinNear
    temp[which(!is.finite(temp))]<-0
    photWarnFlag<-paste0(photWarnFlag,ifelse(temp<=3,"S",""))
  }
  # /*fend*/ }}}
  #Iterative Deblend Warning /*fold*/ {{{
  if (iterate.fluxes) {
    photWarnFlag<-paste0(photWarnFlag,ifelse(iterateLost,"I",""))
  }
  # /*fend*/ }}}
  # /*fend*/ }}}

  #If wanted, output the Results Table /*fold*/ {{{
  if (write.tab) {
    #Output the results table /*fold*/ {{{
    if ((loop.total!=1)&&(length(param.env$tableout.name)!=loop.total)) {
      if (!quiet) { cat(paste('Writing Results Table to ',file.path(path.root,path.work,path.out,paste(sep="",tableout.name,"_file",f,".csv")),'   ')) }
      timer=system.time(write.flux.measurements.table(filename=file.path(path.root,path.work,path.out,paste(sep="",tableout.name,"_file",f,".csv"))) )
    } else {
      if (!quiet) { cat(paste('Writing Results Table to ',file.path(path.root,path.work,path.out,paste(sep="",tableout.name,".csv")),'   ')) }
      timer=system.time(write.flux.measurements.table(filename=file.path(path.root,path.work,path.out,paste(sep="",tableout.name,".csv"))) )
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
     sink(sink.file, type='message')
  }# /*fend*/ }}}
  # /*fend*/ }}}
  # /*fend*/ }}}
  #PART SIX: FINISH /*fold*/ {{{
  #Send Parameters to logfile /*fold*/ {{{
  sink(sink.file, type="output")
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
.executeran.cor<-function() {
cat("Executing ran.cor...\n"); Sys.sleep(2); cat('\n    |   _______   _     _   _     _   _ __      _    |\n    |  |__   __| | |   | | | |   | | |  __ \\   | |   |\n    |     | |    | |___| | | |   | | | |  | |  | |   |\n    |     | |    |  ___  | | |   | | | |  | |  |_|   |\n    |     | |    | |   | | | |___| | | |__| |   _    |\n    |     |_|    |_|   |_| |_______| |_____/   |_|   |\n    |                                                |\n     \\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\\  /\n      \\/  \\/  \\/  \\/  ,------------,  \\/  \\/  \\/  \\/\n         ____        ({ XX      XX })        ____\n        /////|\\      _|   \\\\  //   |_      /|\\\\\\\\\\\n        VVVV | \\____/ { \\  \\\\//  / } \\____/ | VVVV\n        ////_|       ,|V=V=V=V=V=V=|,       |_\\\\\\\\\n     ___\\\\\\\\/_\\~~~~~/_{+^+^+^+^+^+^}_\\~~~~~/_\\////___\n\n'); Sys.sleep(2); cat("... Done, you heartless Jedi Scum\n\n"); Sys.sleep(2)
}
# /*fend*/ }}}
