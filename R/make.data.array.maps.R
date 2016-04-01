make.data.array.maps <-
function(outenv=parent.env(environment()), env=NULL){
#Procedure creates data mask arrays, using input images and
#location data from the catalogue
#Procedure is parallelised to allow scaleability
#
# LIMITS NOMENCLATURE:
#   > ap.lims.data.map is the limits of the aperture stamp in the image.env$im space
#   > ap.lims.data.stamp is the limits of the aperture stamp in the data.stamp space
#   > data.stamp.lims is the limits of the data.stamp in the image.env$im space
#
#   > ap.lims.mask.map is the limits of the aperture stamp in the image.env$imm space
#   > ap.lims.mask.stamp is the limits of the aperture stamp in the mask.stamp space
#   > mask.stamp.lims is the limits of the mask.stamp in the image.env$imm space
#
#   > ap.lims.error.map is the limits of the aperture stamp in the image.env$ime space
#   > ap.lims.error.stamp is the limits of the aperture stamp in the error.stamp space
#   > error.stamp.lims is the limits of the error.stamp in the image.env$ime space
#
#

  if (!quiet) { cat('Make_Data_Masks   ') }
  message('--------------------------Make_Data_Mask---------------------------------')

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

  #Setup sizes {{{
  npos<-length(cat.ra)
  if (length(image.env$imm)!=1) {
    immxpixmax<-length(image.env$imm[,1])
    immypixmax<-length(image.env$imm[1,])
  } else {
    immxpixmax<-1
    immypixmax<-1
  }
  if (length(image.env$ime)!=1) {
    imexpixmax<-length(image.env$ime[,1])
    imeypixmax<-length(image.env$ime[1,])
  } else {
    imexpixmax<-1
    imeypixmax<-1
  }
  #}}}

  #Check to make sure pixel values are finite {{{
  if (length(which(!(is.finite(x.pix) | !(is.finite(y.pix)))))!=0) {
    sink(type="message") ; stop(paste("x.pix[",which(!(is.finite(x.pix))),"] or y.pix[",which(!(is.finite(x.pix))),"] are NA/NaN/Inf in  Make_Data_Mask()", sep=""))
  }#}}}

  #Set Stamp pixel limits {{{
  message("Setting Stamp Limits")
  #Details {{{
  #Stampwidths: aperture width multiplied by a buffer, plus the
  #PSF WIDHT determined by the desired confidence. Total axis length MUST be odd
  #If stamplen==1, then make the stamplen == smallest nonzero aperture by default
  #This happens when we have a point source and are not convolving with PSF
  #If ALL apertures are point sources, then this will have length 5}}}
  if (psf.filt) {
    stamplen<-(floor(ceiling(def.buff*cat.a/arcsec.per.pix)+(ceiling(psf.clip)/2))*2+5)
  } else {
    stamplen<-(floor(ceiling(def.buff*cat.a/arcsec.per.pix))*2+5)
  }
  #Deal with Point Sources {{{
  if      (all(stamplen==5)) { stamplen[which(stamplen==5)]<-5 }
  else if (any(stamplen==5)) { stamplen[which(stamplen==5)]<-floor(min(cat.a[which(cat.a>0)],na.rm=TRUE))*2+5 }
  #}}}
  #Calculate Stamp limits in image-pixel space {{{
  ap.lims.data.stamp<-cbind(x.pix-floor(stamplen/2), x.pix+floor(stamplen/2), y.pix-floor(stamplen/2), y.pix+floor(stamplen/2))
  #}}}
  #}}}

  #Check that stamp limits are within image {{{
  message("Removing Stamps that lie outside the image limits")
  cat.len<-length(cat.x)
  inside.mask<-!((ap.lims.data.stamp[,1]<1)|(ap.lims.data.stamp[,2]>length(image.env$im[,1]))|(ap.lims.data.stamp[,3]<1)|(ap.lims.data.stamp[,4]>length(image.env$im[1,])))
  #}}}

  #Remove any apertures whos stamps would cross the boundary {{{
  if (length(which(inside.mask==TRUE))==0) { sink(type="message") ; stop("No Single Aperture Stamps are entirely inside the image.") }
  x.pix<-x.pix[which(inside.mask)]
  y.pix<-y.pix[which(inside.mask)]
  cat.x<-cat.x[which(inside.mask)]
  cat.y<-cat.y[which(inside.mask)]
  cat.id<-cat.id[which(inside.mask)]
  cat.ra<-cat.ra[which(inside.mask)]
  cat.dec<-cat.dec[which(inside.mask)]
  cat.theta<-cat.theta[which(inside.mask)]
  cat.a<-cat.a[which(inside.mask)]
  cat.b<-cat.b[which(inside.mask)]
  stamplen<-stamplen[which(inside.mask)]
  ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
  if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
  if (filt.contam) { contams<-contams[which(inside.mask)] }
  inside.mask<-inside.mask[which(inside.mask)]
  npos<-length(inside.mask)
  #}}}

  #Setup MPI Options {{{
  chunk.size=ceiling(npos/getDoParWorkers())
  mpi.opts<-list(chunkSize=chunk.size)
  #}}}

  #Calculate image Stamp limits in image-pixel space {{{
  factor<-ifelse(floor(cat.a/arcsec.per.pix)*5>floor(psffwhm*10), floor(cat.a/arcsec.per.pix)*5, floor(psffwhm*10))
  factor<-ifelse(floor(stamplen/2)    >factor           , floor(stamplen/2)    , factor           )
  data.stamp.lims<-cbind(x.pix-factor, x.pix+factor, y.pix-factor, y.pix+factor)
  data.stamp.lims[which(data.stamp.lims<1)]<-1
  data.stamp.lims[which(data.stamp.lims[,4]>length(image.env$im[1,])),4]<-length(image.env$im[1,])
  data.stamp.lims[which(data.stamp.lims[,2]>length(image.env$im[,1])),2]<-length(image.env$im[,1])
  #}}}

  data.stamp<-list(NULL)
  #Determine if cutting up the images is worthwhile {{{
  nchild<-getDoParWorkers()
  #Memory of Cutups {{{
  cutmem<-sum(stamplen^2)*3+sum((factor*2+1)^2)
  #}}}
  if ((cutmem)>(length(image.env$im[,1])*length(image.env$im[1,])*nchild)) {
    cutup<-FALSE
    message("Memory required by cutting up image is more than parsing whole image to children\n>  USE INDIVIDUAL STAMPS: FALSE")
  } else {
    cutup<-TRUE
    message("Memory required by parsing whole image to children is more than cutting into individual stamps\n>  USE INDIVIDUAL STAMPS: TRUE")
  }
  if (!cutup) {
    #total memory from stamps is larger than total memory of moving whole image around; Not worth cutting up {{{
    data.stamp<-NULL
    #}}}
    #Calculate Stamp limits in image-stamp space {{{
    data.stamp.lims<-cbind(rep(1,npos),rep(length(image.env$im[,1]),npos),rep(1,npos),rep(length(image.env$im[1,]),npos))
    ap.lims.data.map<-ap.lims.data.stamp
    ap.lims.data.stamp<-ap.lims.data.stamp
    #}}}
  } else {
    #Total memory of stamps is less that total memory of moving whole image around; Cut up image {{{
    #Create image stamps {{{
    for (i in 1:npos) {
      data.stamp[[i]]<-image.env$im[data.stamp.lims[i,1]:data.stamp.lims[i,2],data.stamp.lims[i,3]:data.stamp.lims[i,4]]
    }
    #}}}
    #}}}
    #Calculate Stamp limits in image-stamp space {{{
    ap.lims.data.map<-ap.lims.data.stamp
    ap.lims.data.stamp<-ap.lims.data.stamp-(data.stamp.lims[,c(1,1,3,3)]-1)
    #}}}
  }
  # }}}

  #Notify {{{
  if (verbose) { message(paste("There are",length(cat.x),"supplied objects have stamps entirely inside the image (",
                                round(((cat.len-length(cat.x))/cat.len)*100, digits=2),"% of supplied lie across the edge of the image )")) }
  #}}}

  #Create Mask Stamps {{{
  #Details {{{
  #For each catalogue entry, make a stamp of the mask at that position
  #If the mask is unity everywhere (or is of length 1)
  #we don't need to use the imm mask at all; so skip}}}
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    #Get Mask Limits {{{
    #Get Pixel Locations in Mask {{{
    if (mask.map=="NONE") {
      #If mask map is NONE, it has been created by the weight.map or error.map
      if (weight.map=="NONE") {
        #Error map
        astr.struc.mask<-read.astrometry(file.path(path.root,path.work,error.map))
      } else {
        #Weight map
        astr.struc.mask<-read.astrometry(file.path(path.root,path.work,weight.map))
      }
    } else {
      astr.struc.mask<-read.astrometry(file.path(path.root,path.work,mask.map))
    }

    if (!(all(astr.struc$CRVAL==astr.struc.mask$CRVAL,na.rm=TRUE)&
          all(astr.struc$CRPIX==astr.struc.mask$CRPIX,na.rm=TRUE)&
          all(astr.struc$CD   ==astr.struc.mask$CD   ,na.rm=TRUE))){
      #Get object locations in pixel space {{{
      gama.pos<-ad.to.xy(cat.ra,cat.dec,astr.struc.mask)
      mask.x.pix<-floor(gama.pos[,1])
      mask.y.pix<-floor(gama.pos[,2])
      #}}}
    } else {
      #Use Image pixel locations {{{
      mask.x.pix<-x.pix
      mask.y.pix<-y.pix
      #}}}
    }
    #}}}
    #}}}

    #Calculate Mask limits in mask-pixel space {{{
    mask.stamp.lims<-cbind(mask.x.pix-factor, mask.x.pix+factor, mask.y.pix-factor, mask.y.pix+factor)
    mask.stamp.lims[which(mask.stamp.lims<1)]<-1
    mask.stamp.lims[which(mask.stamp.lims[,4]>length(image.env$imm[1,])),4]<-length(image.env$imm[1,])
    mask.stamp.lims[which(mask.stamp.lims[,2]>length(image.env$imm[,1])),2]<-length(image.env$imm[,1])
    ap.lims.mask.map<-cbind(mask.x.pix-floor(stamplen/2), mask.x.pix+floor(stamplen/2), mask.y.pix-floor(stamplen/2), mask.y.pix+floor(stamplen/2))
    #}}}

    #Check that stamp limits are within image {{{
    cat.len<-length(cat.x)
    inside.mask<-!((ap.lims.mask.map[,1]<1)|(ap.lims.mask.map[,2]>length(image.env$imm[,1]))|(ap.lims.mask.map[,3]<1)|(ap.lims.mask.map[,4]>length(image.env$imm[1,])))
    #}}}

    #Remove any apertures whos stamps would cross the boundary {{{
    if (any(!inside.mask)) {
      message("Removing Stamps that lie outside the mask image limits")
      if (length(which(inside.mask==TRUE))==0) { sink(type="message") ; stop("No Single Aperture Stamps are entirely inside the image.") }
      factor<-factor[which(inside.mask)]
      x.pix<-x.pix[which(inside.mask)]
      y.pix<-y.pix[which(inside.mask)]
      mask.x.pix<-mask.x.pix[which(inside.mask)]
      mask.y.pix<-mask.y.pix[which(inside.mask)]
      cat.x<-cat.x[which(inside.mask)]
      cat.y<-cat.y[which(inside.mask)]
      cat.id<-cat.id[which(inside.mask)]
      cat.ra<-cat.ra[which(inside.mask)]
      cat.dec<-cat.dec[which(inside.mask)]
      cat.theta<-cat.theta[which(inside.mask)]
      cat.a<-cat.a[which(inside.mask)]
      cat.b<-cat.b[which(inside.mask)]
      stamplen<-stamplen[which(inside.mask)]
      if (length(data.stamp )>1) { data.stamp<-data.stamp[which(inside.mask)] }
      ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
      ap.lims.data.map<-rbind(ap.lims.data.map[which(inside.mask),])
      data.stamp.lims<-rbind(data.stamp.lims[which(inside.mask),])
      ap.lims.mask.map<-rbind(ap.lims.mask.map[which(inside.mask),])
      mask.stamp.lims<-rbind(mask.stamp.lims[which(inside.mask),])
      if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
      if (filt.contam) { contams<-contams[which(inside.mask)] }
      inside.mask<-inside.mask[which(inside.mask)]
      npos<-length(inside.mask)
    }
    #}}}

    #Setup MPI Options {{{
    chunk.size=ceiling(npos/getDoParWorkers())
    mpi.opts<-list(chunkSize=chunk.size)
    #}}}

    #Create Mask Stamps {{{
    message('Creating Mask Stamps')
    mask.stamp<-list(NULL)
    if (!cutup) {
      mask.stamp<-1
      ap.lims.mask.stamp<-ap.lims.mask.map
      mask.stamp.lims<-cbind(rep(1,npos),rep(length(image.env$imm[,1]),npos),rep(1,npos),rep(length(image.env$imm[1,]),npos))
    } else {
      for (i in 1:npos) {
        mask.stamp[[i]]<-image.env$imm[mask.stamp.lims[i,1]:mask.stamp.lims[i,2],mask.stamp.lims[i,3]:mask.stamp.lims[i,4]]
      }
      ap.lims.mask.stamp<-ap.lims.mask.map-(mask.stamp.lims[,c(1,1,3,3)]-1)
    }
  } else {
    mask.x.pix<-x.pix
    mask.y.pix<-y.pix
    mask.stamp<-1
    image.env$imm<-1
    ap.lims.mask.map<-ap.lims.data.map
    ap.lims.mask.stamp<-ap.lims.data.stamp
    mask.stamp.lims<-data.stamp.lims
  }#}}}
  #}}}

  #Create Error Stamps {{{
  #Details {{{
  #For each catalogue entry, make a stamp of the error map at that position
  #If the error map is a single value everywhere (or is of length 1)
  #we don't need to use the ime mask at all; so skip}}}
  if ((length(image.env$ime)!=1)&(any(image.env$ime!=image.env$ime[1]))&(file.exists(file.path(path.root,path.work,error.map)))) {
    #Check Errormap Astrometry {{{
    astr.struc.err<-read.astrometry(file.path(path.root,path.work,error.map))
    if (!(all(astr.struc$CRVAL==astr.struc.err$CRVAL,na.rm=TRUE)&
          all(astr.struc$CRPIX==astr.struc.err$CRPIX,na.rm=TRUE)&
          all(astr.struc$CD   ==astr.struc.err$CD   ,na.rm=TRUE))){
      #Get object locations in pixel space {{{
      gama.pos<-ad.to.xy(cat.ra,cat.dec,astr.struc.err)
      error.x.pix<-floor(gama.pos[,1])
      error.y.pix<-floor(gama.pos[,2])
      #}}}
    } else {
      #Use Image pixel locations {{{
      error.x.pix<-x.pix
      error.y.pix<-y.pix
      #}}}
    }
    #}}}
    #Calculate Mask limits in mask-pixel space {{{
    error.stamp.lims<-cbind(error.x.pix-factor, error.x.pix+factor, error.y.pix-factor, error.y.pix+factor)
    error.stamp.lims[which(error.stamp.lims<1)]<-1
    error.stamp.lims[which(error.stamp.lims[,4]>length(image.env$ime[1,])),4]<-length(image.env$ime[1,])
    error.stamp.lims[which(error.stamp.lims[,2]>length(image.env$ime[,1])),2]<-length(image.env$ime[,1])
    ap.lims.error.map<-cbind(error.x.pix-floor(stamplen/2), error.x.pix+floor(stamplen/2), error.y.pix-floor(stamplen/2), error.y.pix+floor(stamplen/2))
    #}}}

    #Check that stamp limits are within image {{{
    cat.len<-length(cat.x)
    inside.mask<-!((ap.lims.error.map[,1]<1)|(ap.lims.error.map[,2]>length(image.env$ime[,1]))|(ap.lims.error.map[,3]<1)|(ap.lims.error.map[,4]>length(image.env$ime[1,])))
    #}}}

    #Remove any apertures whos stamps would cross the boundary {{{
    if (any(!inside.mask)) {
      message("Removing Stamps that lie outside the error image limits")
      if (length(which(inside.mask==TRUE))==0) { sink(type="message") ; stop("No Single Aperture Stamps are entirely inside the image.") }
      factor<-factor[which(inside.mask)]
      x.pix<-x.pix[which(inside.mask)]
      y.pix<-y.pix[which(inside.mask)]
      mask.x.pix<-mask.x.pix[which(inside.mask)]
      mask.y.pix<-mask.y.pix[which(inside.mask)]
      error.x.pix<-error.x.pix[which(inside.mask)]
      error.y.pix<-error.y.pix[which(inside.mask)]
      cat.x<-cat.x[which(inside.mask)]
      cat.y<-cat.y[which(inside.mask)]
      cat.id<-cat.id[which(inside.mask)]
      cat.ra<-cat.ra[which(inside.mask)]
      cat.dec<-cat.dec[which(inside.mask)]
      cat.theta<-cat.theta[which(inside.mask)]
      cat.a<-cat.a[which(inside.mask)]
      cat.b<-cat.b[which(inside.mask)]
      if (length(data.stamp )>1) { data.stamp<-data.stamp[which(inside.mask)] }
      if (length(mask.stamp)>1) { mask.stamp<-mask.stamp[which(inside.mask)] }
      stamplen<-stamplen[which(inside.mask)]
      ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
      ap.lims.data.map<-rbind(ap.lims.data.map[which(inside.mask),])
      data.stamp.lims<-rbind(data.stamp.lims[which(inside.mask),])
      ap.lims.mask.stamp<-rbind(ap.lims.mask.stamp[which(inside.mask),])
      ap.lims.mask.map<-rbind(ap.lims.mask.map[which(inside.mask),])
      mask.stamp.lims<-rbind(mask.stamp.lims[which(inside.mask),])
      ap.lims.error.map<-rbind(ap.lims.error.map[which(inside.mask),])
      error.stamp.lims<-rbind(error.stamp.lims[which(inside.mask),])
      if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
      if (filt.contam) { contams<-contams[which(inside.mask)] }
      inside.mask<-inside.mask[which(inside.mask)]
      npos<-length(inside.mask)
    }
    #}}}

    #Setup MPI Options {{{
    chunk.size=ceiling(npos/getDoParWorkers())
    mpi.opts<-list(chunkSize=chunk.size)
    #}}}


    #Create Error Stamps {{{
    message('Creating Error Stamps')
    error.stamp<-list(NULL)
    if (!cutup) {
      error.stamp<-NULL
      ap.lims.error.stamp<-ap.lims.error.map
      error.stamp.lims<-cbind(rep(1,npos),rep(length(image.env$ime[,1]),npos),rep(1,npos),rep(length(image.env$ime[1,]),npos))
    } else {
      for (i in 1:npos) {
          error.stamp[[i]]<-image.env$ime[error.stamp.lims[i,1]:error.stamp.lims[i,2],error.stamp.lims[i,3]:error.stamp.lims[i,4]]
      }
      ap.lims.error.stamp<-ap.lims.error.map-(error.stamp.lims[,c(1,1,3,3)]-1)
    }
  } else if (error.map=="NONE") {
    #Use Image pixel locations {{{
    error.x.pix<-x.pix
    error.y.pix<-y.pix
    #}}}
    #Calculate Mask limits in mask-pixel space {{{
    error.stamp.lims<-data.stamp.lims
    ap.lims.error.map<-ap.lims.data.map
    ap.lims.error.stamp<-ap.lims.data.stamp
    #}}}
    #Create Error Stamps {{{
    message('Creating Error Stamps')
    error.stamp<-list(NULL)
    if (!cutup) {
      error.stamp<-NULL
    } else {
      for (i in 1:npos) {
          error.stamp[[i]]<-image.env$ime[error.stamp.lims[i,1]:error.stamp.lims[i,2],error.stamp.lims[i,3]:error.stamp.lims[i,4]]
      }
    }
    #}}}
  } else if (length(image.env$ime)==1) {
    error.x.pix<-x.pix
    error.y.pix<-y.pix
    error.stamp<-image.env$ime
    error.stamp.lims<-data.stamp.lims
    ap.lims.error.map<-ap.lims.data.map
    ap.lims.error.stamp<-ap.lims.data.stamp
  } else if (!any(image.env$ime!=image.env$ime[1])) {
    error.x.pix<-x.pix
    error.y.pix<-y.pix
    error.stamp<-image.env$ime[1]
    error.stamp.lims<-data.stamp.lims
    ap.lims.error.map<-ap.lims.data.map
    ap.lims.error.stamp<-ap.lims.data.stamp
  }#}}}
  #}}}

  #If we've made multiple masks, Check that all arrays are conformable {{{
  check<-FALSE
  if (cutup) {
    if ((length(mask.stamp)!=1)|(length(error.stamp)!=1)) { check<-TRUE }
  } else {
    l1<-length(image.env$im)!=1
    l2<-length(image.env$imm)!=1
    l3<-length(image.env$ime)!=1
    d1<-dim(image.env$im)
    d2<-dim(image.env$imm)
    d3<-dim(image.env$ime)
    if (l1&l2&(any(d1!=d2)|any(x.pix!=mask.x.pix)|any(y.pix!=mask.y.pix))) { check<-TRUE }
    if (l1&l3&(any(d1!=d3)|any(x.pix!=error.x.pix)|any(y.pix!=error.y.pix))) { check<-TRUE }
    if (l2&l3&(any(d2!=d3)|any(mask.x.pix!=error.x.pix)|any(mask.y.pix!=error.y.pix))) { check<-TRUE }
  }

  if (check) {
    for (i in 1:npos) {
      if (!cutup) {
        im.mk<-matrix(0,nrow=length(data.stamp.lims[i,1]:data.stamp.lims[i,2]),ncol=length(data.stamp.lims[i,3]:data.stamp.lims[i,4]))
      } else {
        im.mk<-data.stamp[[i]]
      }
      #Get the masks' dimensions {{{
      dims<-rbind(dim(im.mk))
      index<-1
      if (length(image.env$imm)!=1) {
        if (!cutup) {
          imm.mk<-matrix(0,nrow=length(mask.stamp.lims[i,1]:mask.stamp.lims[i,2]),ncol=length(mask.stamp.lims[i,3]:mask.stamp.lims[i,4]))
        } else {
          imm.mk<-mask.stamp[[i]]
        }
        #If the mask.stamp exists, get its dimension {{{
        dims<-rbind(dims,dim(imm.mk))
        index<-c(index,2)
        #}}}
      }
      if (length(image.env$ime)!=1) {
        if (!cutup) {
          ime.mk<-matrix(0,nrow=length(error.stamp.lims[i,1]:error.stamp.lims[i,2]),ncol=length(error.stamp.lims[i,3]:error.stamp.lims[i,4]))
        } else {
          ime.mk<-error.stamp[[i]]
        }
        #If the error.stamp exists, get its dimension{{{
        dims<-rbind(dims,dim(ime.mk))
        index<-c(index,3)
        #}}}
      }#}}}
      #If any x-dimensions are different {{{
      if (any(dims[,1]!=dims[1,1])) {
        #Determine how many dimensions are different {{{
        len=length(which(dims[,1]!=min(dims[,1])))
        #}}}
        #For the number of different dimensions:
        for (j in 1:len) {
          #Get the matrix with the smallest dimension {{{
          lo<-index[which.min(dims[,1])]
          #}}}
          #Get the matrix with the largest dimension {{{
          hi<-index[which.max(dims[,1])]
          #}}}
          #Determine if the matrix has been cutoff at the low or high end {{{
          if (lo==1) {
            #Image matrix
            end<-ifelse(data.stamp.lims[i,1]==1,FALSE,TRUE)
          } else if (lo==2) {
            #Mask matrix
            end<-ifelse(mask.stamp.lims[i,1]==1,FALSE,TRUE)
          } else {
            #error matrix
            end<-ifelse(error.stamp.lims[i,1]==1,FALSE,TRUE)
          }#}}}
          #Find how much they differ {{{
          differ<-diff(c(dims[which.min(dims[,1]),1],dims[which.max(dims[,1]),1]))
          #}}}
          if (end) {
            #Matrix is cutoff at the high end; modify the high matrix {{{
            if (hi==1) {
             #Image matrix
             im.mk<-im.mk[1:(dims[which(index==hi),1]-differ),]
             if (cutup) { data.stamp[[i]]<-data.stamp[[i]][1:(dims[which(index==hi),1]-differ),] }
             data.stamp.lims[i,2]<-data.stamp.lims[i,2]-differ
             dims[which(index==hi),]<-dim(im.mk)
            } else if (hi==2) {
             #Mask matrix
             imm.mk<-imm.mk[1:(dims[which(index==hi),1]-differ),]
             if (cutup) { mask.stamp[[i]]<-mask.stamp[[i]][1:(dims[which(index==hi),1]-differ),] }
             mask.stamp.lims[i,2]<-mask.stamp.lims[i,2]-differ
             dims[which(index==hi),]<-dim(imm.mk)
            } else {
             #error matrix
             ime.mk<-ime.mk[1:(dims[which(index==hi),1]-differ),]
             if (cutup) { error.stamp[[i]]<-error.stamp[[i]][1:(dims[which(index==hi),1]-differ),] }
             error.stamp.lims[i,2]<-error.stamp.lims[i,2]-differ
             dims[which(index==hi),]<-dim(ime.mk)
            }#}}}
          } else {
          #Matrix is cutoff at the low end {{{
            if (hi==1) {
             #Image matrix
             im.mk<-im.mk[(differ+1):dims[which(index==hi),1],]
             if (cutup) { data.stamp[[i]]<-data.stamp[[i]][(differ+1):dims[which(index==hi),1],] }
             data.stamp.lims[i,1]<-data.stamp.lims[i,1]+differ
             dims[which(index==hi),]<-dim(im.mk)
            } else if (hi==2) {
             #Mask matrix
             imm.mk<-imm.mk[(differ+1):dims[which(index==hi),1],]
             if (cutup) { mask.stamp[[i]]<-mask.stamp[[i]][(differ+1):dims[which(index==hi),1],] }
             mask.stamp.lims[i,1]<-mask.stamp.lims[i,1]+differ
             dims[which(index==hi),]<-dim(imm.mk)
            } else {
             #error matrix
             ime.mk<-ime.mk[(differ+1):dims[which(index==hi),1],]
             if (cutup) { error.stamp[[i]]<-error.stamp[[i]][(differ+1):dims[which(index==hi),1],] }
             error.stamp.lims[i,1]<-error.stamp.lims[i,1]+differ
             dims[which(index==hi),]<-dim(ime.mk)
            }
          }#}}}
        }
      }#}}}
      #If any y-dimensions are different {{{
      if (any(dims[,2]!=dims[1,2])) {
        #Determine how many dimensions are different {{{
        len=length(which(dims[,2]!=min(dims[,2])))
        #}}}
        #For the number of different dimensions:
        for (j in 1:len) {
          #Get the matrix with the smallest dimension #{{{
          lo<-index[which.min(dims[,2])]
          #}}}
          #Get the matrix with the largest dimension #{{{
          hi<-index[which.max(dims[,2])]
          #}}}
          #Determine if the matrix has been cutoff at the low or high end {{{
          if (lo==1) {
            #Image matrix
            end<-ifelse(data.stamp.lims[i,3]==1,FALSE,TRUE)
          } else if (lo==2) {
            #Mask matrix
            end<-ifelse(mask.stamp.lims[i,3]==1,FALSE,TRUE)
          } else {
            #error matrix
            end<-ifelse(error.stamp.lims[i,3]==1,FALSE,TRUE)
          }#}}}
          #Find how much they differ {{{
          differ<-diff(c(dims[which.min(dims[,2]),2],dims[which.max(dims[,2]),2]))
          #}}}
          if (end) {
            #Matrix is cutoff at the high end; modify the high matrix {{{
            if (hi==1) {
             #Image matrix
             im.mk<-im.mk[,1:(dims[which(index==hi),2]-differ)]
             if (cutup) { data.stamp[[i]]<-data.stamp[[i]][,1:(dims[which(index==hi),2]-differ)] }
             data.stamp.lims[i,4]<-data.stamp.lims[i,4]-differ
             dims[which(index==hi),]<-dim(im.mk)
            } else if (hi==2) {
             #Mask matrix
             imm.mk<-imm.mk[,1:(dims[which(index==hi),2]-differ)]
             if (cutup) { mask.stamp[[i]]<-mask.stamp[[i]][,1:(dims[which(index==hi),2]-differ)] }
             mask.stamp.lims[i,4]<-mask.stamp.lims[i,4]-differ
             dims[which(index==hi),]<-dim(imm.mk)
            } else {
             #error matrix
             ime.mk<-ime.mk[,1:(dims[which(index==hi),2]-differ)]
             if (cutup) { error.stamp[[i]]<-error.stamp[[i]][,1:(dims[which(index==hi),2]-differ)] }
             error.stamp.lims[i,4]<-error.stamp.lims[i,4]-differ
             dims[which(index==hi),]<-dim(ime.mk)
            }#}}}
          } else {
          #Matrix is cutoff at the low end #{{{
            if (hi==1) {
             #Image matrix
             im.mk<-im.mk[,(differ+1):dims[which(index==hi),2]]
             if (cutup) { data.stamp[[i]]<-data.stamp[[i]][,(differ+1):dims[which(index==hi),2]] }
             data.stamp.lims[i,3]<-data.stamp.lims[i,3]+differ
             dims[which(index==hi),]<-dim(im.mk)
            } else if (hi==2) {
             #Mask matrix
             imm.mk<-imm.mk[,(differ+1):dims[which(index==hi),2]]
             if (cutup) { mask.stamp[[i]]<-mask.stamp[[i]][,(differ+1):dims[which(index==hi),2]] }
             mask.stamp.lims[i,3]<-mask.stamp.lims[i,3]+differ
             dims[which(index==hi),]<-dim(imm.mk)
            } else {
             #error matrix
             ime.mk<-ime.mk[,(differ+1):dims[which(index==hi),2]]
             if (cutup) { error.stamp[[i]]<-error.stamp[[i]][,(differ+1):dims[which(index==hi),2]] }
             error.stamp.lims[i,3]<-error.stamp.lims[i,3]+differ
             dims[which(index==hi),]<-dim(ime.mk)
            }
          }#}}}
        }
      }#}}}
      if (!cutup){
        #If not cutting up, will be the same for all apertures; Change all and break
        data.stamp.lims<-cbind(rep(data.stamp.lims[i,1],npos),rep(data.stamp.lims[i,2],npos),rep(data.stamp.lims[i,3],npos),rep(data.stamp.lims[i,4],npos))
        mask.stamp.lims<-cbind(rep(mask.stamp.lims[i,1],npos),rep(mask.stamp.lims[i,2],npos),rep(mask.stamp.lims[i,3],npos),rep(mask.stamp.lims[i,4],npos))
        error.stamp.lims<-cbind(rep(error.stamp.lims[i,1],npos),rep(error.stamp.lims[i,2],npos),rep(error.stamp.lims[i,3],npos),rep(error.stamp.lims[i,4],npos))
        break
      }
    }
    #Update Stamp limits in image-stamp space {{{
    ap.lims.data.stamp<-ap.lims.data.map-(data.stamp.lims[,c(1,1,3,3)]-1)
    if (length(image.env$imm)!=1) {
      ap.lims.mask.stamp<-ap.lims.mask.map-(mask.stamp.lims[,c(1,1,3,3)]-1)
    } else {
      mask.stamp.lims<-data.stamp.lims
      ap.lims.mask.stamp<-ap.lims.data.stamp
    }
    if (length(image.env$ime)!=1) {
      ap.lims.error.stamp<-ap.lims.error.map-(error.stamp.lims[,c(1,1,3,3)]-1)
    } else {
      error.stamp.lims<-data.stamp.lims
      ap.lims.error.stamp<-ap.lims.data.stamp
    }
    #}}}

    #Check that no stamps have become too small {{{
     inside.mask<-!((ap.lims.data.stamp[,1]<1)|(ap.lims.data.stamp[,3]<1)|(ap.lims.data.stamp[,2]>as.numeric(abs(diff(t(data.stamp.lims[,c(1,2)])))+1))|(ap.lims.data.stamp[,4]>as.numeric(abs(diff(t(data.stamp.lims[,c(3,4)])))+1)))
    if (length(which(!inside.mask))!=0) {
      warning(paste("Mismatches in the edges of the input images (image, mask.map, error.map) means",length(which(!inside.mask)),"stamps have to be discarded"))
      #Remove stamps {{{
      data.stamp<-data.stamp[which(inside.mask)]
      if (length(mask.stamp)>1) { mask.stamp<-mask.stamp[which(inside.mask)] }
      if (length(error.stamp)>1) { error.stamp<-error.stamp[which(inside.mask)] }
      x.pix<-x.pix[which(inside.mask)]
      y.pix<-y.pix[which(inside.mask)]
      mask.x.pix<-mask.x.pix[which(inside.mask)]
      mask.y.pix<-mask.y.pix[which(inside.mask)]
      error.x.pix<-error.x.pix[which(inside.mask)]
      error.y.pix<-error.y.pix[which(inside.mask)]
      cat.x<-cat.x[which(inside.mask)]
      cat.y<-cat.y[which(inside.mask)]
      cat.id<-cat.id[which(inside.mask)]
      cat.ra<-cat.ra[which(inside.mask)]
      cat.dec<-cat.dec[which(inside.mask)]
      cat.theta<-cat.theta[which(inside.mask)]
      cat.a<-cat.a[which(inside.mask)]
      cat.b<-cat.b[which(inside.mask)]
      stamplen<-stamplen[which(inside.mask)]
      if (length(data.stamp )>1) { data.stamp<-data.stamp[which(inside.mask)] }
      if (length(mask.stamp)>1) { mask.stamp<-mask.stamp[which(inside.mask)] }
      if (length(error.stamp)>1) { error.stamp<-error.stamp[which(inside.mask)] }
      ap.lims.data.stamp<-rbind(ap.lims.data.stamp[which(inside.mask),])
      ap.lims.mask.stamp<-rbind(ap.lims.mask.stamp[which(inside.mask),])
      ap.lims.error.stamp<-rbind(ap.lims.error.stamp[which(inside.mask),])
      ap.lims.data.map<-rbind(ap.lims.data.map[which(inside.mask),])
      ap.lims.mask.map<-rbind(ap.lims.mask.map[which(inside.mask),])
      ap.lims.error.map<-rbind(ap.lims.error.map[which(inside.mask),])
      data.stamp.lims<-rbind(data.stamp.lims[which(inside.mask),])
      mask.stamp.lims<-rbind(mask.stamp.lims[which(inside.mask),])
      error.stamp.lims<-rbind(error.stamp.lims[which(inside.mask),])
      if (length(flux.weight)!=1) { flux.weight<-flux.weight[which(inside.mask)] }
      if (filt.contam) { contams<-contams[which(inside.mask)] }
      inside.mask<-inside.mask[which(inside.mask)]
      npos<-length(inside.mask)
      #}}}
    }
    #}}}
  }
  #}}}

  #Parse Parameter Space {{{
  if (!is.null(env)) { detach(env) }
  assign("cutup",cutup,envir=outenv)
  assign("x.pix",x.pix,envir=outenv)
  assign("y.pix",y.pix,envir=outenv)
  assign("cat.x",cat.x,envir=outenv)
  assign("cat.y",cat.y,envir=outenv)
  assign("cat.id",cat.id,envir=outenv)
  assign("cat.ra",cat.ra,envir=outenv)
  assign("cat.dec",cat.dec,envir=outenv)
  assign("cat.theta",cat.theta,envir=outenv)
  assign("cat.a",cat.a,envir=outenv)
  assign("cat.b",cat.b,envir=outenv)
  assign("stamplen",stamplen,envir=outenv)
  assign("ap.lims.data.map",ap.lims.data.map,envir=outenv)
  assign("ap.lims.mask.map",ap.lims.mask.map,envir=outenv)
  assign("ap.lims.error.map",ap.lims.error.map,envir=outenv)
  assign("ap.lims.data.stamp",ap.lims.data.stamp,envir=outenv)
  assign("ap.lims.mask.stamp",ap.lims.mask.stamp,envir=outenv)
  assign("ap.lims.error.stamp",ap.lims.error.stamp,envir=outenv)
  assign("data.stamp.lims",data.stamp.lims,envir=outenv)
  assign("mask.stamp.lims",mask.stamp.lims,envir=outenv)
  assign("error.stamp.lims",error.stamp.lims,envir=outenv)
  assign("mask.stamp",mask.stamp,envir=outenv)
  assign("error.stamp",error.stamp,envir=outenv)
  assign("flux.weight",flux.weight,envir=outenv)
  if (filt.contam) { assign("contams",contams,envir=outenv) }
  assign("inside.mask",inside.mask,envir=outenv)
  assign("mpi.opts",mpi.opts,envir=outenv)
  #}}}

  message('===========END============Make_SA_MASK=============END=================\n')

  #Return array of apertures {{{
  return=data.stamp
  #}}}
}
