make_data_mask <-
function(outenv=parent.env(environment()), env=NULL){
#Procedure creates data mask arrays, using input images and
#location data from the catalogue
#Procedure is parallelised to allow scaleability

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
  npos<-length(ra_g)
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
  if (length(which(!(is.finite(x_p) | !(is.finite(y_p)))))!=0) {
    sink(type="message") ; stop(paste("x_p[",which(!(is.finite(x_p))),"] or y_p[",which(!(is.finite(x_p))),"] are NA/NaN/Inf in  Make_Data_Mask()", sep=""))
  }#}}}

  #Set Stamp pixel limits {{{
  message("Setting Stamp Limits")
  #Details {{{
  #Stampwidths: aperture width multiplied by a buffer, plus the
  #PSF WIDHT determined by the desired confidence. Total axis length MUST be odd
  #If stamplen==1, then make the stamplen == smallest nonzero aperture by default
  #This happens when we have a point source and are not convolving with PSF
  #If ALL apertures are point sources, then this will have length 3}}}
  stamplen<-(floor((ceiling(defbuff*a_g*2/asperpix)+ceiling(psf.clip))/2)*2+1)
  #Deal with Point Sources {{{
  if      (all(stamplen==1)) { stamplen[which(stamplen==1)]<-3 }
  else if (any(stamplen==1)) { stamplen[which(stamplen==1)]<-floor(min(a_g[which(a_g>0)],na.rm=TRUE))*2+3 }
  #}}}
  #Calculate Stamp limits in image-pixel space {{{
  stamp_lims<-cbind(x_p-floor(stamplen/2), x_p+floor(stamplen/2), y_p-floor(stamplen/2), y_p+floor(stamplen/2))
  #}}}
  #}}}

  #Check that stamp limits are within image {{{
  message("Removing Stamps that lie outside the image limits")
  catlen<-length(x_g)
  insidemask<-!((stamp_lims[,1]<1)|(stamp_lims[,2]>length(image.env$im[,1]))|(stamp_lims[,3]<1)|(stamp_lims[,4]>length(image.env$im[1,])))
  #}}}

  #Remove any apertures whos stamps would cross the boundary {{{
  if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Single Aperture Stamps are entirely inside the image.") }
  x_p<-x_p[which(insidemask)]
  y_p<-y_p[which(insidemask)]
  x_g<-x_g[which(insidemask)]
  y_g<-y_g[which(insidemask)]
  id_g<-id_g[which(insidemask)]
  ra_g<-ra_g[which(insidemask)]
  dec_g<-dec_g[which(insidemask)]
  theta_g<-theta_g[which(insidemask)]
  a_g<-a_g[which(insidemask)]
  b_g<-b_g[which(insidemask)]
  stamplen<-stamplen[which(insidemask)]
  stamp_lims<-rbind(stamp_lims[which(insidemask),])
  if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
  if (filtcontam) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]
  npos<-length(insidemask)
  #}}}

  #Setup MPI Options {{{
  chunkSize=ceiling(npos/getDoParWorkers())
  mpiopts<-list(chunkSize=chunkSize)
  #}}}

  #Calculate image Stamp limits in image-pixel space {{{
  factor<-ifelse(floor(a_g/asperpix)*5>floor(psffwhm*10), floor(a_g/asperpix)*5, floor(psffwhm*10))
  factor<-ifelse(floor(stamplen/2)    >factor           , floor(stamplen/2)    , factor           )
  im_stamp_lims<-cbind(x_p-factor, x_p+factor, y_p-factor, y_p+factor)
  im_stamp_lims[which(im_stamp_lims<1)]<-1
  im_stamp_lims[which(im_stamp_lims[,4]>length(image.env$im[1,])),4]<-length(image.env$im[1,])
  im_stamp_lims[which(im_stamp_lims[,2]>length(image.env$im[,1])),2]<-length(image.env$im[,1])
  #}}}

  #Create image stamps {{{
  im_mask<-list(NULL)
  for (i in 1:npos) {
    im_mask[[i]]<-image.env$im[im_stamp_lims[i,1]:im_stamp_lims[i,2],im_stamp_lims[i,3]:im_stamp_lims[i,4]]
  }
  #}}}

  #Calculate Stamp limits in image-stamp space {{{
  image_lims<-stamp_lims
  stamp_lims<-stamp_lims-(im_stamp_lims[,c(1,1,3,3)]-1)
  #}}}

  #Notify {{{
  if (verbose) { message(paste("There are",length(x_g),"supplied objects have stamps entirely inside the image (",
                                round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied lie across the edge of the image )")) }
  #}}}

  #Create Mask Stamps {{{
  #Details {{{
  #For each catalogue entry, make a stamp of the mask at that position
  #If the mask is unity everywhere (or is of length 1)
  #we don't need to use the imm mask at all; so skip}}}
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    #Get Mask Limits {{{
    #Get Pixel Locations in Mask {{{
    astr_struc.mask<-read.astr(file.path(pathroot,pathwork,maskmap))
    if (!(all(astr_struc$CRVAL==astr_struc.mask$CRVAL)&
          all(astr_struc$CRPIX==astr_struc.mask$CRPIX)&
          all(astr_struc$CD   ==astr_struc.mask$CD   ))){
      #Get object locations in pixel space {{{
      gamapos<-ad2xy(ra_g,dec_g,astr_struc.mask)
      x_mp<-floor(gamapos[,1])
      y_mp<-floor(gamapos[,2])
      #}}}
      #Calculate Mask limits in mask-pixel space {{{
      imm_stamp_lims<-cbind(x_mp-factor, x_mp+factor, y_mp-factor, y_mp+factor)
      imm_stamp_lims[which(imm_stamp_lims<1)]<-1
      imm_stamp_lims[which(imm_stamp_lims[,4]>length(image.env$imm[1,])),4]<-length(image.env$imm[1,])
      imm_stamp_lims[which(imm_stamp_lims[,2]>length(image.env$imm[,1])),2]<-length(image.env$imm[,1])
      mask_lims<-cbind(x_mp-floor(stamplen/2), x_mp+floor(stamplen/2), y_mp-floor(stamplen/2), y_mp+floor(stamplen/2))
      #}}}
      #}}}
    } else {
      #Use Image pixel locations {{{
      x_mp<-x_p
      y_mp<-y_p
      imm_stamp_lims<-im_stamp_lims
      mask_lims<-image_lims
      #}}}
    }#}}}

    #Create Mask Stamps {{{
    message('Creating Mask Stamps')
    imm_mask<-list(NULL)
    for (i in 1:npos) {
      imm_mask[[i]]<-image.env$imm[imm_stamp_lims[i,1]:imm_stamp_lims[i,2],imm_stamp_lims[i,3]:imm_stamp_lims[i,4]]
    }
    mask_lims<-mask_lims-(imm_stamp_lims[,c(1,1,3,3)]-1)
  } else {
    imm_mask<-1
    mask_lims<-image_lims
    imm_stamp_lims<-im_stamp_lims
  }#}}}

  #Create Error Stamps {{{
  #Details {{{
  #For each catalogue entry, make a stamp of the error map at that position
  #If the error map is a single value everywhere (or is of length 1)
  #we don't need to use the ime mask at all; so skip}}}
  if (((length(image.env$ime)!=1)&(length(unique(image.env$ime))!=1)&(errormap!="NONE"))) {
    #Check Mask Astrometry {{{
    astr_struc.err<-read.astr(file.path(pathroot,pathwork,errormap))
    if (!(all(astr_struc$CRVAL==astr_struc.err$CRVAL)&
          all(astr_struc$CRPIX==astr_struc.err$CRPIX)&
          all(astr_struc$CD   ==astr_struc.err$CD   ))){
      #Get object locations in pixel space {{{
      gamapos<-ad2xy(ra_g,dec_g,astr_struc.err)
      x_ep<-floor(gamapos[,1])
      y_ep<-floor(gamapos[,2])
      #}}}
      #Calculate Mask limits in mask-pixel space {{{
      ime_stamp_lims<-cbind(x_ep-factor, x_ep+factor, y_ep-factor, y_ep+factor)
      ime_stamp_lims[which(ime_stamp_lims<1)]<-1
      ime_stamp_lims[which(ime_stamp_lims[,4]>length(image.env$ime[1,])),4]<-length(image.env$ime[1,])
      ime_stamp_lims[which(ime_stamp_lims[,2]>length(image.env$ime[,1])),2]<-length(image.env$ime[,1])
      error_lims<-cbind(x_ep-floor(stamplen/2), x_ep+floor(stamplen/2), y_ep-floor(stamplen/2), y_ep+floor(stamplen/2))
      #}}}
    } else {
      #Use Image pixel locations {{{
      x_ep<-x_p
      y_ep<-y_p
      #}}}
      #Calculate Mask limits in mask-pixel space {{{
      error_lims<-image_lims
      #}}}
    }
    #}}}

    #Create Error Stamps {{{
    message('Creating Error Stamps')
    ime_mask<-list(NULL)
    for (i in 1:npos) {
        ime_mask[[i]]<-image.env$ime[ime_stamp_lims[i,1]:ime_stamp_lims[i,2],ime_stamp_lims[i,3]:ime_stamp_lims[i,4]]
    }
    #ime_mask<-foreach(xlo=error_lims[,1],xup=error_lims[,2], ylo=error_lims[,3],yup=error_lims[,4], .export=c("imexpixmax","imeypixmax","image.env$ime"), .inorder=TRUE) %do% {
    #  if ((xlo < 1) | (xup > imexpixmax) | (ylo < 1) | (yup > imeypixmax)) {
    #    #Stamp Extends beyond mask image limits - set this mask to zero {{{
    #    array(0, dim=c(length(xlo:xhi),length(ylo:yhi)))
    #    #}}}
    #  } else {
    #    #Stamp is within mask image limits - return mask stamp {{{
    #    image.env$ime[xlo:xup,ylo:yup]
    #    #}}}
    #  }
    #}
    #}}}
  } else if (errormap=="NONE") {
    #Use Image pixel locations {{{
    x_ep<-x_p
    y_ep<-y_p
    #}}}
    #Calculate Mask limits in mask-pixel space {{{
    ime_stamp_lims<-im_stamp_lims
    error_lims<-image_lims
    #}}}
    #Create Error Stamps {{{
    message('Creating Error Stamps')
    ime_mask<-list(NULL)
    for (i in 1:npos) {
        ime_mask[[i]]<-image.env$ime[ime_stamp_lims[i,1]:ime_stamp_lims[i,2],ime_stamp_lims[i,3]:ime_stamp_lims[i,4]]
    }
    #}}}
  } else if (length(image.env$ime)==1) {
    ime_mask<-image.env$ime
    ime_stamp_lims<-im_stamp_lims
    error_lims<-image_lims
  } else if (length(unique(image.env$ime))==1) {
    ime_mask<-unique(image.env$ime)
    ime_stamp_lims<-im_stamp_lims
    error_lims<-image_lims
  }#}}}

  #If we've made multiple masks, Check that all arrays are conformable {{{
  if ((length(imm_mask)!=1)|(length(ime_mask)!=1)) {
    for (i in 1:npos) {
      #Get the masks' dimensions {{{
      dims<-rbind(dim(im_mask[[i]]))
      index<-1
      if (length(imm_mask)!=1) {
        #If the imm_mask exists, get its dimension {{{
        dims<-rbind(dims,dim(imm_mask[[i]]))
        index<-c(index,2)
        #}}}
      }
      if (length(ime_mask)!=1) {
        #If the imm_mask exists, get its dimension{{{
        dims<-rbind(dims,dim(ime_mask[[i]]))
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
            end<-ifelse(im_stamp_lims[i,1]==1,FALSE,TRUE)
          } else if (lo==2) {
            #Mask matrix
            end<-ifelse(imm_stamp_lims[i,1]==1,FALSE,TRUE)
          } else {
            #error matrix
            end<-ifelse(ime_stamp_lims[i,1]==1,FALSE,TRUE)
          }#}}}
          #Find how much they differ {{{
          differ<-diff(c(dims[which.min(dims[,1]),1],dims[which.max(dims[,1]),1]))
          #}}}
          if (end) {
            #Matrix is cutoff at the high end; modify the high matrix {{{
            if (hi==1) {
             #Image matrix
             im_mask[[i]]<-im_mask[[i]][1:(dims[which(index==hi),1]-differ),]
             im_stamp_lims[i,2]<-im_stamp_lims[i,2]-differ
            } else if (hi==2) {
             #Mask matrix
             imm_mask[[i]]<-imm_mask[[i]][1:(dims[which(index==hi),1]-differ),]
             imm_stamp_lims[i,2]<-imm_stamp_lims[i,2]-differ
            } else {
             #error matrix
             ime_mask[[i]]<-ime_mask[[i]][1:(dims[which(index==hi),1]-differ),]
             ime_stamp_lims[i,2]<-ime_stamp_lims[i,2]-differ
            }#}}}
          } else {
          #Matrix is cutoff at the low end {{{
            if (hi==1) {
             #Image matrix
             im_mask[[i]]<-im_mask[[i]][(differ+1):dims[which(index==hi),1],]
             im_stamp_lims[i,1]<-im_stamp_lims[i,1]+differ
            } else if (hi==2) {
             #Mask matrix
             imm_mask[[i]]<-imm_mask[[i]][(differ+1):dims[which(index==hi),1],]
             imm_stamp_lims[i,1]<-imm_stamp_lims[i,1]+differ
            } else {
             #error matrix
             ime_mask[[i]]<-ime_mask[[i]][(differ+1):dims[which(index==hi),1],]
             ime_stamp_lims[i,1]<-ime_stamp_lims[i,1]+differ
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
            end<-ifelse(im_stamp_lims[i,3]==1,FALSE,TRUE)
          } else if (lo==2) {
            #Mask matrix
            end<-ifelse(imm_stamp_lims[i,3]==1,FALSE,TRUE)
          } else {
            #error matrix
            end<-ifelse(ime_stamp_lims[i,3]==1,FALSE,TRUE)
          }#}}}
          #Find how much they differ {{{
          differ<-diff(c(dims[which.min(dims[,2]),2],dims[which.max(dims[,2]),2]))
          #}}}
          if (end) {
            #Matrix is cutoff at the high end; modify the high matrix {{{
            if (hi==1) {
             #Image matrix
             im_mask[[i]]<-im_mask[[i]][,1:(dims[which(index==hi),2]-differ)]
             im_stamp_lims[i,4]<-im_stamp_lims[i,4]-differ
            } else if (hi==2) {
             #Mask matrix
             imm_mask[[i]]<-imm_mask[[i]][,1:(dims[which(index==hi),2]-differ)]
             imm_stamp_lims[i,4]<-imm_stamp_lims[i,4]-differ
            } else {
             #error matrix
             ime_mask[[i]]<-ime_mask[[i]][,1:(dims[which(index==hi),2]-differ)]
             ime_stamp_lims[i,4]<-ime_stamp_lims[i,4]-differ
            }#}}}
          } else {
          #Matrix is cutoff at the low end #{{{
            if (hi==1) {
             #Image matrix
             im_mask[[i]]<-im_mask[[i]][,(differ+1):dims[which(index==hi),2]]
             im_stamp_lims[i,3]<-im_stamp_lims[i,3]+differ
            } else if (hi==2) {
             #Mask matrix
             imm_mask[[i]]<-imm_mask[[i]][,(differ+1):dims[which(index==hi),2]]
             imm_stamp_lims[i,3]<-imm_stamp_lims[i,3]+differ
            } else {
             #error matrix
             ime_mask[[i]]<-ime_mask[[i]][,(differ+1):dims[which(index==hi),2]]
             ime_stamp_lims[i,3]<-ime_stamp_lims[i,3]+differ
            }
          }#}}}
        }
      }#}}}
    }
    #Update Stamp limits in image-stamp space {{{
    stamp_lims<-image_lims-(im_stamp_lims[,c(1,1,3,3)]-1)
    #}}}

    #Check that no stamps have become too small {{{
     insidemask<-!((stamp_lims[,1]<1)|(stamp_lims[,3]<1)|(stamp_lims[,2]>as.numeric(abs(diff(t(im_stamp_lims[,c(1,2)])))+1))|(stamp_lims[,4]>as.numeric(abs(diff(t(im_stamp_lims[,c(3,4)])))+1)))
    if (length(which(!insidemask))!=0) {
      warning(paste("Mismatches in the edges of the input images (image, maskmap, errormap) means",length(insidemask),"stamps have to be discarded"))
      #Remove stamps {{{
      if (length(which(insidemask==TRUE))==0) { sink(type="message") ; stop("No Aperture Stamps are aligned enough to be valid.") }
      im_mask<-im_mask[which(insidemask)]
      if (length(imm_mask)>1) { imm_mask<-imm_mask[which(insidemask)] }
      if (length(ime_mask)>1) { ime_mask<-ime_mask[which(insidemask)] }
      x_p<-x_p[which(insidemask)]
      y_p<-y_p[which(insidemask)]
      x_g<-x_g[which(insidemask)]
      y_g<-y_g[which(insidemask)]
      id_g<-id_g[which(insidemask)]
      ra_g<-ra_g[which(insidemask)]
      dec_g<-dec_g[which(insidemask)]
      theta_g<-theta_g[which(insidemask)]
      a_g<-a_g[which(insidemask)]
      b_g<-b_g[which(insidemask)]
      stamplen<-stamplen[which(insidemask)]
      stamp_lims<-rbind(stamp_lims[which(insidemask),])
      image_lims<-rbind(image_lims[which(insidemask),])
      mask_lims<-rbind(mask_lims[which(insidemask),])
      error_lims<-rbind(error_lims[which(insidemask),])
      im_stamp_lims<-rbind(im_stamp_lims[which(insidemask),])
      imm_stamp_lims<-rbind(imm_stamp_lims[which(insidemask),])
      ime_stamp_lims<-rbind(ime_stamp_lims[which(insidemask),])
      if (length(fluxweight)!=1) { fluxweight<-fluxweight[which(insidemask)] }
      if (filtcontam) { contams<-contams[which(insidemask)] }
      insidemask<-insidemask[which(insidemask)]
      npos<-length(insidemask)
      #}}}
    }
    #}}}
  }
  #}}}

  #Parse Parameter Space {{{
  if (!is.null(env)) { detatch(env) }
  assign("x_p",x_p,envir=outenv)
  assign("y_p",y_p,envir=outenv)
  assign("x_g",x_g,envir=outenv)
  assign("y_g",y_g,envir=outenv)
  assign("id_g",id_g,envir=outenv)
  assign("ra_g",ra_g,envir=outenv)
  assign("dec_g",dec_g,envir=outenv)
  assign("theta_g",theta_g,envir=outenv)
  assign("a_g",a_g,envir=outenv)
  assign("b_g",b_g,envir=outenv)
  assign("stamplen",stamplen,envir=outenv)
  assign("image_lims",image_lims,envir=outenv)
  assign("stamp_lims",stamp_lims,envir=outenv)
  assign("im_stamp_lims",im_stamp_lims,envir=outenv)
  assign("imm_stamp_lims",imm_stamp_lims,envir=outenv)
  assign("mask_lims",mask_lims,envir=outenv)
  assign("imm_mask",imm_mask,envir=outenv)
  assign("ime_mask",ime_mask,envir=outenv)
  assign("fluxweight",fluxweight,envir=outenv)
  if (filtcontam) { assign("contams",contams,envir=outenv) }
  assign("insidemask",insidemask,envir=outenv)
  assign("mpiopts",mpiopts,envir=outenv)
  #}}}

  message('===========END============Make_SA_MASK=============END=================\n')

  #Return array of apertures {{{
  return=im_mask
  #}}}
}
