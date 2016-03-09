make_sfa_mask <-
function(outenv=parent.env(environment()), sa_mask,fluxweightin=NULL, immask=NULL, env=NULL,subs=NULL) {
#Make the Single Filtered Aperture mask

  #Print appropriate banner  {{{
  if (length(fluxweightin) > 0 ) {
    if (!quiet) { cat("Make_WSFA_Mask     ")  }
    message('-------------------------Make_WSFA_Mask--------------------------------')
  } else {
    if (!quiet) { cat("Make_SFA_Mask     ")  }
    message('-------------------------Make_SFA_Mask---------------------------------')
  }#}}}

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

  #Setup Sizes {{{
  if (is.null(subs)) {
    npos<-length(id_g)
    subs<-1:npos
  } else {
    npos<-length(subs)
    message(paste("-> There are",npos,"apertures in this subset\n"))
  }
  #}}}

  #Set Fluxweights {{{
  #If no fluxweight - all weights == 1 {{{
  if (length(fluxweightin) == 0) { fluxweight<-array(1, dim=c(npos))
  #}}}
  #If single value suplied, all weights == value {{{
  } else if (length(fluxweightin) == 1) { fluxweight<-array(fluxweightin, dim=c(npos))
  #}}}
  #If vector of values supplied, weights == input {{{
  } else { fluxweight<-fluxweightin }
  #}}}

  #Check that number of Fluxweights == number of subs stamps or number of total stamps {{{
  if (length(fluxweight) != npos) {
    if (length(fluxweight) == length(id_g)) {
      fluxweight<-fluxweight[subs]
    } else {
      sink(type="message") ; stop("Bad FluxWeightIn Array")
    }
  }
  #}}}

  #If we have input fluxweights, Check their type {{{
  if (length(fluxweightin)==0) {
    message("No input Fluxweights; no rescale is required")
  } else if (weightType=="mag") {
    #Fluxweights are Previously determined object Magnitudes {{{
    #Convert back to Flux {{{
    fluxweight<-10^((8.9-fluxweight)/2.5)
    #}}}
    #Get Range {{{
    wgtRange<-range(fluxweight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(fluxweight))) {
      #Bad values. Things with bad fluxweights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(fluxweight))
      fluxweight[ind]<-NA
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(fluxweight<=0)
      fluxweight[ind]<-NA
      #}}}
    }#}}}
    #Adjust Fluxweights by aperture integral
    fluxweight<-foreach(samask=sa_mask[subs],fluxwgt=fluxweight, .inorder=TRUE, .options.mpi=mpiopts, .combine='c') %dopar%{ fluxwgt/sum(samask)  }
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(fluxweight,na.rm=TRUE)<=0) {
      message("WARNING: No fluxweights are > 0! Flux weighting has removed all objects from the image.")
      fluxweight<-rep(0,length(subs))
    } else if (length(which(fluxweight>0))<2) {
      message("WARNING: only one fluxweight is > 0! Flux weighting has removed all other objects from the image.")
      fluxweight[which(fluxweight<=0)]<-NA
      fluxweight<-magmap(fluxweight,lo=0,hi=max(fluxweight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
    } else {
      #Map Fluxweights onto [0,1] {{{
      fluxweight<-magmap(fluxweight,lo=0,hi=max(fluxweight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
      #}}}
    }
    #}}}
    #}}}
  } else if (weightType=="flux") {
    #Fluxweights are Previously determined object Fluxes in Jy {{{
    #Get Range {{{
    wgtRange<-range(fluxweight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(fluxweight))) {
      #Bad values. Things with bad fluxweights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(fluxweight))
      fluxweight[ind]<-NA
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(fluxweight<=0)
      fluxweight[ind]<-NA
      #}}}
    }#}}}
    #Adjust Fluxweights by aperture integral
    fluxweight<-foreach(samask=sa_mask[subs],fluxwgt=fluxweight, .inorder=TRUE, .options.mpi=mpiopts, .combine='c') %dopar%{ fluxwgt/sum(samask)  }
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(fluxweight,na.rm=TRUE)<=0) {
      message("WARNING: No fluxweights are > 0! Flux weighting has removed all objects from the image.")
      fluxweight<-rep(0,length(subs))
    } else if (length(which(fluxweight>0))<2) {
      message("WARNING: only one fluxweight is > 0! Flux weighting has removed all other objects from the image.")
      fluxweight[which(fluxweight<=0)]<-NA
      fluxweight<-magmap(fluxweight,lo=0,hi=max(fluxweight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
    } else {
      #Map Fluxweights onto [0,1] {{{
      fluxweight<-magmap(fluxweight,lo=0,hi=max(fluxweight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
      #}}}
    }
    #}}}
  } else {
    #Weight type must be scale i.e. [0,1] {{{
    #Get Range {{{
    wgtRange<-range(fluxweight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(fluxweight))) {
      #Bad values. Things with bad fluxweights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(fluxweight))
      fluxweight[ind]<-0
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(fluxweight<=0)
      fluxweight[ind]<-0
      #}}}
    }#}}}
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(fluxweight,na.rm=TRUE)<=0) {
      message("WARNING: No fluxweights are > 0! Flux weighting has removed all objects from the image.")
      fluxweight<-rep(0,length(subs))
    } else if (length(which(fluxweight>0))<2) {
      message("WARNING: only one fluxweight is > 0! Flux weighting has removed all other objects from the image.")
    }
    #}}}
  }
  #}}}

  #Perform calculations {{{
  if ((length(which(fluxweight != 1)) != 0)&(!psffilt)) {
    #No Convolution, need Fluxweighting {{{
    #Details {{{
    #If we are not filtering apertures with PSFs, and the fluxweights
    #are not all unity, multiply each stamp by it's fluxweight }}}
    sfa_mask<-foreach(samask=sa_mask[subs],fluxwgt=fluxweight, .inorder=TRUE, .options.mpi=mpiopts) %dopar%{ samask*fluxwgt  }
    #}}}
  } else if (psffilt) {
    #Convolution with PSF & Fluxweighting {{{
    #Details {{{
    #Else if we are filtering apertures with the PSFs, perform the convolution and
    #multiply by the fluxweight simultaneously }}}
    message('--------------------------Convolution----------------------------------')
    sfa_mask<-foreach(samask=sa_mask[subs], slen=stamplen[subs], fluxwgt=fluxweight, i=1:npos, xc=x_g[subs], yc=y_g[subs], .inorder=TRUE,
    .export=c("psf", "diagnostic"), .options.mpi=mpiopts) %dopar% {
      #Use subset of psf conformable with current aperture stamp {{{
      if (slen<length(psf[,1])){
        #Aperture Stamp is smaller than PSF stamp {{{
        #limits from PSF peak - 1/2 stampwidth {{{
        centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
        delta<-floor(slen/2)*c(-1,+1)
        lims<-rbind(centre[1]+delta,centre[2]+delta)
        #}}}
        #Check for -ve indexes or indexes above PSF width {{{
        if (lims[1,1]<1) {
          lims[1,]<-lims[1,]+(1-lims[1,1])
        }
        if (lims[2,1]<1) {
          lims[2,]<-lims[2,]+(1-lims[2,1])
        }
        if (lims[1,2]>length(psf[,1])) {
          lims[1,]<-lims[1,]+(length(psf[,1])-lims[1,2])
        }
        if (lims[2,2]>length(psf[,1])) {
          lims[2,]<-lims[2,]+(length(psf[,1])-lims[2,2])
        }
        #}}}
        #Check that arrays are conformable {{{
        if (length(psf[lims[1]:lims[3],1])!=slen) {
          stop(paste("PSF and Aperture Arrays are non-conformable. Lengths are:",length(psf[lims[1]:lims[3],1]),slen))
        }
        #}}}
        #}}}
      } else if (slen==length(psf[,1])){
      #Aperture stamp is same size as PSF stamp {{{
        lims<-matrix(data=rep(c(1,length(psf[,1])),2),nrow=2, byrow=TRUE)
        #}}}
      } else {
        #Aperture stamp is *larger* than psf stamp - not allowed {{{
        sink(type="message")
        stop("Aperture Stamp is larger than PSF - convolution cannot be performed")
        #}}}
      }
      #}}}
      #Diagnostic {{{
      if (diagnostic) { message(paste("PSF Subset complete in Aperture",i)) }
      #}}}
      #Convolve & Weight {{{
      if (length(which(samask!=0))>1){
        #Aperture {{{
        #Diagnostic {{{
        if (diagnostic) { message(paste("Doing convolution of PSF with Aperture", i)) }
        #}}}
        #Convolve {{{
        ap<-(convolvepsf(psf[lims[1]:lims[3],lims[2]:lims[4]],samask,normalise=TRUE))
        #}}}
        #Remove any Negatives produced by convolution {{{
        ap[which((-log10(abs(ap)))>=(-log10(abs(min(ap)))))]<-0
        #}}}
        #Check for Errors {{{
        if (length(which(ap <0))>0) {
          #If Negatives produced, Error  {{{
          sink(sinkfile,type='output')
          print((summary(as.numeric(ap))))
          print((summary(as.numeric(psf[lims[1]:lims[3],lims[2]:lims[4]]))))
          sink(type='message')
          sink(type='output')
          stop("Unable to remove negatives produced in convolution")
          #}}}
        }
        #}}}
        #Finalise Aperture {{{
        ap<-(ap/max(ap))*fluxwgt
        #}}}
        #Return {{{
        return=ap
        #}}}
        #}}}
      } else {
        #Point Source {{{
        #Recalculate Limits to make sure PSF is centred correctly {{{
        #Aperture Stamp is smaller than PSF stamp {{{
        #limits from PSF peak - 1/2 stampwidth {{{
        centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
        delta<-floor(slen/2)*c(-1,+1)
        lims<-rbind(centre[1]+delta,centre[2]+delta)
        #}}}
        aplims<-rbind(c(1,slen),c(1,slen))
        #Check for -ve indexes or indexes above PSF width {{{
        if (lims[1,1]<1) {
          aplims[1,1]<-aplims[1,1]+(1-lims[1,1])
          lims[1,1]<-1
        }
        if (lims[2,1]<1) {
          aplims[2,1]<-aplims[2,1]+(1-lims[2,1])
          lims[2,1]<-1
        }
        if (lims[1,2]>length(psf[,1])) {
          aplims[1,2]<-(slen-(lims[1,2]-length(psf[1,])))
          lims[1,2]<-length(psf[1,])
        }
        if (lims[2,2]>length(psf[,1])) {
          aplims[2,2]<-(slen-(lims[2,2]-length(psf[,1])))
          lims[2,2]<-length(psf[,1])
        }
        #}}}
        #Check that arrays are conformable {{{
        if (length(psf[lims[1]:lims[3],1])!=length(samask[aplims[1]:aplims[3],1])) {
          stop(paste("PSF and Point Source Aperture Array are non-conformable. Lengths are:",length(psf[lims[1]:lims[3],1]),length(samask[aplims[1]:aplims[3],1])))
        }
        #}}}
        #}}}
        #}}}
        #Reinterpolate the PSF at point source XcenYcen {{{
        lenx<-length(lims[1]:lims[3])
        leny<-length(lims[2]:lims[4])
        #Make grid for psf at old pixel centres {{{
        psf_obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[lims[1]:lims[3],lims[2]:lims[4]])
        #}}}
        #Make expanded grid of new pixel centres {{{
        expanded<-expand.grid(seq(1,lenx),seq(1,leny))
        xnew<-expanded[,1]-xc%%1
        ynew<-expanded[,2]-yc%%1
        #}}}
        #Interpolate {{{
        ap<-matrix(interp2D(xnew, ynew, psf_obj)[,3], ncol=leny,nrow=lenx)
        #}}}
        #}}}
        #Make Final Stamp {{{
        apfin<-matrix(0, nrow=slen, ncol=slen)
        apfin[aplims[1]:aplims[3],aplims[2]:aplims[4]]<-ap
        #}}}
        #Normalise PSF to 1 {{{
        apfin<-(apfin/max(apfin))*fluxwgt
        #}}}
        #Return {{{
        return=apfin
        #}}}
        #}}}
      }
      #}}}
    }
    message('===========END============Convolution==============END=================\n')
    #}}}
  } else {
    #Catch for Errors {{{
    #Details {{{
    #If we are not filtering by psf, and all fluxweights are equal to unity,
    #then we would not have entered this function in the first place. If we
    #arrive here, something has gone wrong with our fluxweights or parsed variables. }}}
    sink(type="message")
    stop("Bad variables and/or fluxweights in production of filtered apertures.")
    #}}}
  }
  #}}}

  #Check that apertures do not cross image mask boundary {{{
  if ((cutup & length(immask)>1) | ((!cutup) & length(immask)!=0 & length(image.env$imm) > 1)) {
    #Check Mask stamps for Aperture Acceptance {{{
    message('Combining Aps with Mask Stamps')
    if (cutup) {
      sff_mask<-foreach(slen=stamplen[subs], smask=sfa_mask,mmask=immask[subs],mxl=mstamp_lims[subs,1],mxh=mstamp_lims[subs,2],myl=mstamp_lims[subs,3],myh=mstamp_lims[subs,4], .export="useMaskLim", .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(mmask[mxl:mxh,myl:myh]*smask,na.rm=TRUE)/sum(smask)
        if (is.na(check)) {
          #All pixels in mask are Na/NaN
          array(0, dim=c(slen,slen))
        } else if (check<useMaskLim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          smask*mmask[mxl:mxh,myl:myh]
          #}}}
        }
        #}}}
      }
    } else {
      sff_mask<-foreach(slen=stamplen[subs], smask=sfa_mask,mxl=mask_lims[subs,1],mxh=mask_lims[subs,2],myl=mask_lims[subs,3],myh=mask_lims[subs,4], .export=c("useMaskLim","image.env"), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(image.env$imm[mxl:mxh,myl:myh]*smask,na.rm=TRUE)/sum(smask)
        if (is.na(check)) {
          #All pixels in mask are Na/NaN
          array(0, dim=c(slen,slen))
        } else if (check<useMaskLim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          smask*image.env$imm[mxl:mxh,myl:myh]
          #}}}
        }
        #}}}
      }
    }
    message('Mask Combine Finished.')
    #}}}
  } else {
    #No image mask, all can be kept {{{
    sff_mask<-sfa_mask
    #}}}
  }#}}}

  #Check for production of NA/NaN/Infs in convolution {{{
  if (diagnostic) { for (i in 1:length(sff_mask)) { if (length(which(is.na(as.numeric(sff_mask[[i]]))))>0) {sink(type="message") ; stop("NAs produced in Convolution")} } }
  if (length(fluxweightin) > 0 ) { message('===========END===========Make_WSFA_Mask============END=================\n')
  } else { message('===========END===========Make_SFA_MASK=============END=================\n') }
  #}}}

  #Return array of Stamps {{{
  if (!is.null(env)) { detatch(env) }
  return=sff_mask
  #}}}
}
