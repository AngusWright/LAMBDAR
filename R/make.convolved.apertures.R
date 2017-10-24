make.convolved.apertures <-
function(outenv=parent.env(environment()), sa_mask,flux.weightin=NULL, immask=NULL, env=NULL,sumsa,subs=NULL) {
#Make the Single Filtered Aperture mask

  #Print appropriate banner  {{{
  if (length(flux.weightin) > 0 ) {
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
    npos<-length(cat.id)
    subs<-1:npos
  } else {
    npos<-length(subs)
    message(paste("-> There are",npos,"apertures in this subset\n"))
  }
  #}}}

  #Check for sumsa (and its size) {{{
  if (!missing(sumsa)) { 
    if (length(sumsa) != npos) {
      if (length(sumsa) == length(cat.id)) {
        sumsa<-sumsa[subs]
      } else {
        sink(type="message") ; stop("Bad SumSA Array")
      }
    }
  }
  #}}}

  #Set Fluxweights {{{
  #If no flux.weight - all weights == 1 {{{
  if (length(flux.weightin) == 0) { flux.weight<-array(1, dim=c(npos))
  #}}}
  #If single value suplied, all weights == value {{{
  } else if (length(flux.weightin) == 1) { flux.weight<-array(flux.weightin, dim=c(npos))
  #}}}
  #If vector of values supplied, weights == input {{{
  } else { flux.weight<-flux.weightin }
  #}}}

  #Check that number of Fluxweights == number of subs stamps or number of total stamps {{{
  if (length(flux.weight) != npos) {
    if (length(flux.weight) == length(cat.id)) {
      flux.weight<-flux.weight[subs]
    } else {
      sink(type="message") ; stop("Bad FluxWeightIn Array")
    }
  }
  #}}}

  #If we have input flux.weights, Check their type {{{
  if (length(flux.weightin)==0) {
    message("No input Fluxweights; no rescale is required")
  } else if (weight.type=="mag") {
    #Fluxweights are Previously determined object magnitudes {{{
    #Convert back to Flux {{{
    flux.weight<-10^((8.9-flux.weight)/2.5)
    #}}}
    #Get Range {{{
    wgtRange<-range(flux.weight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(flux.weight))) {
      #Bad values. Things with bad flux.weights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(flux.weight))
      flux.weight[ind]<-NA
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(flux.weight<=0)
      flux.weight[ind]<-NA
      #}}}
    }#}}}
    #Adjust Fluxweights by aperture integral
    if (missing(sumsa)) {
      flux.weight<-foreach(samask=sa_mask[subs],fluxwgt=flux.weight, .inorder=TRUE, .options.mpi=mpi.opts, .combine='c') %dopar%{ fluxwgt/sum(samask)  }
    } else { 
      flux.weight<-flux.weight/sumsa
    }
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(flux.weight,na.rm=TRUE)<=0) {
      message("WARNING: No flux.weights are > 0! Flux weighting has removed all objects from the image.")
      flux.weight<-rep(0,length(subs))
    } else if (length(which(flux.weight>0))<2) {
      message("WARNING: only one flux.weight is > 0! Flux weighting has removed all other objects from the image.")
      flux.weight[which(flux.weight<=0)]<-NA
      flux.weight<-magmap(flux.weight,lo=0,hi=max(flux.weight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
    } else {
      #Map Fluxweights onto [0,1] {{{
      flux.weight<-magmap(flux.weight,lo=0,hi=max(flux.weight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
      #}}}
    }
    #}}}
  } else if (weight.type=="flux") {
    #Fluxweights are Previously determined object Fluxes in Jy {{{
    #Get Range {{{
    wgtRange<-range(flux.weight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(flux.weight))) {
      #Bad values. Things with bad flux.weights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(flux.weight))
      flux.weight[ind]<-NA
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(flux.weight<=0)
      flux.weight[ind]<-NA
      #}}}
    }#}}}
    #Adjust Fluxweights by aperture integral
    if (missing(sumsa)) {
      flux.weight<-foreach(samask=sa_mask[subs],fluxwgt=flux.weight, .inorder=TRUE, .options.mpi=mpi.opts, .combine='c') %dopar%{ fluxwgt/sum(samask)  }
    } else { 
      flux.weight<-flux.weight/sumsa
    }
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(flux.weight,na.rm=TRUE)<=0) {
      message("WARNING: No flux.weights are > 0! Flux weighting has removed all objects from the image.")
      flux.weight<-rep(0,length(subs))
    } else if (length(which(flux.weight>0))<2) {
      message("WARNING: only one flux.weight is > 0! Flux weighting has removed all other objects from the image.")
      flux.weight[which(flux.weight<=0)]<-NA
      flux.weight<-magmap(flux.weight,lo=0,hi=max(flux.weight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
    } else {
      #Map Fluxweights onto [0,1] {{{
      flux.weight<-magmap(flux.weight,lo=0,hi=max(flux.weight,na.rm=TRUE),bad=0,range=c(0,1),stretch='lin',stretchscale=1, type="num")$map
      #}}}
    }
    #}}}
  } else {
    #Weight type must be scale i.e. [0,1] {{{
    #Get Range {{{
    wgtRange<-range(flux.weight, na.rm=TRUE)
    #}}}
    #Check Bad Values {{{
    if (any(!is.finite(flux.weight))) {
      #Bad values. Things with bad flux.weights will be removed in calcs{{{
      warning("Supplied Fluxweight has values that are not finite (NA/NaN/Inf). These objects will not be fit")
      ind<-which(!is.finite(flux.weight))
      flux.weight[ind]<-0
      #}}}
    }#}}}
    #Check Minima {{{
    if (min(wgtRange)<=0) {
      #Bad minima. Things less than 0 will be removed in calcs{{{
      warning("Supplied Fluxweight has values less than or equal to 0. These objects will not be fit")
      ind<-which(flux.weight<=0)
      flux.weight[ind]<-0
      #}}}
    }#}}}
    #}}}
    #Check for Errors & Fluxweight {{{
    if (max(flux.weight,na.rm=TRUE)<=0) {
      message("WARNING: No flux.weights are > 0! Flux weighting has removed all objects from the image.")
      flux.weight<-rep(0,length(subs))
    } else if (length(which(flux.weight>0))<2) {
      message("WARNING: only one flux.weight is > 0! Flux weighting has removed all other objects from the image.")
    }
    #}}}
  }
  #}}}

  #If grouping weights, recalculate the weights {{{
  if (group.weights) {
    vals<-foreach(sel=levels(factor(groups)),.export=c('flux.weight','groups'),.combine='c') %dopar% { median(flux.weight[which(groups==sel)],na.rm=T) }
    vals<-magmap(vals,range=c(0.1,1),stretch='lin')$map
    names(vals)<-levels(factor(groups))
    for (sel in levels(factor(groups))) {
      flux.weight[which(groups==sel)]<-flux.weight[which(groups==sel)]*vals[sel]
    }
  }
  #}}}
  #}}}

  #Perform calculations {{{
  if ((length(which(flux.weight != 1)) != 0)&(!psf.filt)) {
    #No Convolution, need Fluxweighting {{{
    #Details {{{
    #If we are not filtering apertures with PSFs, and the flux.weights
    #are not all unity, multiply each stamp by it's flux.weight }}}
    sfa_mask<-foreach(samask=sa_mask[subs],fluxwgt=flux.weight, .inorder=TRUE, .options.mpi=mpi.opts) %dopar%{ samask*fluxwgt  }
    #}}}
  } else if (psf.filt) {
    #Convolution with PSF & Fluxweighting {{{
    #Details {{{
    #Else if we are filtering apertures with the PSFs, perform the convolution and
    #multiply by the flux.weight simultaneously }}}
    message('--------------------------Convolution----------------------------------')
    #Check that the PSF is a list of images, and for its length {{{
    if (!is.list(psf)) { 
      psf<-list(psf)
      psf.id<-rep(1,length(cat.x))
    } else if (!exists('psf.id') & length(psf)>1) {
      stop("No PSF Id provided with multiple PSFs") 
    } else if (!exists('psf.id')) {
      psf.id<-rep(1,length(cat.x))
    }
    n.psf<-length(psf)
    #}}}
    #In the event that we have point sources (they are not convolved), check that the PSF is correctly centred {{{
    psf.cen<-psf
    recent<-cbind(rep(NA,n.psf),rep(NA,n.psf))
    for (i in 1:n.psf) { 
      conv<-array(0,dim=dim(psf[[i]]))
      conv[which(psf[[i]]==max(psf[[i]]),arr.ind=T)]<-1
      psf.cen[[i]]<-convolve.psf(psf[[i]],conv)
      psf.resid<-psf[[i]]-psf.cen[[i]]
      if (plot.sample) { 
        PlotPNG(file.path(path.root,path.work,path.out,paste0("PointSourceCentre_test_bin",i,".png")),width=110*8,height=110*4,res=110)
        layout(cbind(1,2,3,4))
        magimage(psf[[i]])
        magimage(psf.cen[[i]])
        magimage(psf.resid)
      }
      if (max(abs(psf.resid),na.rm=T) > 0.001) { 
        recent[i,]<-which(psf.resid==max(psf.resid,na.rm=T),arr.ind=T)[1,]-which(psf[[i]]==max(psf[[i]],na.rm=T),arr.ind=T)[1,]
      } else { 
        recent[i,]<-c(0,0)
      }
      if (any(abs(recent[i,])>1)) { recent[i,which(abs(recent[i,])>1)]<-0 }
      recent<-recent/2
      #Make grid for psf at old pixel centres {{{
      psf.obj<-list(x=seq(1,dim(psf[[i]])[1])+recent[i,1], y=seq(1,dim(psf[[i]])[2])+recent[i,2],z=psf[[i]])
      #}}}
      #Make expanded grid of new pixel centres {{{
      expanded<-expand.grid(seq(1,dim(psf[[i]])[1]),seq(1,dim(psf[[i]])[2]))
      xnew<-expanded[,1]
      ynew<-expanded[,2]
      #}}}
      #Interpolate {{{
      psf.new<-matrix(interp.2d(xnew, ynew, psf.obj)[,3], ncol=dim(psf[[i]])[2],nrow=dim(psf[[i]])[1])
      #}}}
      #Check the new psf residuals {{{
      psf.resid<-psf.new-psf.cen[[i]]
      new.recent<-which(psf.resid==max(psf.resid,na.rm=T),arr.ind=T)[1,]-which(psf[[i]]==max(psf[[i]],na.rm=T),arr.ind=T)[1,]
      if (any(new.recent!=0)) { 
        #The recentering didnt work, just use it as is
        recent[i,]<-c(0,0)
      }
      #}}}
      #}}}
      if(plot.sample) { 
        magimage(psf.resid)
        label('top',lab=paste0('rerecentered: recentre fact = c(',new.recent[1],',',new.recent[2],')'),col='red')
        label('bottom',lab=paste0('diagnosis: recentre fact = c(',recent[i,1],',',recent[i,2],')'),col='red')
        dev.off()
      }
    }
    #}}}
    sfa_mask<-foreach(samask=sa_mask[subs], slen=stamplen[subs], pid=psf.id[subs], fluxwgt=flux.weight, i=1:npos, xc=cat.x[subs], yc=cat.y[subs], .inorder=TRUE,
    .export=c("psf","diagnostic","recent"), .options.mpi=mpi.opts) %dopar% {
      #Use subset of psf conformable with current aperture stamp {{{
      if (slen<length(psf[[pid]][,1])){
        #Aperture Stamp is smaller than PSF stamp {{{
        #limits from PSF peak - 1/2 stampwidth {{{
        centre<-as.numeric(which(psf[[pid]]==max(psf[[pid]]), arr.ind=TRUE))
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
        if (lims[1,2]>length(psf[[pid]][,1])) {
          lims[1,]<-lims[1,]+(length(psf[[pid]][,1])-lims[1,2])
        }
        if (lims[2,2]>length(psf[[pid]][,1])) {
          lims[2,]<-lims[2,]+(length(psf[[pid]][,1])-lims[2,2])
        }
        #}}}
        #Check that arrays are conformable {{{
        if (length(psf[[pid]][lims[1]:lims[3],1])!=slen) {
          stop(paste("PSF and Aperture Arrays are non-conformable. Lengths are:",length(psf[[pid]][lims[1]:lims[3],1]),slen))
        }
        #}}}
        #}}}
      } else if (slen==length(psf[[pid]][,1])){
      #Aperture stamp is same size as PSF stamp {{{
        lims<-matrix(data=rep(c(1,length(psf[[pid]][,1])),2),nrow=2, byrow=TRUE)
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
        ap<-(convolve.psf(psf[[pid]][lims[1]:lims[3],lims[2]:lims[4]],samask,normalise=TRUE))
        #}}}
        #Remove any Negatives produced by convolution {{{
        ap[which((-log10(abs(ap)))>=(-log10(abs(min(ap)))))]<-0
        #}}}
        #Check for Errors {{{
        if (length(which(ap <0))>0) {
          #If Negatives produced, Error  {{{
          sink(sink.file,type='output')
          print((summary(as.numeric(ap))))
          print((summary(as.numeric(psf[[pid]][lims[1]:lims[3],lims[2]:lims[4]]))))
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
        centre<-as.numeric(which(psf[[pid]]==max(psf[[pid]]), arr.ind=TRUE))+recent
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
        if (lims[1,2]>length(psf[[pid]][,1])) {
          aplims[1,2]<-(slen-(lims[1,2]-length(psf[[pid]][1,])))
          lims[1,2]<-length(psf[[pid]][1,])
        }
        if (lims[2,2]>length(psf[[pid]][,1])) {
          aplims[2,2]<-(slen-(lims[2,2]-length(psf[[pid]][,1])))
          lims[2,2]<-length(psf[[pid]][,1])
        }
        #}}}
        #Check that arrays are conformable {{{
        if (length(psf[[pid]][lims[1]:lims[3],1])!=length(samask[aplims[1]:aplims[3],1])) {
          stop(paste("PSF and Point Source Aperture Array are non-conformable. Lengths are:",length(psf[[pid]][lims[1]:lims[3],1]),length(samask[aplims[1]:aplims[3],1])))
        }
        #}}}
        #}}}
        #}}}
        #Reinterpolate the PSF at point source XcenYcen {{{
        lenx<-length(lims[1]:lims[3])
        leny<-length(lims[2]:lims[4])
        #Make grid for psf at old pixel centres {{{
        psf.obj<-list(x=seq(1,lenx), y=seq(1,leny),z=psf[[pid]][lims[1]:lims[3],lims[2]:lims[4]])
        #}}}
        #Make expanded grid of new pixel centres {{{
        expanded<-expand.grid(seq(1,lenx),seq(1,leny))
        xnew<-expanded[,1]-xc%%1
        ynew<-expanded[,2]-yc%%1
        #}}}
        #Interpolate {{{
        ap<-matrix(interp.2d(xnew, ynew, psf.obj)[,3], ncol=leny,nrow=lenx)
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
    #If we are not filtering by psf, and all flux.weights are equal to unity,
    #then we would not have entered this function in the first place. If we
    #arrive here, something has gone wrong with our flux.weights or parsed variables. }}}
    sink(type="message")
    stop("Bad variables and/or flux.weights in production of filtered apertures.")
    #}}}
  }
  #}}}

  #Check that apertures do not cross image mask boundary {{{
  if ((cutup & length(immask)>1) | ((!cutup) & length(immask)!=0 & length(image.env$imm) > 1)) {
    #Check Mask stamps for Aperture Acceptance {{{
    message('Combining Aps with Mask Stamps')
    if (cutup) {
      sff_mask<-foreach(slen=stamplen[subs], smask=sfa_mask,mmask=immask[subs],mxl=ap.lims.mask.stamp[subs,1],mxh=ap.lims.mask.stamp[subs,2],myl=ap.lims.mask.stamp[subs,3],myh=ap.lims.mask.stamp[subs,4], .export="use.mask.lim", .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(mmask[mxl:mxh,myl:myh]*smask,na.rm=TRUE)/sum(smask)
        if (is.na(check)) {
          #All pixels in mask are Na/NaN
          array(0, dim=c(slen,slen))
        } else if (check<use.mask.lim) {
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
      sff_mask<-foreach(slen=stamplen[subs], smask=sfa_mask,mxl=ap.lims.mask.map[subs,1],mxh=ap.lims.mask.map[subs,2],myl=ap.lims.mask.map[subs,3],myh=ap.lims.mask.map[subs,4], .export=c("use.mask.lim","image.env"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(image.env$imm[mxl:mxh,myl:myh]*smask,na.rm=TRUE)/sum(smask)
        if (is.na(check)) {
          #All pixels in mask are Na/NaN
          array(0, dim=c(slen,slen))
        } else if (check<use.mask.lim) {
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
  if (length(flux.weightin) > 0 ) { message('===========END===========Make_WSFA_Mask============END=================\n')
  } else { message('===========END===========Make_SFA_MASK=============END=================\n') }
  #}}}

  #Return array of Stamps {{{
  if (!is.null(env)) { detatch(env) }
  return=sff_mask
  #}}}
}
