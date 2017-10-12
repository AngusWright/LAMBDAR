estimate.psf <-
function (outenv=parent.env(environment()),n.bins=1,n.sources=3e3,bin.type='SNR.quan',lo=20,hi=200,type='num',blend.tolerance=0.5,env=NULL,plot=FALSE) {

  message('--------------------------Estimate_PSF-------------------------------------')
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

  # Identify the point sources
  point.sources<-which(cat.a==0) 
  # Reinterpolate and stack the point sources
  mat<-matrix(0,max(stamplen[point.sources]),max(stamplen[point.sources]))
  im_psf.nomask<-im_psf<-weight<-nomask.n<-list()
  for (i in 1:n.bins) { 
    im_psf.nomask[[i]]<-im_psf[[i]]<-weight[[i]]<-mat
    nomask.n[[i]]<-0
  }
  if (do.sky.est & exists('skylocal')) { 
    pixval<-image.env$im[cbind(cat.x,cat.y)]-skylocal
  } else if (do.sky.est) { 
    if (cutup) {
      timer<-system.time(skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                      cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,data.stamp=data.stamp,mask.stamp=mask.stamp,
                      clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources))
    } else {
      timer<-system.time(skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                      cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                      data.stamp=image.env$im, mask.stamp=image.env$imm.dimim,
                      clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources))
    }
    skylocal<-rep(NA,length(cat.x)) 
    skylocal[point.sources]<-skyest[,'sky']
    skyrms<-rep(NA,length(cat.x)) 
    skyrms[point.sources]<-skyest[,'skyRMS']
    pixval<-image.env$im[cbind(cat.x,cat.y)]
    pixval<-pixval[point.sources]-skylocal
  } else {  
    pixval<-image.env$im[cbind(cat.x,cat.y)]
  }
  
  if (grepl('SNR',bin.type)) {
    if (!(do.sky.est|get.sky.rms)) { stop("Error: cannot SNR bin without RMS estimate!") }
    pixval<-pixval/skyrms
  }

  if (grepl("quan",bin.type)) { 
    if (type=='quan') { 
      #quantile bin limits 
      bin.lim<-quantile(pixval[is.finite(pixval)],seq(lo,hi,length=n.bins+1))
    } else { 
      #Absolute bins limits
      bin.lim<-quantile(pixval[which(pixval>=lo & pixval<=hi)],seq(0,1,length=n.bins+1))
    }
    bin.zero<-bin.lim[1]
    bin.lim<-bin.lim[-1]
  } else { 
    if (type=='quan') { 
      #quantile bin limits 
      quans<-quantile(pixval[is.finite(pixval)],c(lo,hi))
      bin.lim<-seq(quans[1],quans[2],length=n.bins+1)
    } else { 
      #Absolute bins limits
      bin.lim<-seq(lo,hi,length=n.bins+1)
    }
    #Equal spaced bins 
    bin.zero<-bin.lim[1]
    bin.lim<-bin.lim[-1]
  } 
  #If the pixval is outside the bins, skip it 
  keep <- which(pixval[point.sources] >= bin.zero & pixval[point.sources] <= max(bin.lim)) 
  point.sources<-point.sources[keep]
  #Assign the bins
  bin<-rep(-1,length(pixval)) 
  for (i in 1:n.bins) { 
    bin[which(bin<0 & pixval<=bin.lim[i] & pixval >= bin.zero)]<-i
  }
  
  im_psf_cen<-dim(im_psf)/2
  #Remove sources which have masked _datamaps_ (can create artefacts under convolution)
  if (length(image.env$imm) > 1) { 
    maskfrac<-rep(NA,length=length(cat.x))
    for (i in point.sources) {
      maskfrac[i]<-sum(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]]==0)
      maskfrac[i]<-maskfrac[i]/length(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]])
    }
  }
  if (exists('maskfrac')) { 
    keep<-which(maskfrac[point.sources]<=0.2)
    point.sources<-point.sources[keep]
    maskfrac<-maskfrac[point.sources[keep]]
    point.sources<-point.sources[order(maskfrac)]
  } 
  #Calculate the blend fractions
  if (exists('dfa')) { 
    #Use the measured DFAs 
    blendfrac<-rep(NA,length(pixval))
    j=1
    for (i in point.sources) { 
      blendfrac[i]<-sum(dfa[[i]])/sum(sfa[[i]])
    }
    #Remove sources that are too blended
    keep<-which(blendfrac[point.sources]>=blend.tolerance)
    point.sources<-point.sources[keep]
    #Order: least to most blended   
    point.sources<-point.sources[order(blendfrac[point.sources],decreasing=TRUE)]
  } else { 
    #Use pixel-space nearest neighbours 
    match<-nn2(data.frame(cat.x,cat.y),data.frame(cat.x,cat.y)[point.sources,],searchtype='radius',radius=20,k=10)
    #Order by the nearest non-self match (2nd nnd column)
    point.sources<-point.sources[order(match$nn.dists[,2],decreasing=TRUE)]
    nn.dist<-match$nn.dists[order(match$nn.dists[,2],decreasing=TRUE),2]
    #Reject sources that are, assuming at-least Nyquist sampling, within 3sigma overlap of the point source
    point.sources<-point.sources[-which(nn.dist<9)]
  }
  #Remove sources above and beyond what is requested
  for (i in 1:n.bins) { 
    keep<-which(bin[point.sources]==i)
    if (length(keep) > n.sources) { 
      keep<-keep[(n.sources+1):length(keep)]
      point.sources<-point.sources[-keep]
    }
  }
  if (plot) { 
    #show a sample of the PSFs used in the stack
    nsamp<-min(18,length(point.sources))
    sample=sample(point.sources,nsamp)
    layout(matrix(c(1:(nsamp*2),rep(0,length=ceiling(nsamp/2)*2-nsamp)),ncol=6,byrow=T))
    par(mar=c(0,0,0,0),oma=c(2,2,2,2))
  }
  #Loop through the sources remaining
  for (i in point.sources) { 
    slen=stamplen[i]
    xc=cat.x[i]
    yc=cat.y[i]
    if (cutup) {
      im=data.stamp[[i]]
      mask=mask.stamp[[i]]
    } else {
      im=image.env$im
      mask=mask.stamp
      if (length(mask.stamp)==1) { 
        mask=image.env$imm
      }
    }
    xlo=ap.lims.data.stamp[i,1]
    xup=ap.lims.data.stamp[i,2]
    ylo=ap.lims.data.stamp[i,3]
    yup=ap.lims.data.stamp[i,4]
    xmlo=ap.lims.mask.stamp[i,1]
    xmup=ap.lims.mask.stamp[i,2]
    ymlo=ap.lims.mask.stamp[i,3]
    ymup=ap.lims.mask.stamp[i,4]
    if (do.sky.est) { 
      sky=skylocal[i]
      rms=skyrms[i]
    } else { 
      sky=0
      rms=0
    }
    
    #Skip the source if it's saturated 
    if (any(im[xlo:xup,ylo:yup]>=image.env$saturation)) { next } 
    #Remove any sub-pixel centroiding 
    #Make grid for psf at old pixel centres /*fold*/ {{{
    im.obj<-list(x=xlo:xup, y=ylo:yup,z=im[xlo:xup,ylo:yup]-sky)
    # /*fend*/ }}}
    #Make expanded grid of new pixel centres /*fold*/ {{{
    expanded<-expand.grid(xlo:xup,ylo:yup)
    xnew<-expanded[,1]-xc%%1
    ynew<-expanded[,2]-yc%%1
    # /*fend*/ }}}
    #Interpolate /*fold*/ {{{
    im.cen<-matrix(interp.2d(xnew, ynew, im.obj)[,3], ncol=ncol(im.obj$z),nrow=nrow(im.obj$z))
    # /*fend*/ }}}
    #conv<-array(0,dim=dim(im[xlo:xup,ylo:yup]))
    #conv[dim(conv)[1]/2,dim(conv)[2]/2]<-1
    #im.cen<-convolve.psf(im[xlo:xup,ylo:yup]-sky,conv)
    #}
    #Ensure that the central source is unmasked (except when blended)
    #And add to the stack
    if (length(mask)==1) { 
      if (mask==1) { 
        #Add the image to the stack
        im_psf[[bin[i]]]<-(im_psf[[bin[i]]]+(im.cen))
        weight[[bin[i]]]<-(weight[[bin[i]]])+1
      } 
    } else { 
      if (exists('dfa')) { 
        mask[xmlo:xmup,ymlo:ymup][which(zapsmall(dfa[[i]])>= (1-sourcemask.conf.lim))]<-1
      }
      #Add the image to the stack
      im_psf[[bin[i]]]<-(im_psf[[bin[i]]]+(im.cen)*(mask[xmlo:xmup,ymlo:ymup]))
      weight[[bin[i]]]<-(weight[[bin[i]]]+(mask[xmlo:xmup,ymlo:ymup]))
    }
    im_psf.nomask[[bin[i]]]<-(im_psf.nomask[[bin[i]]]+(im.cen))
    nomask.n[[bin[i]]]<-nomask.n[[bin[i]]]+1
    if (plot & i%in%sample) { 
      if (length(mask)==1) { 
        capture=magimage(im[xlo:xup,ylo:yup],axes=FALSE)
        label('topleft',paste0(i,': (1) raw (there is no mask)'),col='red',lwd=2)
        magimage(im.cen,axes=FALSE)
        label('topleft',paste0(i,': (2) recentred (no mask)'),col='red',lwd=2)
        capture=magimage(im.cen-im[xlo:xup,ylo:yup],axes=FALSE)
        label('topleft',paste0(i,': Raw-centred (no mask)'),col='red',lwd=2)
      } else { 
        capture=magimage(im[xlo:xup,ylo:yup]*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': (1) raw (masked)'),col='red',lwd=2)
        magimage(im.cen*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': (2) recentred (masked)'),col='red',lwd=2)
        capture=magimage((im.cen-im[xlo:xup,ylo:yup])*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': Raw-centred (masked)'),col='red',lwd=2)
      }
    }

  }
  # Divide by the per-pixel weights
  for (i in 1:n.bins) { 
    im_psf[[i]]<-im_psf[[i]]/weight[[i]]
    im_psf.nomask[[i]]<-im_psf[[i]]/nomask.n[[i]]
  }


  if (plot & length(point.sources) > 0) { 
    mtext(side=3,outer=T,text=paste0("There are ",length(point.sources),"sources for ",n.bins," bins"))
    layout(matrix(1:(n.bins*3),ncol=n.bins))
    par(mar=c(0,0,0,0),oma=c(2,2,2,2))
    for (j in 1:n.bins) { 
      catch=capture.output(magimage(im_psf[[j]],axes=FALSE))
      label('topleft',paste0('Final PSF estimate'),col='red',lwd=2)
      catch=capture.output(magimage(weight[[j]],axes=FALSE))
      label('topleft',paste0('Masking (Weights)'),col='red',lwd=2)
      catch=capture.output(magimage(im_psf.nomask[[j]],axes=FALSE))
      label('topleft',paste0('Without Mask Estimate'),col='red',lwd=2)
    }
  }

  #Parse Parameter Space & Return{{{
  if (!is.null(env)) { detach(env) }
  #assign("psf.clip" , psf.clip  , envir = outenv)
  #assign("psfwidth" , psfwidth  , envir = outenv)
  #assign("sumpsf"   , sumpsf    , envir = outenv)
  return=list(PSF=im_psf,WEIGHT=weight,NOMASK=im_psf.nomask)
  #}}}
}
