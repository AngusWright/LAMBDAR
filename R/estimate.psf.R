estimate.psf <-
function (outenv=parent.env(environment()),n.bins=1,bloom.bin=TRUE,n.sources=1e3,bin.type='SNR.quan',lo=20,hi=200,type='num',blend.tolerance=0.5,env=NULL,plot=FALSE) {

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

  # Identify the point sources we want to try and stack 
  if (exists('sdfa') & exists('ssfa')) { 
    blendfrac<-sdfa/ssfa
    point.sources<-which(cat.a==min.ap.rad & blendfrac > blend.tolerance) 
  } else { 
    point.sources<-which(cat.a==min.ap.rad) 
  } 
  if (do.sky.est & exists('skylocal')) { 
    pixval.all<-pixval<-image.env$im[cbind(cat.x,cat.y)]-skylocal
  } else if (do.sky.est) { 
    if (cutup) {
      if (quick.sky) { 
        message("Perfoming Fast Sky Estimation")
        skyest<-fast.sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,fit.gauss=fit.sky,
                    cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,data.stamp=data.stamp,mask.stamp=mask.stamp,
                    clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources)
        skyest$sources<-cat.id[point.sources]
        if (fit.sky) { 
          skylocal<-rep(NA,length(cat.x)) 
          skylocal[point.sources]<-skyest[,'skyMu']
          skyrms<-rep(NA,length(cat.x)) 
          skyrms[point.sources]<-skyest[,'skySD']
        } else { 
          skylocal<-rep(NA,length(cat.x)) 
          skylocal[point.sources]<-skyest[,'skyMedian']
          skyrms<-rep(NA,length(cat.x)) 
          skyrms[point.sources]<-skyest[,'skyRMS']
        }
      } else { 
        message("Perfoming Sky Estimation")
        skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                    cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,data.stamp=data.stamp,mask.stamp=mask.stamp,
                    clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources)
        skyest$sources<-cat.id[point.sources]
        skylocal<-rep(NA,length(cat.x)) 
        skylocal[point.sources]<-skyest[,'sky']
        skyrms<-rep(NA,length(cat.x)) 
        skyrms[point.sources]<-skyest[,'skyRMS']
      }
    } else {
      if (quick.sky) { 
        message("Perfoming Fast Sky Estimation") 
        skyest<-fast.sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,fit.gauss=fit.sky,
                    cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                    data.stamp=image.env$im, mask.stamp=image.env$imm.dimim,
                    clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources)
        skyest$sources<-cat.id[point.sources]
        if (fit.sky) { 
          skylocal<-rep(NA,length(cat.x)) 
          skylocal[point.sources]<-skyest[,'skyMu']
          skyrms<-rep(NA,length(cat.x)) 
          skyrms[point.sources]<-skyest[,'skySD']
        } else { 
          skylocal<-rep(NA,length(cat.x)) 
          skylocal[point.sources]<-skyest[,'skyMedian']
          skyrms<-rep(NA,length(cat.x)) 
          skyrms[point.sources]<-skyest[,'skyRMS']
        }
      } else { 
        message("Perfoming Sky Estimation") 
        skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,
                    cutlo=(cat.a/arcsec.per.pix),cuthi=(cat.a/arcsec.per.pix)*5,
                    data.stamp=image.env$im, mask.stamp=image.env$imm.dimim,
                    clipiters=sky.clip.iters,sigma.cut=sky.clip.prob,PSFFWHMinPIX=psffwhm, mpi.opts=mpi.opts,subset=point.sources)
        skyest$sources<-cat.id[point.sources]
        skylocal<-rep(NA,length(cat.x)) 
        skylocal[point.sources]<-skyest[,'sky']
        skyrms<-rep(NA,length(cat.x)) 
        skyrms[point.sources]<-skyest[,'skyRMS']
      }
    }
    pixval.all<-image.env$im[cbind(cat.x,cat.y)]
    pixval<-pixval.all-skylocal
    pixval.all<-pixval.all-median(skylocal,na.rm=TRUE)
  } else {  
    pixval.all<-pixval<-image.env$im[cbind(cat.x,cat.y)]
  }
  
  if (grepl('SNR',bin.type)) {
    if (!(do.sky.est|get.sky.rms)) { 
      message("WARNING: cannot SNR bin without RMS estimate!") 
    } else { 
      pixval<-pixval/skyrms
      pixval.all<-pixval.all/median(skyrms,na.rm=TRUE)
    }
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
  #Do we want an additional blooming bin?
  if (bloom.bin) { 
    new.bin<-image.env$saturation
    if (do.sky.est) { new.bin<-new.bin-median(skylocal,na.rm=TRUE) }
    if (grepl('SNR',bin.type)) { new.bin<-new.bin/median(skyrms,na.rm=TRUE) }
    bin.lim<-c(bin.lim,new.bin) 
    n.bins<-n.bins+1
  }
  #If the pixval is outside the bins, skip it 
  keep <- which(pixval[point.sources] >= bin.zero & pixval[point.sources] <= max(bin.lim)) 
  point.sources<-point.sources[keep]
  #Assign the bins
  bin<-rep(-1,length(pixval)) 
  for (i in 1:n.bins) { 
    bin[point.sources][which(bin[point.sources]<0 & pixval[point.sources]<=bin.lim[i] & pixval[point.sources] >= bin.zero)]<-i
  }
  
  #Remove sources which have masked _datamaps_ (can create artefacts under convolution)
  if (length(image.env$imm.orig) > 1) { 
    maskfrac<-rep(NA,length=length(cat.x))
    for (i in point.sources) {
      maskfrac[i]<-sum(image.env$imm.orig[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]]==0)
      maskfrac[i]<-maskfrac[i]/length(image.env$imm.orig[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]])
      if (exists('dfa')) { 
        maskfrac[i]<-maskfrac[i]-length(which(sfa[[i]]>=1-sourcemask.conf.lim))/length(sfa[[i]])
      }
    }
  } else if (length(image.env$imm) > 1) {
    maskfrac<-rep(NA,length=length(cat.x))
    for (i in point.sources) {
      maskfrac[i]<-sum(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]]==0)
      maskfrac[i]<-maskfrac[i]/length(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]])
      if (exists('dfa')) { 
        maskfrac[i]<-maskfrac[i]-length(which(sfa[[i]]>=1-sourcemask.conf.lim))/length(sfa[[i]])
      }
    }
  }
  if (exists('maskfrac')) { 
    keep<-which(maskfrac[point.sources]<=0.2)
    maskfrac<-maskfrac[point.sources[keep]]
    point.sources<-point.sources[keep]
    point.sources<-point.sources[order(maskfrac)]
  } 
  #Calculate the blend fractions
  if (exists('sdfa') & exists('ssfa')) { 
    if (length(point.sources) > 0) { 
      #Use pixel-space nearest neighbours 
      match<-nn2(data.frame(cat.x,cat.y)[which(blendfrac>blend.tolerance),],data.frame(cat.x,cat.y)[point.sources,],searchtype='radius',radius=20,k=10)
      #Order by the nearest non-self match (2nd nnd column)
      point.sources<-point.sources[order(pixval[point.sources],match$nn.dists[,2],decreasing=TRUE)]
      nn.dist<-match$nn.dists[order(pixval[point.sources],match$nn.dists[,2],decreasing=TRUE),2]
      #Reject sources that are, assuming at-least Nyquist sampling, within 3sigma overlap of the point source
      if (any(nn.dist<9)) { 
        point.sources<-point.sources[-which(nn.dist<9)]
      }
    }
  } else { 
    if (length(point.sources) > 0) { 
      #Use pixel-space nearest neighbours 
      match<-nn2(data.frame(cat.x,cat.y),data.frame(cat.x,cat.y)[point.sources,],searchtype='radius',radius=20,k=10)
      #Order by the nearest non-self match (2nd nnd column)
      point.sources<-point.sources[order(pixval[point.sources],match$nn.dists[,2],decreasing=TRUE)]
      nn.dist<-match$nn.dists[order(pixval[point.sources],match$nn.dists[,2],decreasing=TRUE),2]
      #Reject sources that are, assuming at-least Nyquist sampling, within 3sigma overlap of the point source
      if (any(nn.dist<9)) { 
        point.sources<-point.sources[-which(nn.dist<9)]
      }
    }
  }
  # Initialise the arrays 
  im_psf.nomask<-im_psf<-weight<-nomask.n<-list()
  for (i in 1:n.bins) { 
    if (length(which(bin[point.sources]==i))>0){ 
      mat<-matrix(0,max(stamplen[point.sources][which(bin[point.sources]==i)]),max(stamplen[point.sources][which(bin[point.sources]==i)]))
    } else {
      mat<-matrix(0,min(stamplen),min(stamplen))
    }
    im_psf.nomask[[i]]<-im_psf[[i]]<-weight[[i]]<-mat
    nomask.n[[i]]<-0
  }
  #Remove sources above and beyond what is requested
  for (i in 1:n.bins) { 
    keep<-which(bin[point.sources]==i)
    if (length(keep) > n.sources) { 
      keep<-keep[(n.sources+1):length(keep)]
      point.sources<-point.sources[-keep]
    }
  }
  if (plot & length(point.sources)>0) { 
    #show a sample of the PSFs used in the stack
    nsamp<-min(24,length(point.sources))
    sample=sample(point.sources,nsamp)
    layout(matrix(c(1:(ifelse(nsamp>12,12,nsamp)*3),rep(0,36-ifelse(nsamp>12,12,nsamp)*3)),ncol=6,byrow=T))
    par(mar=c(0,0,0,0),oma=c(2,2,2,2))
  } else { 
    sample<-NULL
  }
  #Loop through the sources remaining
  im=image.env$im
  if (length(image.env$imm.orig) > 1) { 
    mask=image.env$imm.orig
  } else { 
    mask=image.env$imm
  }
  for (i in point.sources) { 
    xc=cat.x[i]
    yc=cat.y[i]
    xlo=ap.lims.data.map[i,1]
    xup=ap.lims.data.map[i,2]
    ylo=ap.lims.data.map[i,3]
    yup=ap.lims.data.map[i,4]
    xmlo=ap.lims.mask.map[i,1]
    xmup=ap.lims.mask.map[i,2]
    ymlo=ap.lims.mask.map[i,3]
    ymup=ap.lims.mask.map[i,4]
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
        mask[xmlo:xmup,ymlo:ymup][which(dfa[[i]]>= (1-sourcemask.conf.lim))]<-1
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
      if (!all(is.na(im_psf[[j]]))) { 
        catch=capture.output(magimage(im_psf[[j]],axes=FALSE))
      } else { 
        image(axes=FALSE,matrix(1),col='white')
        label('centre','NO ESTIMATE')
      }
      label('topleft',paste0('Final PSF estimate'),col='red',lwd=2)
      if (!all(is.na(weight[[j]]))) { 
        catch=capture.output(magimage(weight[[j]],axes=FALSE))
      } else { 
        image(axes=FALSE,matrix(1),col='white')
        label('centre','NO ESTIMATE')
      }
      label('topleft',paste0('Masking (Weights)'),col='red',lwd=2)
      if (!all(is.na(im_psf.nomask[[j]]))) { 
        catch=capture.output(magimage(im_psf.nomask[[j]],axes=FALSE))
      } else { 
        image(axes=FALSE,matrix(1),col='white')
        label('centre','NO ESTIMATE')
      }
      label('topleft',paste0('Without Mask Estimate'),col='red',lwd=2)
    }
  }
  #Put the PSF onto the big stamp {{{ 
  centre<-dim(im_psf[[1]])/2
  #}}}
  #Update bin definition {{{ 
  #Used bin limits 
  bin<-rep(-1,length(pixval.all)) 
  bin[which(pixval.all < bin.zero)]<-1
  bin[which(pixval.all > max(bin.lim))]<-n.bins
  for (i in 1:n.bins) { 
    bin[which(bin<0 & pixval.all < bin.lim[i])]<-i
  }
  #}}}
  #Parse Parameter Space & Return{{{
  if (!is.null(env)) { detach(env) }
  if (exists('skyest')) { 
    assign("tmp.skyest" , skyest  , envir = outenv)
  }
  assign("tmp.psf.id" , bin, envir = outenv)
  assign("tmp.psfest.val",pixval.all, envir = outenv)
  #assign("psfwidth" , psfwidth  , envir = outenv)
  #assign("sumpsf"   , sumpsf    , envir = outenv)
  return=list(PSF=im_psf,WEIGHT=weight,NOMASK=im_psf.nomask,centre=centre)
  #}}}
}
