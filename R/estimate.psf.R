estimate.psf <-
function (outenv=parent.env(environment()),n.bins=1,bloom.bin=FALSE,n.sources=5e2,onlyContams=TRUE,bin.type='SNR.quan',
          lo=20,hi=200,type='num',check.one.sky=length(point.sources)>5*n.sources,blend.tolerance=0.5,
          mask.tolerance=0.0,radial.tolerance=25,all.limit=0.15,env=NULL,plot=FALSE) {

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
    blendfrac<-1-sdfa/ssfa
    if (onlyContams & filt.contam) { 
      point.sources<-which(cat.a==min.ap.rad & blendfrac <= blend.tolerance & contams==1) 
    } else { 
      point.sources<-which(cat.a==min.ap.rad & blendfrac <= blend.tolerance) 
    } 
  } else { 
    blendfrac<-rep(0,length(cat.a))
    if (onlyContams & filt.contam) { 
      point.sources<-which(cat.a==min.ap.rad & contams==1) 
    } else { 
      point.sources<-which(cat.a==min.ap.rad) 
    }
  } 

  #Remove things that are blended 
  if (exists('sdfa') & exists('ssfa')) { 
    if (length(point.sources) > 1) { 
      #Use pixel-space nearest neighbours 
      match<-nn2(data.frame(cat.x,cat.y)[point.sources,][which(blendfrac[point.sources]<=blend.tolerance),],data.frame(cat.x,cat.y)[point.sources,],searchtype='radius',
                 radius=radial.tolerance*2.0,k=min(10,length(which(blendfrac[point.sources]<=blend.tolerance))))
      #Order by the nearest non-self match (2nd nnd column)
      point.sources<-point.sources[order(match$nn.dists[,2],decreasing=TRUE)]
      nn.dist<-match$nn.dists[order(match$nn.dists[,2],decreasing=TRUE),2]
      #Reject sources that are, assuming at-least Nyquist sampling, within 3sigma overlap of the point source
      if (length(which(nn.dist<1e4))/length(nn.dist)<all.limit) { 
        #Just remove all the blends
        point.sources<-point.sources[-which(nn.dist<1e4)]
        nn.dist<-nn.dist[-which(nn.dist<1e4)]
      } else if (any(nn.dist<radial.tolerance)) { 
        point.sources<-point.sources[-which(nn.dist<radial.tolerance)]
        nn.dist<-nn.dist[-which(nn.dist<radial.tolerance)]
      }
    } else { 
      nn.dist<-NULL
    }
  } else { 
    if (length(point.sources) > 0) { 
      #Use pixel-space nearest neighbours 
      match<-nn2(data.frame(cat.x,cat.y),data.frame(cat.x,cat.y)[point.sources,],searchtype='radius',radius=radial.tolerance*2.0,k=min(length(cat.x),10))
      #Order by the nearest non-self match (2nd nnd column)
      nn.dist<-match$nn.dists[order(match$nn.dists[,2],decreasing=TRUE),2]
      point.sources<-point.sources[order(match$nn.dists[,2],decreasing=TRUE)]
      #Reject sources that are, assuming at-least Nyquist sampling, within 3sigma overlap of the point source
      if (length(which(nn.dist<1e4))/length(nn.dist)<all.limit) { 
        #Just remove all the blends
        point.sources<-point.sources[-which(nn.dist<1e4)]
        nn.dist<-nn.dist[-which(nn.dist<1e4)]
      } else if (any(nn.dist<radial.tolerance)) { 
        point.sources<-point.sources[-which(nn.dist<radial.tolerance)]
        nn.dist<-nn.dist[-which(nn.dist<radial.tolerance)]
      }
    } else { 
      nn.dist<-NULL
    }
  }

  if (do.sky.est & exists('skylocal')) { 
    pixval.all<-pixval<-image.env$im[cbind(cat.x,cat.y)]-skylocal
  } else if (do.sky.est) { 
    if (check.one.sky) { 
      #Remove things with pixel values far outside what is requested
      if (length(image.env$imm)>1) { 
        skypix<-image.env$im
        skypix[which(image.env$imm==0)]<-NA
        skypix<-skypix[-(cat.x + nrow(image.env$im) * (cat.y - 1))]
        skypix<-skypix[which(is.finite(skypix))]
      } else { 
        skypix<-image.env$im[-(cat.x + nrow(image.env$im) * (cat.y - 1))]
      }
      skypix<-skypix[which(abs(skypix-median(skypix,na.rm=T))<10*mad(skypix,na.rm=T))]
      onesky<-try(fit.gauss2low(skypix))
      if (class(onesky)=='try-error') { 
        onesky<-data.frame(mu=median(skypix,na.rm=TRUE),sd=mad(skypix,na.rm=T))
      }
      pixval<-image.env$im[cbind(cat.x[point.sources],cat.y[point.sources])] - onesky$mu
      if (grepl("SNR",bin.type)) { 
        pixval<-pixval/onesky$sd
      }
      if (grepl("quan",bin.type)) { 
        if (type=='quan') { 
          #quantile bin limits 
          bin.lim<-quantile(pixval[is.finite(pixval)],c(max(c(0,lo-0.1)),min(c(hi+0.1,1))))
        } else { 
          #Absolute bins limits
          bin.lim<-c(lo*0.9,hi*1.1)
        }
      } else { 
        if (type=='quan') { 
          #quantile bin limits 
          quans<-quantile(pixval[is.finite(pixval)],c(max(c(0,lo-0.1)),min(c(hi+0.1,1))))
          bin.lim<-c(quans[1],quans[2])
        } else { 
          #Absolute bins limits
          bin.lim<-c(lo*0.9,hi*1.1)
        }
      } 
      if (bloom.bin) { 
        new.bin<-image.env$saturation
        if (do.sky.est) { new.bin<-new.bin-onesky$mu }
        if (grepl('SNR',bin.type)) { new.bin<-new.bin/onesky$sd }
        bin.lim[2]<-new.bin 
      }
      keep<-which(pixval >= bin.lim[1] & pixval <= bin.lim[2])
      point.sources<-point.sources[keep]
      nn.dist<-nn.dist[keep]
    }
    if (cutup) {
      if (quick.sky) { 
        message("Perfoming Fast Sky Estimation")
        skyest<-fast.sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,fit.gauss=fit.sky,saturate=image.env$saturation,
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
        skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,saturate=image.env$saturation,
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
        skyest<-fast.sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,fit.gauss=fit.sky,saturate=image.env$saturation,
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
        skyest<-sky.estimate(cat.x=cat.x,cat.y=cat.y,data.stamp.lims=data.stamp.lims,saturate=image.env$saturation,
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
      message("WARNING: cannot SNR bin without RMS estimate! 
              Using MAD of all pixels without a source centred on them (i.e. im[-(cat.x,cat.y)]) ") 
      tmprms<-mad(image.env$im[-(floor(cat.x)+nrow(image.env$im)*(floor(cat.y)-1))],na.rm=T)
      pixval<-pixval/tmprms
      pixval.all<-pixval.all/tmprms
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
    if (grepl('SNR',bin.type) & (do.sky.est | get.sky.rms)) { 
      new.bin<-new.bin/median(skyrms,na.rm=TRUE) 
    } else if (grepl('SNR',bin.type) & !(do.sky.est | get.sky.rms)) { 
      new.bin<-new.bin/tmprms 
    } 
    bin.lim<-c(bin.lim,new.bin) 
    n.bins<-n.bins+1
  }
  #If the pixval is outside the bins, skip it 
  keep <- which(pixval[point.sources] >= bin.zero & pixval[point.sources] <= max(bin.lim)) 
  point.sources<-point.sources[keep]
  nn.dist<-nn.dist[keep]
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
    }
  } else if (length(image.env$imm) > 1) {
    maskfrac<-rep(NA,length=length(cat.x))
    for (i in point.sources) {
      maskfrac[i]<-sum(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]]==0)
      maskfrac[i]<-maskfrac[i]/length(image.env$imm[ap.lims.mask.map[i,1]:ap.lims.mask.map[i,2],ap.lims.mask.map[i,3]:ap.lims.mask.map[i,4]])
    }
  }
  if (exists('maskfrac')) { 
    keep<-which(maskfrac[point.sources]<=mask.tolerance)
    point.sources<-point.sources[keep]
    nn.dist<-nn.dist[keep]
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
    blend.tolerance.tmp<-blend.tolerance
    mask.tolerance.tmp<-mask.tolerance
    radial.tolerance.tmp<-radial.tolerance
    #Check if we should iteratively can clean the sample 
    if (length(which(nn.dist[keep] > radial.tolerance*2.0 & maskfrac[point.sources[keep]] == 0 & blendfrac[point.sources[keep]] == 0)) < n.sources) { 
      radial.tolerance.use<-radial.tolerance
      mask.tolerance.use<-mask.tolerance
      blend.tolerance.use<-blend.tolerance
      while (length(keep) > n.sources) { 
        #Grow the distance tolerance
        radial.tolerance.use<-radial.tolerance.tmp
        radial.tolerance.tmp<-radial.tolerance.tmp+0.5
        #Reduce the masking tolerance
        mask.tolerance.use<-mask.tolerance.tmp
        mask.tolerance.tmp<-max(c(0,mask.tolerance.tmp-0.05))
        #Reduce the blending tolerance
        blend.tolerance.use<-blend.tolerance.tmp
        blend.tolerance.tmp<-max(c(0,blend.tolerance.tmp-0.05))
        #Calculate the new sample size
        keep<-which(bin[point.sources]==i & nn.dist>=radial.tolerance.tmp & maskfrac[point.sources]<=mask.tolerance.tmp & 
                    blendfrac[point.sources]<=blend.tolerance.tmp)
      }
      keep<-which(bin[point.sources]==i & nn.dist>=radial.tolerance.use & maskfrac[point.sources]<=mask.tolerance.use & 
                  blendfrac[point.sources]<=blend.tolerance.use)
      throw<-which(bin[point.sources]==i)
      throw<-throw[which(!throw%in%keep)]
      if (length(throw)!=0) { 
        point.sources<-point.sources[-throw]
        nn.dist<-nn.dist[-throw]
      }
    } else { 
      #There are more pure PSF sources than the number requested; just pick the first N.sources of them...
      keep<-which(bin[point.sources]==i & nn.dist > radial.tolerance*2.0 & maskfrac[point.sources] == 0 & blendfrac[point.sources] == 0) 
      if (n.sources+1 < length(keep)) { 
        throw<-sample(keep,size=length(keep)-n.sources,replace=FALSE)
        #throw<-keep[(n.sources+1):length(keep)]
        point.sources<-point.sources[-throw]
        nn.dist<-nn.dist[-throw]
      } else { 
        point.sources<-point.sources[keep]
        nn.dist<-nn.dist[keep]
      }
    }
  }
  if (plot & length(point.sources)>0) { 
    diagnostic<-TRUE
    if (diagnostic) { 
      #show all of the PSFs used in the stack
      nsamp<-length(point.sources)
      sample<-point.sources
    } else { 
      #show a sample of the PSFs used in the stack
      nsamp<-min(24,length(point.sources))
      sample=sample(point.sources,nsamp)
    }
    laymat<-(matrix(c(1:(ifelse(nsamp>12,12,nsamp)*4),rep(0,12*4-ifelse(nsamp>12,12,nsamp)*4)),ncol=8,byrow=T))
    layout(laymat)
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
    #Skip the source if there's a much brighter source that isn't this source 
    if (max(im[xlo:xup,ylo:yup])-sky > 5*pixval[i]) { 
      if (sqrt(sum(abs(which(im[xlo:xup,ylo:yup]==max(im[xlo:xup,ylo:yup],na.rm=T),arr.ind=T)-cbind(xc-xlo,yc-ylo))^2))>=3) { next } 
    }
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
        capture=magimage(im_psf[[bin[i]]]/weight[[bin[i]]],axes=FALSE)
        label('topleft',paste0(i,': PSF-stack (no mask)'),col='red',lwd=2)
      } else { 
        capture=magimage(im[xlo:xup,ylo:yup]*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': (1) raw (masked)'),col='red',lwd=2)
        magimage(im.cen*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': (2) recentred (masked)'),col='red',lwd=2)
        capture=magimage((im.cen-im[xlo:xup,ylo:yup])*mask[xmlo:xmup,ymlo:ymup],axes=FALSE)
        label('topleft',paste0(i,': Raw-centred (masked)'),col='red',lwd=2)
        capture=magimage(im_psf[[bin[i]]]/weight[[bin[i]]],axes=FALSE)
        label('topleft',paste0(i,': PSF-stack (masked)'),col='red',lwd=2)
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
