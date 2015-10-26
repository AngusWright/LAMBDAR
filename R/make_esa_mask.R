make_esa_mask <-
function(outenv=parent.env(environment()), env=NULL,ObsParm,padGals,col.corr=0,confuse=FALSE){
#Procedure creates exponential galaxy profiles, using input parameters
#from the catalogue, and places them (in order) onto stamps
#Procedure is parallelised to allow scaleability

  message("Creating Simulated Galaxies from Input Catalogue")
  message("Assumptions are:\nCatalogue Galaxy Apertures are 2.5*Kron Apertures\nCatalogue Fluxweights are Jy fluxes for use in Simulation")
  if (!quiet) { cat('Make_ESA_Mask   ') }
  message('--------------------------Make_ESA_Mask---------------------------------')

  # Load Parameter Space
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  if(is.null(outenv)&!is.null(env)) { outenv<-env }
  else if (is.null(outenv)) {
    warning("Output Environment cannot be NULL; using parent env")
    outenv<-parent.env(environment())
  }
  #

  #Convert semi-major (Kron) axis in arcsec to Major Effective Radius in pixels
  Reff_pix<-(a_g)*(1/(2.5*kronrad(1)))/(asperpix)
  #

  #Aperture Axis Ratio in [0,1]
  axrat<-b_g/a_g
  axrat[which(!is.finite(axrat))]<-1
  #

  #Get input Magnitudes
  inputmags<- -2.5*log10(fluxweight)+8.9
  #

  #Remove Galaxies with non-finite input Magnitudes
  index<-which(is.finite(inputmags))
  inputmags<-inputmags[index]
  a_g<-a_g[index]
  b_g<-b_g[index]
  ra_g<-ra_g[index]
  dec_g<-dec_g[index]
  x_g<-x_g[index]
  y_g<-y_g[index]
  theta_g<-theta_g[index]
  Reff_pix<-Reff_pix[index]
  fluxweight<-fluxweight[index]
  axrat<-axrat[index]
  id_g<-id_g[index]
  contams<-contams[index]
  #

  if (padGals) {
    #Pad the input catalogue with stellar contaminants down to below the magnitude limit
    dens<-density(inputmags,bw=0.5,kernel='rect',na.rm=TRUE)
    mag.mode<-dens$x[which.max(dens$y)]
    #Determine Number Counts using Driver et al r-band relation
    ngal<-length(which(inputmags<=(mag.mode) & inputmags>=(mag.mode-0.5)))/(diff(range(ra_g))*diff(range(dec_g)))
    NDriver<-10^(0.38*(mag.mode-0.25+col.corr)-4.78)
    norm<-min(ngal/NDriver,1) #Normalise to represent correct area
    cat("Normalisation of Driver Relation is:",norm,"\n")
    #Magnitude Bin Values
    x<-seq(quantile(inputmags,c(0.01)),mag.mode+3.25,by=0.5)
    #Number of gals in each Bin
    ndraws<-norm*10^(0.38*(x+col.corr)-4.78)
    cat("Source Density: ",(sum(ndraws)+length(inputmags))/(diff(range(ra_g))*diff(range(dec_g))),' Objects/sqDegree\n')
    if (confuse) {
      #If we want confusion, increase this number so that
      #source density ~= psf size
      #Source Density / deg^2
      sourceDens<-(sum(ndraws)+length(inputmags))/(diff(range(ra_g))*diff(range(dec_g)))
      #Get source Density per PSF
      sDensPSF<-sourceDens*((gauss_fwhm_as/3600)^2)
      #If source density is less than probability of finding gal withing FWHM:
      if (sDensPSF<1.5) {
        add=0
        while (sDensPSF<1.5) {
          add=add+0.5
          #Get another bin of draws
          x<-seq(quantile(inputmags,c(0.01)),mag.mode+3.25+add,by=0.5)
          #Get new draws
          ndraws<-norm*10^(0.38*(x+col.corr)-4.78)
          #Source Density / deg^2
          sourceDens<-(sum(ndraws)+length(inputmags))/(diff(range(ra_g))*diff(range(dec_g)))
          #Get source Density per PSF
          sDensPSF<-sourceDens*((gauss_fwhm_as/3600)^2)
        }
        #ndraws<-ndraws*(1.5/sDensPSF)
      }
      cat("Confused Source Density: ",(sum(ndraws)+length(inputmags))/(diff(range(ra_g))*diff(range(dec_g))),' Objects/sqDegree\n')
    }
    #Do Magnitude Draws
    y<-foreach(draws=ndraws,lo=x-0.25,hi=x+0.25,.combine='c',.export='inputmags',.inorder=FALSE)%dopar%{
      n=draws-length(which(inputmags>lo & inputmags<hi))
      if (n<=0) {
        NULL
      } else {
        runif(draws-length(which(inputmags>lo & inputmags<hi)),min=lo,max=hi)
      }
    }
    y<-y+rnorm(length(y),0,0.2)
    photCount<-function(abmag,exp,area,Weff,lamEff) { return=10^((8.9-abmag)/2.5+7)*(1/1.51)*exp*area*(Weff/lamEff) }
    N=photCount(x,ObsParm$exp,ObsParm$area,ObsParm$Weff,ObsParm$lamEff)
    cat("Photon Counts in each x-bin:\nBin Centre: ",x,"\nN Phot: ",N,"\n")
    #y<-NULL
    #for( i in 1:length(x)) {
    #  draws=ndraws[i]
    #  lo=x[i]-0.25
    #  hi=x[i]+0.25
    #  y<-c(y,runif(draws-length(which(inputmags>lo & inputmags<hi)),min=lo,max=hi))
    #}
    mag.padgals<-y
    ra.padgals<-runif(length(mag.padgals),min=min(ra_g),max=max(ra_g))
    dec.padgals<-runif(length(mag.padgals),min=min(dec_g),max=max(dec_g))
    a_g.padgals<-a_g[runif(length(mag.padgals),min=1,max=length(a_g))]
    axrat.padgals<-runif(length(mag.padgals),min=min(axrat,na.rm=TRUE),max=max(axrat,na.rm=TRUE))
    b_g.padgals<-a_g.padgals*axrat.padgals
    Reff.padgals<-(a_g.padgals)*(1/(2.5*kronrad(1)))/(asperpix)
    theta.padgals<-runif(length(mag.padgals),min=min(theta_g),max=max(theta_g))
    #Add padding galaxies to the catalogue
    id_g<-c(id_g,9999000+1:length(mag.padgals))
    ra_g<-c(ra_g,ra.padgals)
    dec_g<-c(dec_g,dec.padgals)
    a_g<-c(a_g,a_g.padgals)
    b_g<-c(b_g,b_g.padgals)
    gamapos<-ad2xy(ra_g,dec_g,astr_struc)
    x_g<-gamapos[,1]
    y_g<-gamapos[,2]
    axrat<-c(axrat,axrat.padgals)
    theta_g<-c(theta_g,theta.padgals)
    Reff_pix<-c(Reff_pix,Reff.padgals)
    inputmags<-c(inputmags,mag.padgals)
    fluxweight.new<-10^c((8.9-inputmags)/2.5)
    if (any((zapsmall(fluxweight)-zapsmall(fluxweight.new[1:length(fluxweight)]))!=0,na.rm=TRUE)) {
      cat("\nFLUXWEIGHTS ARE NOT EQUAL\n")
    } else {
      fluxweight<-fluxweight.new
    }
    contams<-c(contams,rep(1,length(mag.padgals)))
  }

  #If needed, correct angluar coordinates
  #Details
  #Exponential Generation function uses N90E0 angular inputs
  ### This is OPPOSITE the aperture function used in LAMBDAR
  #If the input angles are not N90E0, correct for this offset
  if (!angoffset) { theta_off<-90-theta_g } else {theta_off<-theta_g}
  #Correct for any reversal of fits image
  if (astr_struc$CD[1,1]*astr_struc$CD[2,2]<0){ theta_off<-theta_off*-1 }


  #Reorder for faster parallelisation
  #if (MPIBackend) {
  #  nchild=getDoParWorkers()
  #} else {
    nchild<-ncores
  #}
  ind<-order(inputmags)
  matdim<-c(ceiling(length(inputmags)/nchild),nchild)
  ind<-c(ind,rep(NA,matdim[1]*matdim[2]-length(ind)))
  ind<-matrix(ind,ncol=matdim[2],nrow=matdim[1],byrow=T)
  ind<-as.numeric(ind)
  ind<-ind[which(!is.na(ind))]

  Reff_pix<-Reff_pix[ind]
  theta_off<-theta_off[ind]
  inputmags<-inputmags[ind]
  axrat<-axrat[ind]
  contams<-contams[ind]
  ra_g<-ra_g[ind]
  dec_g<-dec_g[ind]
  x_g<-x_g[ind]
  y_g<-y_g[ind]
  id_g<-id_g[ind]
  a_g<-a_g[ind]
  b_g<-b_g[ind]
  axrat<-axrat[ind]
  theta_g<-theta_g[ind]
  fluxweight<-fluxweight[ind]

  #Create Aperture Masks
  message("Creating Aperture Masks")
  psfsigma.pix<-psffwhm.pix/(2*sqrt(2*log(2)))
  if (confidence==1) {
    warning("Cannot have PSF Confidence = Unity when analytically deriving PSF; Using Maximum double value (1-1E-16)")
    confidence=1-1E-16
  }
  nsig<-sqrt(qchisq(confidence,df=2)) #Determine nsigma for desired confidence using chisq distribution
  psf.clip<-ceiling(nsig*psfsigma.pix)
  psf.clip<-psf.clip*2+1 # convert psf.clip from radius to diameter (and make sure it's odd)
  print(paste("PSFClip",psf.clip))
  #es_mask<-foreach(mag=inputmags, theta=theta_off,axr=axrat,Reff=Reff_pix,xdelt=(x_g%%1),ydelt=(y_g%%1), .export=c("psf.clip","psfsigma.pix","psffwhm.pix",'saturation','gain'), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
  es_mask<-NULL
  pb<-txtProgressBar(min=1,max=length(inputmags),style=3)
  for (i in 1:length(inputmags)) {
      setTxtProgressBar(pb,i)
      mag=inputmags[i]
      theta=theta_off[i]
      axr=axrat[i]
      Reff=Reff_pix[i]
      xdelt=(x_g[i]%%1)
      ydelt=(y_g[i]%%1)

      #Calculate Required number of Photons
      #Formula for number of photons given magnitude, exposure time (s),
      #telescope area (m^2), and filter Effective Width and Wavelength
      photCount<-function(abmag,exp,area,Weff,lamEff) { return=10^((8.9-abmag)/2.5+7)*(1/1.51)*exp*area*(Weff/lamEff) }
      #
      #Calculate N Photons given magnitude and SDSS r-band stats
      bn=1.678
      #Get N photons
      N=min(floor(photCount(mag,ObsParm$exp,ObsParm$area,ObsParm$Weff,ObsParm$lamEff)), 1E6)
      N<-1E6
      if (N==1E6) {
      #Use Analytic {{{
        stamplen<-(floor((ceiling(Reff*9)+ceiling(psf.clip))/2)*2+5)
        xy<-expand.grid(1:(stamplen*2),1:(stamplen*2))
        tempxy<-xy
        xy[,1]<-xy[,1]-stamplen-0.5
        xy[,2]<-xy[,2]-stamplen-0.5
        tempxy[,1]<-tempxy[,1]-stamplen-0.5-xdelt
        tempxy[,2]<-tempxy[,2]-stamplen-0.5-ydelt
        if (Reff!=0) {
          #Galaxy
          tempxy<-rotdata2d(tempxy[,1],tempxy[,2],90-theta)
          tempxy[,1]<-tempxy[,1]/axr
          I<-exp(-bn*(sqrt((xy[,1])^2+(xy[,2])^2)/Reff-1))
          I[which(I<max(I)*1E-5)]<-0
          im<-zapsmall(matrix(interp2D(tempxy[,1],tempxy[,2], list(x=1:(stamplen*2)-stamplen-0.0,y=1:(stamplen*2)-stamplen-0.0,z=matrix(zapsmall(I),stamplen*2,stamplen*2))), ncol=stamplen*2,nrow=stamplen*2,byrow=F))
          im<-im[(ceiling(stamplen*0.5)):(floor(stamplen*1.5)),(ceiling(stamplen*0.5)):(floor(stamplen*1.5))]
          #Convert to Jy
          im<-im/sum(im)*10^((8.9-mag)/2.5)
        }
        if (psffilt || Reff==0) {
        #Point Source
          ps<-matrix(dnorm(sqrt((xy[,1])^2+(xy[,2])^2),mean=0,sd=psfsigma.pix),stamplen*2,stamplen*2)
          ps<-ps[(ceiling(stamplen*0.5)):(floor(stamplen*1.5)),(ceiling(stamplen*0.5)):(floor(stamplen*1.5))]
          if (Reff==0) {
            #Convert to Jy
            im<-ps/sum(ps)*10^((8.9-mag)/2.5)
          } else {
            #Convert to Jy
            im<-convolvepsf(ps,im)
            im<-im/sum(im)*10^((8.9-mag)/2.5)
          }
        }
      #}}}
      } else {
        #Generate Profile with Shot noise & Convolution
        if (Reff!=0) {
          #Galaxy
          #Get N random photons from desired exponential profile
          tempr<-rexp(N, rate=bn/Reff)
          tempang<-runif(N, min=0, max=2*pi)

          #Convert from photon arrivals from Sperical to Cartesian Coords
          tempxy<-sph2car(long=tempang, lat=0, radius=tempr, deg=FALSE)

          #Stretch y-axis by ellipticity of profile
          tempxy[,2]<-(tempxy[,2])*axr

          #Rotate to desired theta
          tempxy<-rotdata2d(tempxy[,1],tempxy[,2],theta)

          #Filter through Gaussian PSF
          if (psffilt) {
            tempxy[1:length(tempxy)]<-tempxy[1:length(tempxy)]+rnorm(length(tempxy),mean=0,sd=psfsigma.pix)
          }
        } else {
          #Point Source; Use Gaussian with PSF FWHM
          tempxy<-matrix(rnorm(2*N,mean=0,sd=psfsigma.pix),ncol=2)
        }
        #Centroid correctly
        stamplen<-(floor((ceiling(Reff*9)+ceiling(psf.clip))/2)*2+5)
        tempx<-(tempxy[,1]+xdelt+stamplen/2)+0.5
        tempy<-(tempxy[,2]+ydelt+stamplen/2)+0.5
        k<-kde2d(tempx,tempy,lims=c(1,stamplen,1,stamplen),n=stamplen)
        im<-matrix(k$z,stamplen,stamplen)
        #for (i in 1:stamplen) { for (j in 1:stamplen) {  im[i,j]<-length(which(tempx==i & tempy==j)) } }

        #Convert from ADU to Jy
        im<-im/sum(im)*10^((8.9-mag)/2.5)
      }
      #return=im
      es_mask<-c(es_mask, list(im))
  }
  close(pb)
  message("Aperture Creation Complete")
  #

  #Check for Bad objects after generation
  es_mask<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
    #Check masking to determine if Aperture is acceptable
    if (any(is.na(smask))) {
      warning("Mask contains NA values. Setting those pixels to 0")
      smask[which(is.na(smask))]<-0
    }
    smask
  }
  stamplen<-foreach(smask=es_mask, .inorder=TRUE, .options.mpi=mpiopts, .combine='c') %dopar% { length(smask[,1]) }

  #Check that apertures do not cross image mask boundary
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    #Check Mask stamps for Aperture Acceptance
    message('Combining Aps with Mask Stamps')
    esa_mask<-foreach(slen=stamplen, smask=es_mask,mmask=imm_mask, .export="useMaskLim", .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
      #Check masking to determine if Aperture is acceptable
      if (any(is.na(smask))) {
        warning("Mask contains NA values. Setting those pixels to 0")
        smask[which(is.na(smask))]<-0
      }
      if (mean(mmask)<useMaskLim) {
        #Too much. Skip
        array(0, dim=c(slen,slen))
        #
      } else {
        #Not too much. Keep
        smask
        #
      }
      #
    }
    message('Mask Combine Finished.')
    #
  } else {
    #No image mask, all can be kept
    esa_mask<-es_mask
    #
  }#

  #Determine Stamp Limits in Image Space
  factor<-rep(0,length(esa_mask))
  for (i in 1:length(esa_mask)) { factor[i]<-((length(esa_mask[[i]][,1])-1)/2) }
  image_lims<-cbind(floor(x_g)-factor, floor(x_g)+factor, floor(y_g)-factor, floor(y_g)+factor)
  ind<-which(image_lims[,1]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][(image_lims[i,1]-1):0,]
    image_lims[i,1]<-1
  }
  ind<-which(image_lims[,3]<1)
  for (i in ind) {
    esa_mask[[i]]<-esa_mask[[i]][,(image_lims[i,3]-1):0]
    image_lims[i,3]<-1
  }
  ind<-which(image_lims[,2]>length(image.env$im[,1]))
  for (i in ind) {
    over<-image_lims[i,2]-length(image.env$im[,1])
    lim<-length(esa_mask[[i]][,1])-over
    esa_mask[[i]]<-esa_mask[[i]][1:lim,]
    image_lims[i,2]<-length(image.env$im[,1])
  }
  ind<-which(image_lims[,4]>length(image.env$im[1,]))
  for (i in ind) {
    over<-image_lims[i,4]-length(image.env$im[1,])
    lim<-length(esa_mask[[i]][1,])-over
    esa_mask[[i]]<-esa_mask[[i]][,1:lim]
    image_lims[i,4]<-length(image.env$im[1,])
  }
  #

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(esa_mask)))/length(esa_mask)*100,digits=2),"% of the esa_mask matrix are NA"))
  } #

  #Parse Parameter Space
  if (!is.null(env)) { detach(env) }
  assign("id_g",id_g,envir=outenv)
  assign("ra_g",ra_g,envir=outenv)
  assign("dec_g",dec_g,envir=outenv)
  assign("a_g",a_g,envir=outenv)
  assign("b_g",b_g,envir=outenv)
  assign("x_g",x_g,envir=outenv)
  assign("y_g",y_g,envir=outenv)
  assign("contams",contams,envir=outenv)
  assign("theta_g",theta_g,envir=outenv)
  assign("Reff_pix",Reff_pix,envir=outenv)
  assign("axrat",axrat,envir=outenv)
  assign("inputmags",inputmags,envir=outenv)
  assign("fluxweight",fluxweight,envir=outenv)
  assign("image_lims",image_lims,envir=outenv)
  #

  message('===========END============Make_ESA_MASK=============END=================\n')

  #Return array of apertures
  return=esa_mask
  #
}
