make_sfa_mask <-
function(env=NULL,sa_mask,fluxweightin=NULL, outenv=NULL) {
  #Make the Single Filtered Aperture mask
 
  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env } 
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))
  

  #Are we weighting the apertures?
  if (length(fluxweightin) > 0 ) { 
    if (!quiet) { cat("Make_WSFA_Mask     ")  }
    message('-------------------------Make_WSFA_Mask--------------------------------')
  } else { 
    if (!quiet) { cat("Make_SFA_Mask     ")  }
    message('-------------------------Make_SFA_Mask---------------------------------')
  }

  #Setup Sizes
  npos<-length(id_g)

  #Set Fluxweights
  #If no fluxweight - all weights == 1
  if (length(fluxweightin) == 0) { fluxweight<-array(1, dim=c(npos)) 
  #If single value suplied, all weights == value
  } else if (length(fluxweightin) == 1) { fluxweight<-array(fluxweightin, dim=c(npos)) 
  #If vector of values supplied, weights == input
  } else { fluxweight<-fluxweightin }

  #Check that number of Fluxweights == number of stamps
  if (!(length(fluxweight) == npos)) { sink(type="message") ; stop("Bad FluxWeightIn Array") } 

  #Perform calculations
  if ((length(which(fluxweight != 1)) != 0)&(nopsf)) {
    #If we are not filtering apertures with PSFs, and the fluxweights 
    #are not all unity, multiply each stamp by it's fluxweight
    sfa_mask<-foreach(i=1:npos, .inorder=TRUE)%dopar%{ sa_mask[[i]]*fluxweight[i]  }
  } else if (!(nopsf)) {
    #Else if we are filtering apertures with the PSFs, perform the convolution and 
    #multiply by the fluxweight simultaneously
    message('--------------------------Convolution----------------------------------')
    sfa_mask<-foreach(i=1:npos, xc=x_g, yc=y_g, .inorder=TRUE) %dopar% { 
      #Use square subset of psf array so that arrays are conformable
      #in convolution 
      if (length(sa_mask[[i]][,1])<length(psf[,1])){
        lims<-c(1,length(sa_mask[[i]][,1]))
        if (length(psf[lims[1]:lims[2],1])!=length(sa_mask[[i]][,1])) { print(paste(length(psf[lims[1]:lims[2],1]),length(sa_mask[[i]][,1]))) }
      } else if (length(sa_mask[[i]][,1])==length(psf[,1])){
        lims<-c(1,length(psf[,1]))
      } else {
        #sa_mask is LARGER than psf image
        sink(type="message")
        stop("Aperture Stamp is larger than PSF - convolution cannot be performed")
      }
      if (length(which(sa_mask[[i]]!=0))>1){
      #Convolve and weight
        ap=(convolvepsf(psf[lims[1]:lims[2],lims[1]:lims[2]],sa_mask[[i]])/sum(psf))*fluxweight[i]
        return=ap
      } else {
        #We have a point source - reinterpolate the PSF at point source XcenYcen, and return that
        if ((xc%%1!=0.5)|(yc%%1!=0.5)) {
          if (verbose) { message("Reinterpolating the Point Source PSF") }
          len<-length(psf[,1])
          #Make object of psf at old pixel centres
          psf_obj<-list(x=seq(1,len)+0.5, y=seq(1,len)+0.5,z=psf)
          #Make expanded grid of new pixel centres
          expanded<-expand.grid(seq(1,len),seq(1,len))
          xnew<-expanded[,1]+xc%%1
          ynew<-expanded[,2]+yc%%1
          #Call Interpolation
          ap=matrix(interp2D(xnew, ynew, psf_obj)*fluxweight[i], ncol=dim(psf)[1])[lims[1]:lims[2],lims[1]:lims[2]]
          #return(ap/max(ap,na.rm=TRUE))
          return=ap
        } else { 
          ap=psf[lims[1]:lims[2],lims[1]:lims[2]]
          #return(ap/max(ap,na.rm=TRUE))
          return=ap
	}
      }  
    }
    message('===========END============Convolution==============END=================\n')
  } else {
    #If we are not filtering by psf, and all fluxweights are equal to unity, 
    #then we would not have entered this function in the first place. If we 
    #arrive here, something has gone wrong with our fluxweights or global variables.
    sink(type="message")
    stop("Bad variables and/or fluxweights in production of filtered apertures.")
  }

  #Check for production of NA/NaN/Infs in convolution
  if (length(which(is.na(sfa_mask)))>0) {sink(type="message") ; stop("NAs produced in Convolution")} 
  if (length(fluxweightin) > 0 ) { message('===========END===========Make_WSFA_Mask============END=================\n')
  } else { message('===========END===========Make_SFA_MASK=============END=================\n') }
  #Return array of Stamps
  detach(env)
  sfa_mask
}
