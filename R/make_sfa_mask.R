make_sfa_mask <-
function(env=NULL,sa_mask,fluxweightin=NULL, outenv=NULL) {
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
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #}}}

  #Setup Sizes {{{
  npos<-length(id_g)
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
  #}}}

  #Check that number of Fluxweights == number of stamps {{{
  if (!(length(fluxweight) == npos)) { sink(type="message") ; stop("Bad FluxWeightIn Array") }
  #}}}

  #Perform calculations {{{
  if ((length(which(fluxweight != 1)) != 0)&(nopsf)) {
    #No Convolution, need Fluxweighting {{{
    #Details {{{
    #If we are not filtering apertures with PSFs, and the fluxweights
    #are not all unity, multiply each stamp by it's fluxweight }}}
    sfa_mask<-foreach(i=1:npos, .inorder=TRUE)%dopar%{ sa_mask[[i]]*fluxweight[i]  }
    #}}}
  } else if (!(nopsf)) {
    #Convolution with PSF & Fluxweighting {{{
    #Details {{{
    #Else if we are filtering apertures with the PSFs, perform the convolution and
    #multiply by the fluxweight simultaneously }}}
    message('--------------------------Convolution----------------------------------')
    sfa_mask<-foreach(i=1:npos, xc=x_g, yc=y_g, .inorder=TRUE) %dopar% {
      #Use subset of psf conformable with current aperture stamp {{{
      if (length(sa_mask[[i]][,1])<length(psf[,1])){
        #Aperture Stamp is smaller than PSF stamp {{{
        #limits from PSF peak - 1/2 stampwidth {{{
        centre<-as.numeric(which(psf==max(psf), arr.ind=TRUE))
        sa.len<-length(sa_mask[[i]][,1])
        if (sa.len%%2) { sa.len<-sa.len-1 }
        delta<-(sa.len)*c(-0.5,+0.5)
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
        if (length(psf[lims[1]:lims[3],1])!=length(sa_mask[[i]][,1])) {
          stop(paste("PSF and Aperture Arrays are non-conformable. Lengths are:",length(psf[lims[1]:lims[3],1]),length(sa_mask[[i]][,1])))
        }
        #}}}
        #}}}
      } else if (length(sa_mask[[i]][,1])==length(psf[,1])){
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
      if (length(which(sa_mask[[i]]!=0))>1){
        #Aperture {{{
        #Diagnostic {{{
        if (diagnostic) { message(paste("Doing convolution of PSF with Aperture", i)) }
        #}}}
        #Convolve {{{
        ap<-(convolvepsf(psf[lims[1]:lims[3],lims[2]:lims[4]],sa_mask[[i]]))
        #}}}
        #Check for Negatives produced by convolution {{{
         if (length(which(ap <0))>0) {
          #If Negatives produced, rezap with required 'digits' level {{{
          message(paste("Negatives Produced in PSF convolution with Aperture",i))
          zapdig<-floor(-log10(abs(min(ap))))
          ap<-(convolvepsf(psf[lims[1]:lims[3],lims[2]:lims[4]],sa_mask[[i]],zapdig=zapdig))
          message(paste("Attempting rezap with zapdigit",zapdig))
          if (length(which(ap<0))>0) {
            #if still negatives, error {{{
            stop("Rezap unable to remove negatives in convolution")
            #}}}
          }
          #}}}
        }
        #}}}
        #Normalise and Weight apertures {{{
        ap<-ap/max(ap, na.rm=TRUE)*fluxweight[i]
        #}}}
        #Return {{{
        return=ap
        #}}}
        #}}}
      } else {
      #Point Source {{{
        if ((xc%%1!=0.5)|(yc%%1!=0.5)) {
        #Reinterpolate the PSF at point source XcenYcen {{{
          len<-length(psf[,1])
          #Make grid for psf at old pixel centres {{{
          psf_obj<-list(x=seq(1,len)+0.5, y=seq(1,len)+0.5,z=psf)
          #}}}
          #Make expanded grid of new pixel centres {{{
          expanded<-expand.grid(seq(1,len),seq(1,len))
          xnew<-expanded[,1]+xc%%1
          ynew<-expanded[,2]+yc%%1
          #}}}
          #Interpolate {{{
          ap<-matrix(interp2D(xnew, ynew, psf_obj), ncol=length(psf[,1]))
          #}}}
          #Make integral of interpolated source same as psf integral {{{
          ap<-(ap*(sum(psf)/sum(ap))*fluxweight[i])[lims[1]:lims[3],lims[2]:lims[4]]
          #}}}
          #Return {{{
          return=ap
          #}}}
          #}}}
        } else {
          #Assign PSF as object aperture {{{
          if (diagnostic) { message(paste("Using the Uninterpolated PSF for Point Source Aperture",i)) }
          ap<-psf[lims[1]:lims[3],lims[2]:lims[4]]*fluxweight[i]
          return=ap
          #}}}
	      }
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
    #arrive here, something has gone wrong with our fluxweights or global variables. }}}
    sink(type="message")
    stop("Bad variables and/or fluxweights in production of filtered apertures.")
    #}}}
  }
  #}}}

  #Check for production of NA/NaN/Infs in convolution {{{
  if (diagnostic) { for (i in length(sfa_mask)) { if (length(which(is.na(as.numeric(sfa_mask[[i]]))))>0) {sink(type="message") ; stop("NAs produced in Convolution")} } }
  if (length(fluxweightin) > 0 ) { message('===========END===========Make_WSFA_Mask============END=================\n')
  } else { message('===========END===========Make_SFA_MASK=============END=================\n') }
  #}}}

  #Return array of Stamps {{{
  detach(env)
  return=sfa_mask
  #}}}
}
