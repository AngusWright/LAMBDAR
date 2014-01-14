make_sa_mask <-
function(env=NULL, outenv=NULL){
  # Procedure creates apertures, using input parameters
  # from the catalogue, and places them (in order) onto stamps
  # Procedure is parallelised to allow scaleability

  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))

  if (!quiet) { cat('Make_SA_Mask   ') }
  message('--------------------------Make_SA_Mask---------------------------------')

  #Setup sizes
  npos<-length(ra_g)
  imxpixmax<-length(image.env$im[,1])
  imypixmax<-length(image.env$im[1,])

  #-----Diagnostic-----#
  #Check to make sure pixel values are finite
  if (length(which(!(is.finite(x_p) | !(is.finite(y_p)))))!=0) {
    sink(type="message") ; stop(paste("x_p[",which(!(is.finite(x_p))),"] or y_p[",which(!(is.finite(x_p))),"] are NA/NaN/Inf in  Make_SA_Mask()", sep=""))
  }

  message("Setting Stamp Limits")
  #Stampwidths: aperture width multiplied by a buffer, plus the
  #PSF WIDHT determined by the desired confidence. Total axis length MUST be odd
  stamplen<-(floor((ceiling(defbuff*a_g*2/asperpix)+ceiling(psf.clip))/2)*2+1)
  #If stamplen==1, then make the stamplen == smallest nonzero aperture by default
  #This happens when we have a point source and are not convolving with PSF
  #If ALL apertures are point sources, then this will have length 3 
  if      (all(stamplen==1)) { stamplen[which(stamplen==1)]<-3 }
  else if (any(stamplen==1)) { stamplen[which(stamplen==1)]<-min(a_g[which(a_g>0)],na.rm=TRUE)*2+3 }
  #Calculate Stamp limits in image-pixel space; parallelised
  stamp_lims_list<-foreach(i=1:npos,width=floor(stamplen/2),.inorder=TRUE) %dopar% {
    #Determine stamp limits in sourcemask array-space
    xl<-x_p[i]-width+1
    xu<-x_p[i]+width+1
    yl<-y_p[i]-width+1
    yu<-y_p[i]+width+1
    return(cbind(xl,xu,yl,yu))
  }
  #convert stamp limits back into a [4,npos] array
  stamp_lims<-do.call(rbind, stamp_lims_list)

  #Check that stamp limits are within image
  message("Removing Stamps that lie outside the image limits")
  catlen<-length(x_g)
  insidemask<-!((stamp_lims[,1]<1)|(stamp_lims[,2]>length(image.env$im[,1]))|(stamp_lims[,3]<1)|(stamp_lims[,4]>length(image.env$im[1,])))
  #remove any apertures whos stamps would cross the boundary
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
  stamp_lims<-stamp_lims[which(insidemask),]
  if (filtcontam) { contams<-contams[which(insidemask)] }
  insidemask<-insidemask[which(insidemask)]
  if (verbose) { message(paste("There are",length(x_g),"supplied objects have stamps entirely inside the image (",
                                round(((catlen-length(x_g))/catlen)*100, digits=2),"% of supplied lie across the edge of the image )")) }

  #-----Diagnostic-----#
  if (diagnostic) {
    #Check Stamp Lengths are correct
    check<-foreach(i=1:length(stamplen), .combine='c') %dopar% { (stamplen[i]==diff(stamp_lims[i,1],stamp_lims[i,2]))&(stamplen[i]==diff(stamp_lims[i,1],stamp_lims[i,2])) }
    message(paste("There are",length(which(!check)),"stamps with lengths (from limits) not equal to their desired stamp length"))
    #Check for errors in stamp creation
    if (length(which(!(is.finite(stamp_lims))))!=0) {
      sink(type="message") ; stop(paste("stamp_lims[",which(!(is.finite(stamp_lims))),"] are NA/NaN/Inf in Make_SA_Mask()", sep=""))
    }
  }

  #Convert semi-major axis in arcsec to semi-major axis in pixels
  a_g_pix<-a_g/(asperpix)

  #Correct for Angluar Coordinate System - aperture function uses N0E90 angular inputs
  #If the input angles are not N0E90, correct for this offset
  if (angoffset) { theta_off<-90-theta_g } else {theta_off<-theta_g}
  #Correct for any reversal of fits image
  if (astr_struc$CD[1,1]*astr_struc$CD[2,2]<0){ theta_off<-theta_off*-1 }
  #Aperture Axis Ratio in [0,1]
  axrat<-b_g/a_g

  #-----Diagnostic-----#
  if (diagnostic) {
    sink(file=sinkfile, type='output')
    ls.str(envir=env)
    sink(type='output')
  }

  message("Creating Aperture Masks")
  ##Put apertures on stamps masks
  ##Make grid for apertures
  #Determine Aperture-Pixel offset : (actual pos - pixel pos) - halfpix width
  s_mask<-foreach(slen=stamplen,axa=theta_off,axr=axrat,maj=a_g_pix,xdelt=(x_g%%1),ydelt=(y_g%%1),.inorder=TRUE) %dopar% {
      #Setup Grid for use in aperture function
      grid<-expand.grid(seq((1.5),(slen+0.5), by=1),seq((1.5),(slen+0.5), by=1))
      if (any(is.na(grid))){ stop(paste("NAs produced in Expand Grid. Stamplen=",stamplen)) }
      #For each stamp, place down the relevant aperture
      expanded<-iterapint(x=grid[,1],y=grid[,2],xstep=1,ystep=1,xcen=ceiling(slen/2)+xdelt,ycen=ceiling(slen/2)+ydelt,axang=axa,
                          axrat=axr,majax=maj,upres=upres,itersteps=itersteps, peakscale=TRUE)
      matrix(expanded[,3],ncol=slen,byrow=TRUE)
  }
  #covert list of stamps into array of stamps
  message("Aperture Creation Complete")

  ##For each catalogue entry, make a stamp of the mask at that position
  #If the mask is unity everywhere (or is of length 1)
  #we don't need to use the imm mask at all; so skip
  if (((length(image.env$imm)!=1)&(length(which(image.env$imm!=1))!=0))) {
    message('Combining with Mask')
    #create mask stamps
    sa_mask<-foreach(i=1:npos,stampsize=stamplen,xlo=stamp_lims[,1],xup=stamp_lims[,2], ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE) %dopar% {
      if ((xlo < 1) | (xup > imxpixmax) | (ylo < 1) | (yup > imypixmax)) {
        #Stamp Extends beyond sourcemask limits - set this mask to zero
        array(0, dim=c(stampsize,stampsize))
      } else {
        if (mean(image.env$imm[xlo:xup,ylo:yup])<useMaskLim) {
          array(0, dim=c(stampsize,stampsize))
        } else {
          s_mask[[i]]
        }
      }
    }
    #sa_mask<-array(unlist(sa_mask), dim=c(dim(sa_mask[[1]]),length(sa_mask)))
    message('Mask Combine Finished.')
  } else {
    sa_mask<-s_mask
  }

  #-----Diagnostic-----#
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(sa_mask)))/length(sa_mask)*100,digits=2),"% of the sa_mask matrix are NA"))
  }
  message('===========END============Make_SA_MASK=============END=================\n')
  #Parse Parameter Space
  detach(env)
  assign("x_p",x_p,envir=outenv)
  assign("y_p",y_p,envir=outenv)
  assign("x_g",x_g,envir=outenv)
  assign("y_g",y_g,envir=outenv)
  assign("id_g",id_g,envir=outenv)
  assign("ra_g",ra_g,envir=outenv)
  assign("dec_g",dec_g,envir=outenv)
  assign("theta_g",theta_g,envir=outenv)
  assign("a_g",a_g,envir=outenv)
  assign("a_g_pix",a_g_pix,envir=outenv)
  assign("b_g",b_g,envir=outenv)
  assign("stamplen",stamplen,envir=outenv)
  assign("stamp_lims",stamp_lims,envir=outenv)
  if (filtcontam) { assign("contams",contams,envir=outenv) }
  assign("insidemask",insidemask,envir=outenv)
  #Return array of apertures
  sa_mask
}
