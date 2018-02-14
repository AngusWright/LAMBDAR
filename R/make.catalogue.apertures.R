make.catalogue.apertures <-
function(outenv=parent.env(environment()), env=NULL){
#Procedure creates apertures, using input parameters
#from the catalogue, and places them (in order) onto stamps
#Procedure is parallelised to allow scaleability

  if (!quiet) { cat('Make_SA_Mask   ') }
  message('--------------------------Make_SA_Mask---------------------------------')

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

  #Convert semi-major axis in arcsec to semi-major axis in pixels {{{
  cat.a.pix<-cat.a/(arcsec.per.pix)
  #}}}

  #If needed, correct angluar coordinates {{{
  #Details {{{
  #Aperture function uses N0E90 angular inputs
  #If the input angles are not N0E90, correct for this offset }}}
  if (ang.offset) { theta.offset<-90-cat.theta } else {theta.offset<-cat.theta}
  #Correct for any reversal of fits image
  if (astr.struc$CD[1,1]*astr.struc$CD[2,2]<0){ theta.offset<-theta.offset*-1 }
  #}}}

  #Aperture Axis Ratio in [0,1] {{{
  axrat<-cat.b/cat.a
  #}}}

  #Check if we want to resample.upres {{{
  if (!resample.aperture) { resample.iterations<-0 }
  #}}}

  #Create Aperture Masks {{{
  message("Creating Aperture Masks")
  if (!careful) { 
    s_mask<-foreach(slen=stamplen,axa=theta.offset,axr=axrat,maj=cat.a.pix,stamplim.x=ap.lims.data.map[,1],stamplim.y=ap.lims.data.map[,3],cat.x=cat.x,cat.y=cat.y, .export=c("resample.iterations","resample.upres"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Setup Grid for use in aperture generation {{{
        grid<-expand.grid(seq((1.5),(slen+0.5), by=1),seq((1.5),(slen+0.5), by=1))
        if (any(is.na(grid))){ stop(paste("NAs produced in Expand Grid. Stamplen=",slen)) }
        #}}}
        #For each stamp, place down the relevant aperture {{{
        expanded<-generate.aperture(x=grid[,1],y=grid[,2],xstep=1,ystep=1,xcen=cat.x-stamplim.x,ycen=cat.y-stamplim.y,axang=axa,
                            axrat=axr,majax=maj,resample.upres=resample.upres,resample.iterations=ifelse(maj==0,0,resample.iterations), peakscale=TRUE)
        matrix(expanded[,3],ncol=slen,byrow=FALSE)
        #}}}
    }
  } else {
    s_mask<-foreach(slen=stamplen,axa=theta.offset,axr=axrat,maj=cat.a.pix,xdelt=(cat.x%%1),ydelt=(cat.y%%1), .export=c("resample.iterations","resample.upres"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Setup Grid for use in aperture generation {{{
        grid<-expand.grid(seq((1.5),(slen+0.5), by=1),seq((1.5),(slen+0.5), by=1))
        if (any(is.na(grid))){ stop(paste("NAs produced in Expand Grid. Stamplen=",slen)) }
        #}}}
        #For each stamp, place down the relevant aperture {{{
        expanded<-generate.aperture(x=grid[,1],y=grid[,2],xstep=1,ystep=1,xcen=ceiling(slen/2)+0.5+xdelt,ycen=ceiling(slen/2)+0.5+ydelt,axang=axa,
                            axrat=axr,majax=maj,resample.upres=resample.upres,resample.iterations=ifelse(maj==0,0,resample.iterations), peakscale=TRUE)
        matrix(expanded[,3],ncol=slen,byrow=FALSE)
        #}}}
    }
  } 
  message("Aperture Creation Complete")
  #}}}

  #Check that apertures do not cross image mask boundary {{{
  if ((cutup & length(mask.stamp)==length(s_mask)) | (!cutup & length(image.env$imm) > 1)) {
    #Check Mask stamps for Aperture Acceptance {{{
    message('Combining Aps with Mask Stamps')
    if (cutup) {
      sa_mask<-foreach(slen=stamplen, smask=s_mask,mmask=mask.stamp,mxl=ap.lims.mask.stamp[,1],mxh=ap.lims.mask.stamp[,2],myl=ap.lims.mask.stamp[,3],myh=ap.lims.mask.stamp[,4],.export=c("use.mask.lim","psf.filt"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(mmask[mxl:mxh,myl:myh]*smask)/sum(smask)
        if (check<use.mask.lim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          if (psf.filt) {
            smask
          } else {
            smask*mmask[mxl:mxh,myl:myh]
          }
          #}}}
        }
        #}}}
      }
    } else {
      sa_mask<-foreach(slen=stamplen, smask=s_mask,mxl=ap.lims.mask.map[,1],mxh=ap.lims.mask.map[,2],myl=ap.lims.mask.map[,3],myh=ap.lims.mask.map[,4],.export=c("use.mask.lim","image.env","psf.filt"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(image.env$imm[mxl:mxh,myl:myh]*smask)/sum(smask)
        if (check<use.mask.lim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          if (psf.filt) {
            smask
          } else {
            smask*image.env$imm[mxl:mxh,myl:myh]
          }
          #}}}
        }
        #}}}
      }
    }
    message('Mask Combine Finished.')
    #}}}
  } else {
    #No image mask, all can be kept {{{
    sa_mask<-s_mask
    #}}}
  }#}}}

  #Check for Foreach errors {{{
  if ((length(s_mask) > 0) && class(s_mask[[1]])=="try-error") {
    stop(paste("Fatal Error in Foreach's Multi-Core Application.\n",
    "This seems to at random, but is more likely\n",
    "to happen with increasing numbers of objects\n",
    "per thread. Try increasing the number of\n",
    "doParallel cores being used by the program.",sep=""))
  }
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(sa_mask)))/length(sa_mask)*100,digits=2),"% of the sa_mask matrix are NA"))
  } #}}}

  #Parse Parameter Space {{{
  if (!is.null(env)) { detach(env) }
  assign("theta.offset",theta.offset,envir=outenv)
  assign("cat.theta",cat.theta,envir=outenv)
  assign("cat.a.pix",cat.a.pix,envir=outenv)
  #}}}

  message('===========END============Make_SA_MASK=============END=================\n')

  #Return array of apertures {{{
  return=sa_mask
  #}}}
}
