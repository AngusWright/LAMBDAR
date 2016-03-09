make_sa_mask <-
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
  a_g_pix<-a_g/(asperpix)
  #}}}

  #If needed, correct angluar coordinates {{{
  #Details {{{
  #Aperture function uses N0E90 angular inputs
  #If the input angles are not N0E90, correct for this offset }}}
  if (angoffset) { theta_off<-90-theta_g } else {theta_off<-theta_g}
  #Correct for any reversal of fits image
  if (astr_struc$CD[1,1]*astr_struc$CD[2,2]<0){ theta_off<-theta_off*-1 }
  #}}}

  #Aperture Axis Ratio in [0,1] {{{
  axrat<-b_g/a_g
  #}}}

  #Check if we want to upres {{{
  if (!resampleaperture) { itersteps<-0 }
  #}}}

  #Create Aperture Masks {{{
  message("Creating Aperture Masks")
  s_mask<-foreach(slen=stamplen,axa=theta_off,axr=axrat,maj=a_g_pix,xdelt=(x_g%%1),ydelt=(y_g%%1), .export=c("itersteps","upres"), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
      #Setup Grid for use in aperture generation {{{
      grid<-expand.grid(seq((1.5),(slen+0.5), by=1),seq((1.5),(slen+0.5), by=1))
      if (any(is.na(grid))){ stop(paste("NAs produced in Expand Grid. Stamplen=",slen)) }
      #}}}
      #For each stamp, place down the relevant aperture {{{
      expanded<-iterapint(x=grid[,1],y=grid[,2],xstep=1,ystep=1,xcen=ceiling(slen/2)+0.5+xdelt,ycen=ceiling(slen/2)+0.5+ydelt,axang=axa,
                          axrat=axr,majax=maj,upres=upres,itersteps=ifelse(maj==0,0,itersteps), peakscale=TRUE)
      matrix(expanded[,3],ncol=slen,byrow=FALSE)
      #}}}
  }
  message("Aperture Creation Complete")
  #}}}

  #Check that apertures do not cross image mask boundary {{{
  if ((cutup & length(imm_mask)>1) | (!cutup & length(image.env$imm) > 1)) {
    #Check Mask stamps for Aperture Acceptance {{{
    message('Combining Aps with Mask Stamps')
    if (cutup) {
      sa_mask<-foreach(slen=stamplen, smask=s_mask,mmask=imm_mask,mxl=mstamp_lims[,1],mxh=mstamp_lims[,2],myl=mstamp_lims[,3],myh=mstamp_lims[,4],.export=c("useMaskLim","psffilt"), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(mmask[mxl:mxh,myl:myh]*smask)/sum(smask)
        if (check<useMaskLim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          if (psffilt) {
            smask
          } else {
            smask*mmask[mxl:mxh,myl:myh]
          }
          #}}}
        }
        #}}}
      }
    } else {
      sa_mask<-foreach(slen=stamplen, smask=s_mask,mxl=mask_lims[,1],mxh=mask_lims[,2],myl=mask_lims[,3],myh=mask_lims[,4],.export=c("useMaskLim","image.env","psffilt"), .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        check<-sum(image.env$imm[mxl:mxh,myl:myh]*smask)/sum(smask)
        if (check<useMaskLim) {
          #Too much. Skip {{{
          array(0, dim=c(slen,slen))
          #}}}
        } else {
          #Not too much. Keep {{{
          if (psffilt) {
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
  assign("theta_off",theta_off,envir=outenv)
  assign("theta_g",theta_g,envir=outenv)
  assign("a_g_pix",a_g_pix,envir=outenv)
  #}}}

  message('===========END============Make_SA_MASK=============END=================\n')

  #Return array of apertures {{{
  return=sa_mask
  #}}}
}
