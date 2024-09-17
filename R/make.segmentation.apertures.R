make.segmentation.apertures <-
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

  #Pixel Ratio {{{
  pixrat<-seg.pix/arcsec.per.pix
  #}}}

  #Check if we want to resample.upres {{{
  if (!resample.aperture) { resample.iterations<-0 }
  #}}}

  #Create Aperture Masks {{{
  message("Creating Aperture Masks")
  #if ("EBImage" %in% rownames(installed.packages())) { 
  #  require(EBImage) 
  #  #Using the Affine transformation {{{
  #  s_mask<-foreach(slen=floor(stamplen/2),seglen=ceiling(stamplen/pixrat/2),x=seg.x,y=seg.y,
  #                  id=cat.id,xc=cat.x,yc=cat.y, 
  #                  .export=c("segmentation","pixrat"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
  #     #For each stamp, place down the relevant aperture {{{
  #     lims<-c(x+c(-seglen,seglen),y+c(-seglen,seglen))
  #     if (any(lims[1:2]>seg.astr$NAXIS[1]) | any(lims[3:4]>seg.astr$NAXIS[2]) | any(lims<1)) { 
  #       return(matrix(0,ncol=slen*2+1,nrow=slen*2+1))
  #     } else {
  #       seg.all<-segmentation[lims[1]:lims[2],lims[3]:lims[4]]
  #       seg.all[which(seg.all!=id)]<-0
  #       seg.all[which(seg.all==id)]<-1
  #       #place the aperture down on the new grid
  #       return(resize(seg.all,slen*2+1,slen*2+1))
  #     } 
  #     #}}}
  #  }
  #  #}}}
  #} else { 
    #Using 2D interpolation {{{
    message("WARNING: generating low-resolution apertures using the 2D interpolation, 
            instead of using the Affine Transformation. This can lead to _pathologically_ 
            bad behaviour for sources with widths at the ~pixel scale. For the Affine
            Transform, make sure that the EBImage package has been installed (from Bioconductor; 
            see the INSTALL document in the package git).")
    if (careful) {
      s_mask<-foreach(slen=zfloor(stamplen/2),seglen=ceiling(stamplen/pixrat/2),x=seg.x,y=seg.y,
                      id=cat.id,xc=cat.x,yc=cat.y, 
                      .export=c("segmentation","pixrat"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
         #For each stamp, place down the relevant aperture {{{
         lims<-c(x+c(-seglen,seglen),y+c(-seglen,seglen))
         if (any(lims[1:2]>seg.astr$NAXIS[1]) | any(lims[3:4]>seg.astr$NAXIS[2]) | any(lims<1)) { 
           return(matrix(0,ncol=slen*2+1,nrow=slen*2+1))
         } else { 
           seg.all<-segmentation[lims[1]:lims[2],lims[3]:lims[4]]
           seg.all[which(seg.all!=id)]<-0
           seg.all[which(seg.all==id)]<-1
           #place the aperture down on the new grid {{{
           #Make grid for aperture at old segmetation scale {{{
           seg.obj<-list(x=xc+seq(-seglen,seglen)*pixrat,y=yc+seq(-seglen,seglen)*pixrat,z=seg.all)
           #}}}
           #Make expanded grid of new pixel centres {{{
           expanded<-expand.grid(xc+seq(-slen,slen),yc+seq(-slen,slen))
           #}}}
           #Interpolate {{{
           seg.sing=matrix(interp.2d(expanded[,1], expanded[,2], seg.obj)[,3], ncol=slen*2+1,nrow=slen*2+1)
           seg.sing[slen+1,slen+1]<-1
           return(seg.sing)
           #}}}
           #}}}
         }
         #}}}
      }
    } else { 
      message(paste(seg.astr$NAXIS,collapse=' '))
      s_mask<-foreach(llim.x=ap.lims.data.map[,1],ulim.x=ap.lims.data.map[,2],llim.y=ap.lims.data.map[,3],ulim.y=ap.lims.data.map[,4],x=seg.x,y=seg.y,
                      id=cat.id,xc=cat.x,yc=cat.y, 
                      .export=c("segmentation","pixrat"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
         #For each stamp, place down the relevant aperture {{{
         #limits of stamp in segmentation image space 
         lims<-zfloor(c(x,x,y,y)+c(-1,1,-1,1)*(c(xc-llim.x,ulim.x-xc,yc-llim.y,ulim.y-yc)/pixrat))
         lims[1]<-max(c(lims[1],1))
         lims[3]<-max(c(lims[3],1))
         lims[2]<-min(c(lims[2],length(segmentation[,1])))
         lims[4]<-min(c(lims[4],length(segmentation[1,])))
         seg.all<-segmentation[lims[1]:lims[2],lims[3]:lims[4]]
         seg.all[which(seg.all!=id)]<-0
         seg.all[which(seg.all==id)]<-1
         #place the aperture down on the new grid {{{
         #Make grid for aperture at old segmetation scale {{{
         seg.obj<-list(x=(seq(lims[1],lims[2])-x)*pixrat,y=(seq(lims[3],lims[4])-y)*pixrat,z=seg.all)
         #}}}
         #Make expanded grid of new pixel centres {{{
         expanded<-expand.grid(seq(llim.x,ulim.x)-xc,seq(llim.y,ulim.y)-yc)
         #}}}
         #Interpolate {{{
         seg.sing=matrix(interp.2d(expanded[,1], expanded[,2], seg.obj)[,3], ncol=length(seq(llim.x,ulim.x)),nrow=length(seq(llim.y,ulim.y)))
         seg.sing[zfloor(xc-llim.x),zfloor(yc-llim.y)]<-1
         return(seg.sing)
         #}}}
         #}}}
         #}}}
      }
    } 
    #}}}
  #}
  message("Aperture Creation Complete")
  #}}}

  #Check that apertures do not cross image mask boundary {{{
  if ((cutup & length(mask.stamp)>1) | (!cutup & length(image.env$imm) > 1)) {
    #Check Mask stamps for Aperture Acceptance {{{
    message('Combining Aps with Mask Stamps')
    if (cutup) {
      sa_mask<-foreach(slen=stamplen, smask=s_mask,mmask=mask.stamp,mxl=ap.lims.mask.stamp[,1],mxh=ap.lims.mask.stamp[,2],myl=ap.lims.mask.stamp[,3],myh=ap.lims.mask.stamp[,4],.export=c("use.mask.lim","psf.filt"), .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
        #Check masking to determine if Aperture is acceptable {{{
        sumap<-sum(smask)
        check<-sum(mmask[mxl:mxh,myl:myh]*smask)/sumap
        if (sumap==0 || check<use.mask.lim) {
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
  if ((length(s_mask) > 0) && class(s_mask[[1]])[1]=="try-error") {
    stop(paste("Fatal Error in Foreach's Multi-Core Application.\n",
    "This seems to happen at random, but is more likely\n",
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
