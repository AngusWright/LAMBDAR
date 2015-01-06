make_a_mask <-
function (outenv=parent.env(environment()), masks, fullmaskdim, env=NULL){
#Procedure makes the Full Mask of all apertures

  if (!quiet) { cat("Make_A_Mask   ") }
  message('--------------------------Make_A_Mask----------------------------------')

  ## Load Parameter Space {{{
  if(!is.null(env)) {
    attach(env, warn.conflicts=FALSE)
  }
  if(is.null(outenv)&!is.null(env)) { outenv<-env }
  else if (is.null(outenv)) {
    warning("Output Environment cannot be NULL; using parent env")
    outenv<-parent.env(environment())
  }
  ##}}}

  #Setup sizes {{{
  npos<-length(id_g)
  #}}}

  #Initialise Array {{{
  a_mask<-array(0,dim=fullmaskdim)
  #}}}

  #Position each individual mask stamp at the correct position above the full mask, and add it on {{{
  for (i in 1:npos) {
    #Check dimensions match {{{
    if (any(dim(a_mask[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]])!=dim(masks[[i]]))) {
      #Don't match, so don't add it & Notify {{{
      message(paste("Dimensions of the ",i,"th entries do not match. They are [",paste(dim(a_mask[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]),collapse=','),"] and [",paste(dim(masks[[i]]),collapse=','),"]", sep=""))
      #}}}
    }else {
      #Add the stamp to the full mask {{{
      a_mask[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]<-
      a_mask[image_lims[i,1]:image_lims[i,2],image_lims[i,3]:image_lims[i,4]]+masks[[i]]
      #}}}
    }#}}}
  }#}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    message(paste("After assignment",round(length(which(is.na(a_mask)))/length(a_mask)*100,digits=2),"% of the a_mask matrix are NA"))
    message(paste("          -> NA values set to Zero"))
  }#}}}

  #Remove any NAs {{{
  a_mask[which(is.na(a_mask))]<-0
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) { message(paste("After assignment",round(length(which(a_mask>0))/length(a_mask)*100,digits=2),"% of the a_mask matrix is > 0")) }
  message('===========END=============Make_A_MASK=============END=================\n')
  #}}}

  #Return the Full Mask {{{
  if (!is.null(env)) { detatch(env) }
  return=a_mask
  #}}}
}
