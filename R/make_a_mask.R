make.aperture.map <-
function (outenv=parent.env(environment()), masks, fullmaskdim, env=NULL, subs=NULL){
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
  if (is.null(subs)) {
    npos<-length(cat.id)
    subs<-1:npos
  } else {
    npos<-length(subs)
  }
  #}}}

  #Initialise Array {{{
  a_mask<-array(0,dim=fullmaskdim)
  #}}}

  #Position each individual mask stamp at the correct position above the full mask, and add it on {{{
  for (i in subs) {
    #Check dimensions match {{{
    if (any(dim(a_mask[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]])!=dim(masks[[i]]))) {
      #Don't match, so don't add it & Notify {{{
      message(paste("Dimensions of the ",i,"th entries do not match. They are [",paste(dim(a_mask[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]),collapse=','),"] and [",paste(dim(masks[[i]]),collapse=','),"]", sep=""))
      #}}}
    }else {
      #Add the stamp to the full mask {{{
      a_mask[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]<-
      a_mask[ap.lims.data.map[i,1]:ap.lims.data.map[i,2],ap.lims.data.map[i,3]:ap.lims.data.map[i,4]]+masks[[i]]
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
  if (!is.null(env)) { detach(env) }
  return=a_mask
  #}}}
}
