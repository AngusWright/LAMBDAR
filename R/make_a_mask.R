make_a_mask <-
function (env=NULL, masks, fullmaskdim, outenv=NULL){
#Procedure makes the Full Mask of all apertures

  if (!quiet) { cat("Make_A_Mask   ") }
  message('--------------------------Make_A_Mask----------------------------------')

  # Load Parameter Space {{{
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #}}}

  #Setup sizes {{{
  npos<-length(id_g)
  #}}}

  #Initialise Array {{{
  a_mask<-array(0,dim=fullmaskdim)
  #}}}

  #Position each individual mask stamp at the correct position above the full mask, and add it on {{{
  for (i in 1:npos) {
    #Check dimensions match {{{
    if (any(dim(a_mask[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]])!=dim(masks[[i]]))) {
      #Don't match, so don't add it & Notify {{{
      message(paste("Dimensions of the",i,"th entries do not match. They are",dim(a_mask[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]])[1],"and",dim(masks[[i]])[1]))
      #}}}
    }else {
      #Add the stamp to the full mask {{{
      a_mask[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]]<-
      a_mask[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]]+masks[[i]]
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
  detach(env)
  return=a_mask
  #}}}
}
