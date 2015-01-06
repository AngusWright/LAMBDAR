sourcesubtraction <-
function(image,sfa_models,stamp_lims,fluxes,outputfile,outputheader,beamarea,insidemask,diagnostic=FALSE,verbose=FALSE) {
  # Procedure subtracts stamps from the input image.

  message('------------------------Source_Subtraction-----------------------------')
  #Initialise Residual Image {{{
  resid_im<-image  #image is currently in Jy/pix
  npos<-length(sfa_models)
  if (length(fluxes)==1) { fluxes<-replicate(npos, fluxes) }
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    if (length(which(is.na(image)))!=0) {        message(paste("Input Parameter 'image' contains NA elements")) }
    if (length(which(is.na(sfa_models)))!=0) {   message(paste("Input Parameter 'sfa_models' contains NA elements")) }
    if (length(which(is.na(stamp_lims)))!=0) {   message(paste("Input Parameter 'stamp_lims' contains NA elements")) }
    if (length(which(is.na(fluxes)))!=0) {       message(paste("Input Parameter 'fluxes' contains NA elements")) }
    if (length(which(is.na(outputfile)))!=0) {   message(paste("Input Parameter 'outputfile' contains NA elements")) }
    if (length(which(is.na(outputheader)))!=0) { message(paste("Input Parameter 'outputheader' contains NA elements")) }
    if (length(which(is.na(beamarea)))!=0) {     message(paste("Input Parameter 'beamarea' contains NA elements")) }
    if (length(which(is.na(insidemask)))!=0) {   message(paste("Input Parameter 'insidemask' contains NA elements")) }
    message(paste("At initialisation",round(length(which(is.na(resid_im)))/length(resid_im)*100,digits=2),"% of the resid_im matrix are NA"))
  }
  #}}}

  #Get MPI options {{{
  chunkSize=ceiling(npos/getDoParWorkers())
  mpiopts<-list(chunkSize=chunkSize)
  #}}}

  #Perform Subtraction {{{
  modelsource<-foreach (i=1:npos, flux=fluxes, inside=insidemask, model=sfa_models, .inorder=TRUE, .options.mpi=mpiopts) %dopar% {
    if ((sum(abs(model))==0)|(!(inside))) {
      #If stamp is empty, skip it {{{
      return=NULL
      #}}}
    } else {
      #Otherwise, calculate normalised source model {{{
      return=model*flux/sum(model)
      #}}}
    }
  }
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    #Check Limits
    if (!(is.finite(stamp_lims))) { sink(type="message") ; stop("Bad limits in Source Subtraction") }
  }
  #}}}
  #Subtract model from image {{{
  for (i in 1:npos) {
    if (!is.null(modelsource[[i]])) {
      resid_im[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]]<-
      resid_im[stamp_lims[i,1]:stamp_lims[i,2],stamp_lims[i,3]:stamp_lims[i,4]]-modelsource[[i]]
    }
  }
  #}}}
  #Convert image units back into Jy/bm {{{
  resid_im<-resid_im*beamarea
  #}}}
  if (diagnostic) { message(paste("After assignment",round(length(which(is.na(resid_im)))/length(resid_im)*100,digits=2),"% of the resid_im matrix are NA")) }
  #Output Residual Image {{{
  writefitsout(outputfile,resid_im,outputheader, nochange=TRUE,verbose=verbose,diagnostic=diagnostic)
  #}}}
  message('=========END==========Source_Subtraction===========END=================\n')
  return=NULL
}
