source.subtraction <-
function(image,sfa_models,ap.lims.data.stamp,fluxes,outputfile,outputheader,beamarea,inside.mask,diagnostic=FALSE,verbose=FALSE,fluxmap.outputfile="fluxmap.fits") {
  # Procedure subtracts stamps from the input image.

  message('------------------------Source_Subtraction-----------------------------')
  #Initialise Residual Image {{{
  resid_im<-image  #image is currently in Jy/pix
  flux.map<-array(0,dim(image))
  npos<-length(sfa_models)
  if (length(fluxes)==1) { fluxes<-replicate(npos, fluxes) }
  #}}}

  #-----Diagnostic-----# {{{
  if (diagnostic) {
    if (length(which(is.na(image)))!=0) {        message(paste("Input Parameter 'image' contains NA elements")) }
    if (length(which(is.na(sfa_models)))!=0) {   message(paste("Input Parameter 'sfa_models' contains NA elements")) }
    if (length(which(is.na(ap.lims.data.stamp)))!=0) {   message(paste("Input Parameter 'ap.lims.data.stamp' contains NA elements")) }
    if (length(which(is.na(fluxes)))!=0) {       message(paste("Input Parameter 'fluxes' contains NA elements")) }
    if (length(which(is.na(outputfile)))!=0) {   message(paste("Input Parameter 'outputfile' contains NA elements")) }
    if (length(which(is.na(outputheader)))!=0) { message(paste("Input Parameter 'outputheader' contains NA elements")) }
    if (length(which(is.na(beamarea)))!=0) {     message(paste("Input Parameter 'beamarea' contains NA elements")) }
    if (length(which(is.na(inside.mask)))!=0) {   message(paste("Input Parameter 'inside.mask' contains NA elements")) }
    message(paste("At initialisation",round(length(which(is.na(resid_im)))/length(resid_im)*100,digits=2),"% of the resid_im matrix are NA"))
  }
  #}}}

  #Get MPI options {{{
  chunk.size=ceiling(npos/getDoParWorkers())
  mpi.opts<-list(chunkSize=chunk.size)
  message("Number of objects per thread:",chunk.size)
  #}}}

  #Perform Subtraction {{{
  modelsource<-foreach (i=1:npos, flux=fluxes, inside=inside.mask, model=sfa_models, .inorder=TRUE, .options.mpi=mpi.opts) %dopar% {
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
    if (!(is.finite(ap.lims.data.stamp))) { sink(type="message") ; stop("Bad limits in Source Subtraction") }
  }
  #}}}
  #Subtract model from image, and add to flux map {{{
  for (i in 1:npos) {
    if (!is.null(modelsource[[i]])) {
      resid_im[ap.lims.data.stamp[i,1]:ap.lims.data.stamp[i,2],ap.lims.data.stamp[i,3]:ap.lims.data.stamp[i,4]]<-
      resid_im[ap.lims.data.stamp[i,1]:ap.lims.data.stamp[i,2],ap.lims.data.stamp[i,3]:ap.lims.data.stamp[i,4]]-modelsource[[i]]
      flux.map[ap.lims.data.stamp[i,1]:ap.lims.data.stamp[i,2],ap.lims.data.stamp[i,3]:ap.lims.data.stamp[i,4]]<-
      flux.map[ap.lims.data.stamp[i,1]:ap.lims.data.stamp[i,2],ap.lims.data.stamp[i,3]:ap.lims.data.stamp[i,4]]+modelsource[[i]]
    }
  }
  #}}}
  #Convert image units back into Jy/bm {{{
  resid_im<-resid_im*beamarea
  #}}}
  if (diagnostic) { message(paste("After assignment",round(length(which(is.na(resid_im)))/length(resid_im)*100,digits=2),"% of the resid_im matrix are NA")) }
  #Output Residual Image & Flux Map {{{
  write.fits.image.file(outputfile,resid_im,outputheader, nochange=TRUE,verbose=verbose,diagnostic=diagnostic)
  write.fits.image.file(fluxmap.outputfile,flux.map,outputheader, nochange=TRUE,verbose=verbose,diagnostic=diagnostic)
  #}}}
  message('=========END==========Source_Subtraction===========END=================\n')
  return=NULL
}
