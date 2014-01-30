make_deblended_weightmap <-
function(env=NULL,outenv=NULL) {
#Procedure produces an array of stamps with apertures that
#have been weighted according to their amount of blending.
#This provides a launchpad for the final determination of
#per-object fluxes.

  if (!quiet) { cat("Make_Deblended_Weightmap    ") }
  message('--------------------Make_Deblended_Weightmap---------------------------')

  # Load Parameter Space {{{
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #}}}

  #Setup Sizes {{{
  ngal<-length(id_g)
  #}}}

  #Produce Deblended Stamps {{{
  #Details {{{
  #For each stamp, determine the limits of the stamp
  #in full (weighted) mask array-space, and divide stamp
  #by the full value. }}}
  #Perform Calculation and assignment {{{
  dbw_map<-foreach(i=1:ngal,xlo=stamp_lims[,1],xup=stamp_lims[,2],ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE)%dopar%{
      #Check for errors and create deblended stamps {{{
      if ((sum(wsfa[[i]]) > 0)&&(sum(image.env$wfa[xlo:xup,ylo:yup])>0)) {
        #Aperture stamps are error free {{{
        #Check for errors in array sizes {{{
        if (dim(wsfa[[i]])[1] != dim(image.env$wfa[xlo:xup,ylo:yup])[1]) {
          print(paste(dim(wsfa[[i]]),"!=",dim(image.env$wfa[xlo:xup,ylo:yup])))
          print(paste("[",xlo,":",xup,",",ylo,":",yup,"]"))
          sink(type="message")
          stop("Dimensions in matrix division are not equal")
        }#}}}
        #Calculate Deblended Value and assign {{{
        fmap<-wsfa[[i]]/image.env$wfa[xlo:xup,ylo:yup]
        #}}}
        #Check for production of NA/NaN/Infs {{{
        if (length(which((is.na(fmap)))) != 0 ) {
          #warning("NA values produced in map. They will be zero-ed.")
          fmap[which((is.na(fmap)))]<-1.0
        }#}}}
        #}}}
      } else {
      #Aperture stamps may have errors {{{
        #If sums are equal to zero, short-cut by assigning deblended map to zero {{{
        if ((sum(wsfa[[i]])==0)|(sum(image.env$wfa[xlo:xup,ylo:yup])==0)) { fmap<-array(0, dim=dim(wsfa[[i]])) }
        #}}}
        #If sums are -ve, there may be an error {{{
        #Details {{{
        #The arrays should be +ve ± noise, or 0. Produces warning, and set to zero. }}}
        if ((sum(wsfa[[i]])<0)|(sum(image.env$wfa[xlo:xup,ylo:yup])<0)) {
          warning("The sum value of the wsfa and wfa matricies are negative
          - check for an offest error in FITS image or calculations")
          fmap<-array(0, dim=dim(wsfa[[i]]))
        }#}}}
        #}}}
      }
      fmap
      #}}}
  }#}}}
  #}}}

  message('=========END=======Make_Deblended_Weightmap========END=================\n')

  #Return {{{
  detach(env)
  #Return Deblended Stamps
  return=dbw_map
  #}}}
}
