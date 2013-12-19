make_deblended_weightmap <-
function(env=NULL,outenv=NULL) {
  # Procedure produces an array of stamps with apertures that
  # have been weighted according to their amount of blending.
  # This provides a launchpad for the final determination of
  # per-object fluxes.

  # Load Parameter Space
  if(is.null(env)) {
    stop("No Parameter Space Environment Specified in function call")
  }
  if(is.null(outenv)) { outenv<-env }
  attach(env, warn.conflicts=FALSE)
  #on.exit(detach('env'))


  if (!quiet) { cat("Make_Deblended_Weightmap    ") }
  message('--------------------Make_Deblended_Weightmap---------------------------')
  #Setup Sizes
  ngal<-length(id_g)

  ##Produce Deblended Stamps
  #For each stamp, determine the limits of the stamp
  #in full (weighted) mask array-space, and divide stamp
  #by the full value.
  dbw_map<-foreach(i=1:ngal,xlo=stamp_lims[,1],xup=stamp_lims[,2],ylo=stamp_lims[,3],yup=stamp_lims[,4],.inorder=TRUE)%dopar%{
      ##Perform Calculation and assignment
      #If the sums of the arrays are both nonzero
      if ((sum(wsfa[[i]]) > 0)&&(sum(image.env$wfa[xlo:xup,ylo:yup])>0)) {
        #Check for errors in array sizes
        if (dim(wsfa[[i]])[1] != dim(image.env$wfa[xlo:xup,ylo:yup])[1]) {
          print(paste(dim(wsfa[[i]]),"!=",dim(image.env$wfa[xlo:xup,ylo:yup])))
          print(paste("[",xlo,":",xup,",",ylo,":",yup,"]"))
          sink(type="message")
          stop("Dimensions in matrix division are not equal")
        }
        #Calculate Deblended Value and assign
        fmap<-wsfa[[i]]/image.env$wfa[xlo:xup,ylo:yup]
        #Check for production of NA/NaN/Infs
        if (length(which((is.na(fmap)))) != 0 ) {
          #warning("NA values produced in map. They will be zero-ed.")
          fmap[which((is.na(fmap)))]<-1.0
        }
      } else {
        #If sums are equal to zero, short-cut by assigning deblended map to zero
        if ((sum(wsfa[[i]])==0)|(sum(image.env$wfa[xlo:xup,ylo:yup])==0)) { fmap<-array(0, dim=dim(wsfa[[i]])) }
        #If sums are -ve, there may be an error: The arrays should be +ve ± noise, or 0.
        #Produces warning, and set to zero.
        if ((sum(wsfa[[i]])<0)|(sum(image.env$wfa[xlo:xup,ylo:yup])<0)) {
          warning("The sum value of the wsfa and wfa matricies are negative
          - check for an offest error in FITS image or calculations")
          fmap<-array(0, dim=dim(wsfa[[i]]))
        }
      }
      fmap
  }
  #Check for production of NA/NaN/Infs
  #for (i in 1:ngal) {
  #  if (length(which((is.na(dbw_map[[i]])))) != 0 ) {
  #    warning(paste("NA values produced in dbw map",i,". They will be zero-ed."))
  #    message(paste("NA values produced in dbw map",i,". They will be zero-ed."))
  #    dbw_map[[i]][which((is.na(dbw_map[[i]])))]<-0.0
  #  }
  #}
  #-----Diagnostic-----#
  #if (diagnostic) { message(paste("After assignment",round(length(which(is.na(dbw_map)))/length(dbw_map)*100,digits=2),"% of the dbw_map matrix are NA"))
  #                  #message(paste("After assignment",round(length(which(dbw_map==0))/length(dbw_map)*100,digits=2),"% of the dbw_map matrix are zero"))
  #}
  message('=========END=======Make_Deblended_Weightmap========END=================\n')

  detach(env)
  #Return Deblended Stamps
  dbw_map
}
