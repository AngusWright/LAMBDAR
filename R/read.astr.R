read.astrometry <-
function(fitsname, hdu=1){
#Proceedure Reads the header of <fitsname> and
#Returns the astrometry structure

  #Initialise List {{{
  astr.struc<-NULL
  #}}}
  if (hdu==0) { 
    warning("header hdu in read.astrometry should be specificly >= 1") 
    hdu<-1
  }
  oldwarn<-options(warn=-1)
  #Check for Typical Keywords of known Lengths {{{
  astr.struc<-c(astr.struc, list(BITPIX  =read.fitskey("BITPIX",fitsname,hdu=hdu)))
  astr.struc<-c(astr.struc, list(CTYPE   =(c(read.fitskey(c("CTYPE1","CTYPE2","CTYPE3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(NAXIS   =as.numeric(c(read.fitskey(c("NAXIS1","NAXIS2","NAXIS3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(PC      =as.numeric(c(read.fitskey(c("PC1_1","PC1_2","PC1_3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(CDELT   =as.numeric(c(read.fitskey(c("CDELT1","CDELT2","CDELT3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(CRPIX   =as.numeric(c(read.fitskey(c("CRPIX1","CRPIX2","CRPIX3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(CRVAL   =as.numeric(c(read.fitskey(c("CRVAL1","CRVAL2","CRVAL3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(CD=rbind(as.numeric(c(read.fitskey(c("CD1_1","CD1_2","CD_1_3"),fitsname,hdu=hdu))),
                                          as.numeric(c(read.fitskey(c("CD2_1","CD2_2","CD_2_3"),fitsname,hdu=hdu))),
                                          as.numeric(c(read.fitskey(c("CD3_1","CD3_2","CD_3_3"),fitsname,hdu=hdu))))))
  astr.struc<-c(astr.struc, list(CROTA   =as.numeric(c(read.fitskey(c("CROTA1","CROTA2","CROTA3"),fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(LONGPOLE=as.numeric(c(read.fitskey("LONGPOLE",fitsname,hdu=hdu)))))
  astr.struc<-c(astr.struc, list(LATPOLE=as.numeric(c(read.fitskey("LATPOLE",fitsname,hdu=hdu)))))
  #}}}
  # Check for PV2 Keywords (can be up to 21) {{{
  if ((is.finite(read.fitskey("PV2_1",fitsname,hdu=hdu)))) {
    PV2<-foreach (i=1:21, .combine="cbind") %do% {
      as.numeric(read.fitskey(paste("PV2_",i,sep=""),fitsname,hdu=hdu))
    }
    PV2<-PV2[which(is.finite(PV2))]
    astr.struc<-c(astr.struc, list(PV2=PV2))
  }#}}}
  #Convert to standard Format -> use CD matrix {{{
  if (!is.finite(astr.struc$CD[1,1])) {
  astr.struc$CD<-rbind(as.numeric(c(astr.struc$CDELT[1],astr.struc$CROTA[1],astr.struc$CROTA[1])),
                       as.numeric(c(astr.struc$CROTA[2],astr.struc$CDELT[2],astr.struc$CROTA[2])),
                       as.numeric(c(astr.struc$CROTA[3],astr.struc$CROTA[3],astr.struc$CDELT[3])))
  }#}}}
  options(oldwarn)
  #Return {{{
  return(astr.struc)
  #}}}
}
