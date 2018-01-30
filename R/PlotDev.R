#
# Simple wrapper to various plot devices  
#

PlotDev<-function(file='default.X11',height=6,width=6,units='in',res=120,...) {
   
  ulist<-c('in','px','cm','mm')
  if (!units%in%ulist) { stop(paste0("Unit '",units,"' not recognised (must be one of in/px/cm/mm)")) } 
  conv<-c(1,1/res,1/2.54,1/25.4)
  if (grepl('x11',tolower(file),fixed=TRUE)) { 
    #Open X11
    if (units!='in') {
      height=conv[which(ulist==units)]*height
      width=conv[which(ulist==units)]*width
    }
    opened=X11(height=height,width=width,...)
  } else if (grepl('.pdf',tolower(file),fixed=TRUE)) { 
    #Open PDF
    if (units!='in') {
      height=conv[which(ulist==units)]*height
      width=conv[which(ulist==units)]*width
    }
    opened=pdf(file=file,height=height,width=width,...)
  } else if (grepl('.png',tolower(file),fixed=TRUE)) { 
    #Open PNG
    opened=png(filename=file,height=height,width=width,units=units,res=res,...)
  }
  
  return=opened
}
