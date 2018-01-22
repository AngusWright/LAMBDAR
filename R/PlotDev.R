#
# Simple wrapper to CairoPNG and png that 
# determines which to use:
#   Preferentially CairoPNG if available.
#

PlotDev<-function(file='default.X11',height=6,width=6,units='in',res=120,...) {
   
  #If Cairo is available
  if (requireNamespace('Cairo',quietly=TRUE)) {
    if (grepl('X11',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoX11
      opened=Cairo::Cairo(type='X11',height=height,width=width,units=units,...)
    } else if (grepl('.pdf',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoPDF
      opened=Cairo::Cairo(type='PDF',file=file,height=height,width=width,units=units,...)
    } else if (grepl('.png',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoPNG
      opened=Cairo::Cairo(type='PNG',filename=file,height=height,width=width,units=units,res=res,...)
    }
  } else {
    if (grepl('X11',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoX11
      opened=X11(height=height,width=width,units=units,...)
    } else if (grepl('.pdf',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoPDF
      opened=pdf(file=file,height=height,width=width,...)
    } else if (grepl('.png',file,fixed=TRUE,ignore.case=TRUE)) { 
      #Open CairoPNG
      opened=png(filename=file,height=height,width=width,units=units,res=res,...)
    }
  #}
  return=opened
}
