#
# Wrapper to readFITS for reading a fits catalogue
#

read.fits.cat<-function(filename,data.table=FALSE,stringsAsFactors=FALSE,quiet=TRUE) { 
    
  zz <- file(description = filename, open = "rb")
  output<-capture.output(header0 <- readFITSheader(zz,fixHdr='substitute')) # read primary header
  output<-capture.output(header <- readFITSheader(zz,fixHdr='substitute')) # read extension header
  cat <- readFITSbintable(zz, header)
  close(zz)
  if (data.table) { 
    out<-as.data.table(cat$col,stringsAsFactors=stringsAsFactors)
  } else { 
    out<-as.data.frame(cat$col,stringsAsFactors=stringsAsFactors)
  }
  colnames(out)<-cat$colNames
  return=out
}


