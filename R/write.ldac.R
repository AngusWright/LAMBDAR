#
# Script to read LDAC catalogues in R 
#

write.ldac<-function(file,object,col.comments,col.units,asctoldac=system('which asctoldac',intern=T),table='OBJECTS',
                    options='',asciifile=paste0(file,'.asc'),configfile=paste0(file,'.conf'),
                    force=TRUE,clean=TRUE,diagnostic=FALSE,...) { 
  require(data.table)
  if (asctoldac == "") { 
    stop("There is no asctoldac binary!") 
  }
  if (missing(col.comments)){ 
    col.comments<-colnames(object)
  } else if (length(col.comments)!=length(colnames(object))) { 
    stop("col.comments must be of the same length as the number of columns in 'object'")
  }
  if (missing(col.units)){ 
    col.units<-rep('',length(colnames(object)))
  } else if (length(col.units)!=length(colnames(object))) { 
    stop("col.comments must be of the same length as the number of columns in 'object'")
  }
  #Make the ASCII file
  if (!file.exists(asciifile) | force) { 
    if (diagnostic) { cat("Writting the object to ascii:\n") }
    write.table(file=asciifile, object,col.names=FALSE,row.names=FALSE,...)
  }
  #Create the configuration file 
  for (var in colnames(object)) { 
    coltype=class(object[,var])
    if (coltype=='numeric') {
      coltype<-"FLOAT"
      htype<-"FLOAT"
      col.depth<-1
    } else if (coltype=='integer') {
      coltype<-"INTEGER"
      htype<-"LONG"
      col.depth<-1
    } else if (coltype=='character') {
      coltype<-"STRING"
      htype<-"STRING"
      object[,var]<-strfix(object[,var])
      col.depth<-nchar(object[1,var])
    } else if (coltype=='factor') {
      coltype<-"DOUBLE"
      htype<-"FLOAT"
      col.depth<-1
    } else { 
      cat("WARNING: unknown column type ",coltype,". Assuming is of type FLOAT!\n")
      coltype<-"FLOAT"
      htype<-"FLOAT"
      col.depth<-1
    }
    index<-which(colnames(object)==var)
    #output the configuration for this parameter
    sink(type='output',append='TRUE',file=configfile)
    cat(paste0("#\n"))
    cat(paste0("COL_NAME = ",colnames(object)[index],"\n"))
    cat(paste0("COL_TTYPE = ",coltype,"\n"))
    cat(paste0("COL_HTYPE = ",htype,"\n"))
    cat(paste0("COL_COMM = \"",col.comments[index],"\"\n"))
    cat(paste0("COL_UNIT = \"",col.units[index],"\"\n"))
    cat(paste0("COL_DEPTH = ",col.depth,"\n"))
    sink(type='output')
  } 
  #Convert the ASCII to LDAC
  system(paste(asctoldac,'-i',asciifile,'-t',table,'-o',file,'-c',configfile,ifelse(diagnostic,'','2> /dev/null')))


  #Clean
  if (clean) { 
    if (diagnostic) { cat("Removing temporary files\n") }
    system(paste("rm",asciifile,configfile))
  }
  if (diagnostic) { cat("Returning!\n") }
  #return
  return=cat
}


