#
# Script to read LDAC catalogues in R 
#

read.ldac<-function(file,ldactoasc=system('which ldactoasc',intern=T),table='OBJECTS',
                    options='',outfile=paste0(file,'_ldac.csv'),
                    headcat=paste0(file,'_head.csv'),bodycat=paste0(file,'_body.csv'),
                    force=TRUE,clean=TRUE,diagnostic=FALSE,...) { 
  require(data.table)
  if (ldactoasc == "") { 
    stop("There is no ldactoasc binary!") 
  }
  #Make the ASCII file
  if (!file.exists(outfile) | force) { 
    if (diagnostic) { cat("Converting ldac to ascii:\n\n") }
    if (diagnostic) { cat(paste(ldactoasc,'-i',file,' -t',table,options,'>',outfile,ifelse(diagnostic,"",'2> /dev/null'),'\n\n')) }
    system(paste(ldactoasc,'-i',file,' -t',table,options,'>',outfile,ifelse(diagnostic,"",'2> /dev/null')))
  }
  if (!file.exists(headcat) | force) {
    if (diagnostic) { cat("Reading the header:\n") }
    system(paste("grep \"^#\"",outfile,"| awk '{print $2,$3}' > ",headcat))
  }
  ##Read the data
  #If the catalogue contains strings, use csv. Otherwise, use ascii
  if (grepl("-s",options)) { 
    if (!file.exists(bodycat) | force) {
      if (diagnostic) { cat("Getting main body of data:\n") }
      #Convert to CSV 
      system(paste0("grep -v \"^#\" ",outfile," | sed 's/ /,/g' > ",bodycat))
    }
    cat<-try(fread(bodycat,header=FALSE,sep=',',...))
  } else { 
    system(paste0("grep -v \"^#\" ",outfile," > ",bodycat))
    cat<-try(fread(bodycat,header=FALSE,...))
  }
  if (class(cat)=='try-error') { 
    system(paste0("grep -v \"^#\" ",outfile," > ",bodycat))
    cat<-fread(bodycat,header=FALSE,...)
  }
  if (diagnostic) { cat(paste("Data has dimensions: ",paste(collapse=':',dim(cat)),"\n")) }
  ##Read the header 
  head<-fread(headcat,header=FALSE,...)
  if (diagnostic) { cat(paste("Header has length: ",length(head$V2),"\n")) }
  if (length(head$V2)!=length(cat)) {
    if (diagnostic) { cat("Reading the header:\n") }
    fullhead<-rep("",rev(head$V1)[1])
    if (length(fullhead)!=length(cat)) {
      stop("Cannot determine header for catalogue because of a mix between strings and vectors")
    }
    count<-1
    for (i in 1:length(head$V1)) { 
      vec.count<-0
      while (head$V1[i]>count) {
        fullhead[count]<-paste(head$V2[i-1],vec.count,sep='.')
        vec.count<-vec.count+1
        count<-count+1
      }
      fullhead[count]<-head$V2[i]
      count<-count+1
    }
  } else { 
    fullhead<-head$V2
  } 
  ##update the header names
  if (diagnostic) { cat("Assigning header:\n") }
  colnames(cat)<-fullhead
  #Clean
  if (clean) { 
    if (diagnostic) { cat("Removing temporary files\n") }
    system(paste("rm",outfile,headcat,bodycat))
  }
  if (diagnostic) { cat("Returning!\n") }
  #return
  return=cat
}


