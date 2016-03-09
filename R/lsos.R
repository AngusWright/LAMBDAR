lsos <-
function (pattern="", envir=NULL, pos = 1, order.by="Size", decreasing=TRUE, head=TRUE, n=10,noShow="",exact=FALSE) {
  #Get variable names {{{
  if (is.null(envir)) {
    #Get variable names using 'pos' {{{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos,inherits=TRUE)))
    names <- ls(pos = pos, pattern = pattern)
    #}}}
  } else {
    #Get variable names using 'envir' {{{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, envir = envir,inherits=TRUE)))
    names <- ls(envir = envir, pattern = pattern)
    #}}}
  }#}}}
  if (exact) {
    names<-names[which(names == pattern)]
  }
  noShow<-tolower(noShow)
  if (length(names)!=0) {
    #There are objects, get details {{{
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) { capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    if (!('vals'%in%noShow)) {
      obj.vals <-
      (napply(names, function(x) {
        if (length(x)==0) {
          #object is NULL
          capture.output(cat("NULL"))
        } else if (length(class(x))>1) {
          if (class(x)[1]=="file") {
            #object is a file connection
            capture.output(cat("<connection>"))
          } else {
            #object is unknown
            capture.output(cat("UNKNOWN TYPE"))
          }
        } else if ((class(x)=="numeric")) {
          #object is numeric
          if (length(x)>5) {
            capture.output(cat(paste(paste(round(x[1:5],digits=4),collapse=" "),"...")))
          } else {
            capture.output(cat(paste(round(x[1:length(x)],digits=4), collapse=" ")))
          }
        } else if (class(x)=="data.frame") {
            #object is a data frame
            capture.output(cat("Data-Frame of length ",length(x)))
        } else if (class(x)=="density") {
            #object is a list
            capture.output(cat("Density Object"))
        } else if (class(x)=="list") {
            #object is a list
            capture.output(cat("List of length ",length(x)))
        } else if (class(x)=="function") {
            #object is a function
            capture.output(cat("<function>"))
        } else if (class(x)=="environment") {
            #object is an environment
            capture.output(print(x))
        } else if (class(x)=="character") {
          #object is a character string
          capture.output(cat(gsub("\n","",paste(strtrim(paste(x,collapse=" "),35),"..."))))
        } else if ((class(x)=="matrix")||class(x)=="array") {
          if (class(x[1])=="numeric") {
            #object is a numeric matrix/array
            if (length(x)>5) {
              capture.output(cat(paste(paste(round(as.numeric(x[1:5]),digits=4),collapse=" "),"...")))
            } else {
              capture.output(cat(paste(round(as.numeric(x[1:length(x)]),digits=4),collapse=" ")))
            }
          } else {
            #object is a character matrix/array
            if (length(x)>5) {
              capture.output(cat(paste(paste(x[1:5],collapse=" "),"...")))
            } else {
              capture.output(cat(paste(x[1:length(x)],collapse=" ")))
            }
          }
        } else {
          #object is logical or integer
          if (length(x)>5) {
            capture.output(cat(paste(paste(x[1:5],collapse=" "),"...")))
          } else {
            capture.output(cat(paste(x[1:length(x)],collapse=" ")))
          }
        }
      } ))
    } else {
      obj.vals <- rep(NA,length(names))
    }
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim, obj.vals)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns", "Value")
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
    out <- out[,which(!c("type", "size", "pretty", "dim", "dim", "vals") %in% noShow)]
    if (head) { out <- head(out, n) }
    #}}}
  } else {
    #There are no matching object, return null {{{
    out<-as.character(NULL)
    #}}}
  }
  #return {{{
  return(out)
  #}}}
}
