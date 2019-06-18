

strfix<-function(vec,len=max(nchar(vec))){
  nch<-nchar(vec)
  max<-max(nch)
  if (any(nch!=max)) { 
    vec[nch!=max]<-paste0(strtrim(rep(paste0(rep(" ",max),collapse=''),length(which(nch!=max))),(max-nch)[nch!=max]),vec[nch!=max])
  }
  return=vec
}
