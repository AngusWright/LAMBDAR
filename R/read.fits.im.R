#
# Wrapper to readFITS for reading a fits image (and outputting it like read.fits)
#

read.fits.im<-function(...,comments=TRUE,strip='') { 
  dat<-readFITS(...) 
  im<-list(hdr=list(NULL),dat=list(dat$imDat))
  im$hdr[[1]]<-.parsehdr(paste0(dat$header,collapse=''),comments=comments,strip=strip)
  return=im
}

.parsehdr<-function(rawhead,comments=TRUE,strip=""){
  key = {}
  value = {}
  comment = {}
  for (i in 1:ceiling(nchar(rawhead)/80)) {
    row = substr(rawhead, ((80 * (i - 1)) + 1), (80 * i))
    k = paste(strsplit(substr(row, 1, 8), " +")[[1]], collapse = "", 
              sep = "")
    if (k == "HIERARCH") {
      charend = (min(which(strsplit(row, "")[[1]] == "=")) - 
                 2)
      k = strip(paste(strsplit(row, "")[[1]][10:charend], 
                      collapse = "", sep = ""), strip = strip)
      leftover = substr(row, charend + 3, 80)
    }
    else {
      leftover = substr(row, 10, 80)
    }
    key = c(key, k)
    if (k == "COMMENT" | k == "HISTORY") {
      v = ""
    }
    else {
      poss = strsplit(leftover, "/")[[1]]
      while (((length(grep("'", strsplit(poss[1], "")[[1]]))%%2) != 
              0) & (length(poss) > 2)) {
        poss[1] = paste(poss[1], "/", poss[2], collapse = "", 
                        sep = "")
        poss = poss[-2]
      }
      temp = strip(poss[1], strip = strip)
      v = paste(temp, collapse = " ")
    }
    value = c(value, v)
    if (comments) {
      if (k == "COMMENT" | k == "HISTORY") {
        m = strip(substr(row, 9, 80), strip = " ")
      }
      else {
        m = strip(paste(poss[2:length(poss)], collapse = "/"), 
                  strip = " ")
      }
      comment = c(comment, m)
    }
  }
  if (comments) {
    hdr = cbind(key, value, comment)
  }
  else {
    hdr = cbind(key, value)
  }
  return(hdr)
}
