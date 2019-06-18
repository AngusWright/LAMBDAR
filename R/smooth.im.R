
smooth.im<-function(im,filter.sd.pix,normalise=TRUE) { 
  require(LAMBDAR)
  if (length(filter.sd.pix)==1) { 
    filter.sd.pix<-rep(filter.sd.pix,2)
  }
  psf.x = matrix(1:ncol(im), nrow = nrow(im), ncol = ncol(im),byrow=T)
  psf.y = matrix(1:nrow(im), nrow = nrow(im), ncol = ncol(im))
  psf <- exp(-(((psf.x - ncol(im)/2)^2/(2 * filter.sd.pix[1]^2)) + ((psf.y - 
      nrow(im)/2)^2/(2 * filter.sd.pix[2]^2))))
  conv<-convolve.psf(psf,im)
  if (normalise) { 
    conv<-conv/sum(psf)
  }
  return=conv
}
