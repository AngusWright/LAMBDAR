get.confidence<-
function(zdist,confidence=0.95, nsteps=100, value=FALSE){
#function returns the radius from the
#*image centre*, in pixels, that contains
#<confidence> proportion of the image
#if value=TRUE, then the value of the PSF
#at that confidence radius is returned instead
  
  im.len.x<-dim(zdist)[1]
  im.len.y<-dim(zdist)[2]
  x = seq(floor(-im.len.x/2), floor(im.len.x/2), length=im.len.x)
  y = seq(floor(-im.len.y/2), floor(im.len.y/2), length=im.len.y)
  xy = expand.grid(x,y)

  r=sqrt(xy[,1]^2+xy[,2]^2)

  tmp.order<-order(r)
  zvec<-as.numeric(zdist)
  zcumsum<-cumsum(zvec[tmp.order])
  rcut<-max(abs(r[tmp.order][zcumsum<=sum(zdist)*confidence]))
  if (value) {
    return=zvec[which(r==rcut)]
  } else {
    return=ceiling(rcut)
  }
}

