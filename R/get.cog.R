get.cog<-
function(zdist, centre=NULL,sample=NULL,proj=NULL){
#Details {{{
#function returns the FHWM of the image from
#*maxima*, in pixels.
#}}}

  #Setup Radius-map {{{
  if (is.null(centre)) { centre<-as.numeric(which(zdist==max(zdist), arr.ind=TRUE)) }
  im.rad.x<-min(centre[1],length(zdist[,1])-centre[1])-1
  im.rad.y<-min(centre[2],length(zdist[1,])-centre[2])-1
  lim<-c(centre[1]-im.rad.x, centre[1]+im.rad.x, centre[2]-im.rad.y, centre[2]+im.rad.y)
  x = seq(floor(-im.rad.x), floor(im.rad.x), length=im.rad.x*2+1)
  y = seq(floor(-im.rad.y), floor(im.rad.y), length=im.rad.y*2+1)
  xy = expand.grid(x,y)
  r=sqrt(xy[,1]^2+xy[,2]^2)
  if (!is.null(proj)) {
    if (length(proj)!=2) { stop("Projection parameters must be length 2: c(Axrat,PA)") }
    ang = atan2(xy[,1], xy[,2]) - proj[2] * pi/180
    xy[,1]=r*sin(ang)
    xy[,2]=r*cos(ang)
    xy[,1]=xy[,1]/proj[1]
    r=sqrt(xy[,1]^2+xy[,2]^2)
  }
  #}}}
  #Calculate COG {{{
  tmp.order<-order(r)
  zvec<-as.numeric(zdist[lim[1]:lim[2],lim[3]:lim[4]])
  zbin<-zvec[tmp.order]
  cog<-data.frame(x=r[tmp.order],y=cumsum(zbin))
  if (!is.null(sample) && sample > length(cog$y)) {
    sample<-seq.int(1,length(cog$y),length.out=sample)
    cog<-cog[sample,]
  }
  #}}}
  #Return FWHM {{{
  return=cog
  #}}}
}

