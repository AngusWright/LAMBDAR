checkgrid <-
function(x,y,xstep,ystep,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T){
if(any(is.na(x))|any(is.na(y))){stop(paste("Checkgrid returns NA at input.\nis.na(x)=",length(which(is.na(x))),"of",length(x),"\nis.na(y)=",length(which(is.na(y))),"of",length(y),"\nInputs are:\nxcen=",xcen,"\nycen=",ycen,"\naxrat=",axrat,"\naxang=",axang,"\nmajax=",majax,"\ndeg=",deg))}
MxMy=inap(x,y,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
full= LxLy & LxHy & HxLy & HxHy & MxMy
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
part= ((LxLy | LxHy | HxLy | HxHy) & full==F)
empty= (! full) & (! part)

#if(all(part==FALSE)){(part=xcen<x+xstep/2 & xcen>x-xstep/2 & ycen<y+ystep/2 & ycen>y-ystep/2)} #this is to catch situations where the aperture is smaller than the pixel grid
Npart= LxLy + LxHy + HxLy + HxHy
return=list(logic=cbind(part=part,empty=empty,full=full,mid=MxMy),Npart=Npart)
}

