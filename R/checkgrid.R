checkgrid <-
function(x,y,xstep,ystep,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=T){
#Check for calling errors {{{
if(any(is.na(x))|any(is.na(y))){stop(paste("Checkgrid returns NA at input.\nis.na(x)=",length(which(is.na(x))),"of",length(x),"\nis.na(y)=",length(which(is.na(y))),"of",length(y),"\nInputs are:\nxcen=",xcen,"\nycen=",ycen,"\naxrat=",axrat,"\naxang=",axang,"\nmajax=",majax,"\ndeg=",deg))}
#}}}
#Determine if pixel corners and centre are covered by aperture {{{
MxMy=inap(x,y,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax,deg=deg)
#}}}
#Check if whole pixel is coveres {{{
full= LxLy & LxHy & HxLy & HxHy & MxMy
#}}}
#The proportion of the pixel that is covered {{{
Npart= LxLy + LxHy + HxLy + HxHy
#}}}
#Check if aperture breaches the pixel at all {{{
LxLy=inap(x-xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
LxHy=inap(x-xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxLy=inap(x+xstep/2,y-ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
HxHy=inap(x+xstep/2,y+ystep/2,xcen=xcen,ycen=ycen,axrat=axrat,axang=axang,majax=majax+max(c(xstep,ystep)*sqrt(2)/2),deg=deg)
part= ((LxLy | LxHy | HxLy | HxHy) & full==F)
#}}}
#Check if the pixel is empty {{{
empty= (! full) & (! part)
#}}}
#Return {{{
return=list(logic=cbind(part=part,empty=empty,full=full,mid=MxMy),Npart=Npart)
#}}}
}
