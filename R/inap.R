inap <-
function(x,y,xcen=0,ycen=0,axrat=1,axang=0,majax=1,deg=TRUE){
  if(majax>0){
    #Error catch {{{
    if(any(is.na(x))|any(is.na(y))){stop(paste("Inap returns NA at input.\nis.na(x)=",length(which(is.na(x))),"of",length(x),"\nis.na(y)=",length(which(is.na(y))),"of",length(y),"\nInputs are:\nxcen=",xcen,"\nycen=",ycen,"\naxrat=",axrat,"\naxang=",axang,"\nmajax=",majax,"\ndeg=",deg))}
    #}}}
    x=x-xcen;y=y-ycen;r=sqrt(x^2+y^2)
    #Error catch {{{
    if(any(is.na(x))|any(is.na(y))){stop(paste("Inap returns NA at step 1.\nis.na(x)=",length(which(is.na(x))),"of",length(x),"\nis.na(y)=",length(which(is.na(y))),"of",length(y),"\nInputs are:\nxcen=",xcen,"\nycen=",ycen,"\naxrat=",axrat,"\naxang=",axang,"\nmajax=",majax,"\ndeg=",deg))}
    #}}}
    if(deg){ang=atan2(x,y)-axang*pi/180}else{ang=atan2(x,y)-axang}
    x=r*sin(ang);y=r*cos(ang)
    x=x/axrat
    x=x/majax;y=y/majax
    #Error catch {{{
    if(any(is.na(x))|any(is.na(y))){stop(paste("Inap returns NA at step 4.\nwhich(is.na(x))=",which(is.na(x)),"\nwhich(is.na(y))=",which(is.na(y)),"\nInputs are:\nxcen=",xcen,"\nycen=",ycen,"\naxrat=",axrat,"\naxang=",axang,"\nmajax=",majax,"\ndeg=",deg))}
    #}}}
    #Return {{{
    return=(x^2+y^2<=1)
    #}}}
  }else{
    return=(x==xcen & y==ycen)
  }
}

