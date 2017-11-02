get.cog<-
function(zdist, centre=NULL,sample=NULL,proj=NULL,SNR=FALSE,poly.degree=4,n.per.bin=5,flexible=TRUE,na.rm=TRUE){
#Details {{{
#function returns the FHWM of the image from
#*maxima*, in pixels.
#}}}

  #Setup Radius-map {{{
  if (is.null(centre)) { centre<-as.numeric(which(zdist==max(zdist), arr.ind=TRUE)) }
  moving<-TRUE
  while (moving) { 
    x = seq(1,dim(zdist)[1])
    y = seq(1,dim(zdist)[2])
    xy = expand.grid(x,y)
    r=sqrt((xy[,1]-centre[1]+0.5)^2+(xy[,2]-centre[2]+0.5)^2)
    #im.rad.x<-min(centre[1],length(zdist[,1])-centre[1])-1
    #im.rad.y<-min(centre[2],length(zdist[1,])-centre[2])-1
    #lim<-c(centre[1]-im.rad.x, centre[1]+im.rad.x, centre[2]-im.rad.y, centre[2]+im.rad.y)
    #x = seq(floor(-im.rad.x), floor(im.rad.x), length=im.rad.x*2+1)
    #y = seq(floor(-im.rad.y), floor(im.rad.y), length=im.rad.y*2+1)
    #xy = expand.grid(x,y)
    #r=sqrt(xy[,1]^2+xy[,2]^2)
    if (!is.null(proj)) {
      if (length(proj==2) & is.finite(proj[1])) {
        #Proj is present, correctly used, and the source is resolved (Axrat != 0/0)
        #if (length(proj)!=2) { stop("Projection parameters must be length 2: c(Axrat,PA)") }
        ang = atan2(xy[,1], xy[,2]) - proj[2] * pi/180
        r=sqrt((r*sin(ang)/proj[1])^2+(r*cos(ang))^2)
      }
    }
    if (flexible) { 
      val<-zdist[cbind(xy[,1],xy[,2])]*(1-r/max(r))^4
      new.cen<-as.numeric(xy[which(val==max(val,na.rm=TRUE))[1],])
      if (all(new.cen==centre)) { 
        moving<-FALSE
      } else { 
        centre<-new.cen
      }
    } else { 
      moving<-FALSE
    }
  }
  #}}}
  #Calculate COG {{{
  tmp.order<-order(r)
  #zvec<-as.numeric(zdist[lim[1]:lim[2],lim[3]:lim[4]])
  zvec<-as.numeric(zdist)
  zbin<-zvec[tmp.order]
  r<-r[tmp.order]
  if (na.rm) { 
    good<-which(is.finite(zbin))
  } else {
    good<-1:length(zbin)
  }
  cog<-data.frame(x=r[good],y=cumsum(zbin[good]))
  tab=factor(cog$x)
  lev<-levels(tab)
  #Sample if wanted
  if (!is.null(sample) && length(lev) > sample) { 
    #If there are more radius bins than wanted samples
    #Randomly sample n bins without duplication
    lev<-sample(lev,sample)
    ind<-which(cog$x %in% lev) 
    cog<-cog[ind,]
    tab=factor(cog$x)
    lev<-levels(tab)
  } 
  if (!is.null(sample) && (length(lev) <= sample & length(cog$x)>n.per.bin*sample)) { 
    #If there are fewer radius bins than wanted samples
    #Randomly thin the whole sample so that there are at most n.per.bin samples per bin
    shuffle<-order(runif(length(cog$x))) 
    tx<-cog$x[shuffle]
    for (i in 1:n.per.bin) { 
      #select the unduplicated elements
      ind<-which(!duplicated(tx))
      #Remove them 
      tx<-tx[-ind]
      shuffle<-shuffle[-ind]
    }
    #shuffle now contains only things that were (randomly) 
    #duplicated in each of the n.per.bin loops
    cog<-cog[-shuffle,]
    tab=factor(cog$x)
    lev<-levels(tab)
  }
  if (any(as.numeric(tab)>1)) { 
    avg.x<-as.numeric(lev)
    avg.y<-rep(NA,length(avg.x))
    tx<-cog$x
    ty<-cog$y
    for (val in 1:length(avg.x)) { 
      ind=which(tx==lev[val])
      avg.y[val]<-mean(ty[ind])
      tx=tx[-ind]
      ty=ty[-ind]
    }
    avg<-data.frame(x=avg.x,y=avg.y)
  } else { 
    avg<-cog
  }
  #}}}

  #Do various useful fits to the CoG {{{ 
  if (length(which(is.finite(avg$y)))>poly.degree) { 
    r<-count<-0
    while (r < 0.999 & count < 10) { 
      fit<-lm(avg$y ~ poly(avg$x, degree=poly.degree+count, raw=TRUE))
      r<-summary(fit)$r.squared
      count<-count+1
    }
    
    avg$fitted <- fitted(fit) 
    #Calculate the derivative
    deriv_coef<-function(x) {
      x <- coef(x)
      stopifnot(names(x)[1]=="(Intercept)")
      y <- x[-1]
      stopifnot(all(grepl("^poly", names(y))))
      px <- as.numeric(gsub("poly\\(.*\\)","",names(y)))
      rr <- setNames(c(y * px, 0), names(x))
      rr[is.na(rr)] <- 0
      rr
    }
    avg$slope <- model.matrix(fit) %*% matrix(deriv_coef(fit), ncol=1)
    #Calculate the second derivative
    fit<-lm(avg$slope ~ poly(avg$x, degree=poly.degree, raw=TRUE))
    avg$concav <- model.matrix(fit) %*% matrix(deriv_coef(fit), ncol=1)
  } else { 
    #Not enough data pts vs degrees of freedom 
    avg$fitted<-rep(NA,length(avg$y))
    avg$slope<-rep(NA,length(avg$y))
    avg$concav<-rep(NA,length(avg$y))
  }
  #}}}

  if (SNR) {
    cog$y<-cog$y/sqrt((pi*cog$x^2))
    avg$y<-avg$y/sqrt((pi*avg$x^2))
  }
  if (!is.null(sample) && sample > length(cog$y)) {
    sample<-seq.int(1,length(cog$y),length.out=sample)
    cog<-cog[sample,]
  }
  if (!is.null(sample) && sample > length(avg$y)) {
    sample<-seq.int(1,length(avg$y),length.out=sample)
    avg<-avg[sample,]
  }
  #}}}
  #Return FWHM {{{
  return=list(all=cog,avg=avg,cen=centre)
  #}}}
}

