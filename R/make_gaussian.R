make_gauss <-
function(fwhm,asperpix,xlen,ylen,nsig=6, height=1) {

    sigma_as=fwhm/(2.*sqrt(2.*log(2.)))
    sigma_p=sigma_as/asperpix
    x0=floor(nsig*sigma_p)
    y0=floor(nsig*sigma_p)

    psflim<-floor(2*nsig*sigma_p)
    
    zp<-foreach (xi=1:psflim, .combine='cbind') %:%
      foreach (yi=1:psflim, .combine='c') %do% {
        height*exp(-(((xi-x0)^2/(2*sigma_p^2))+((yi-y0)^2/(2*sigma_p^2))))
      }

    psf<-array(0, dim=c(xlen,ylen))
    psf[1:psflim,1:psflim]<-zp
    return(psf)
}
