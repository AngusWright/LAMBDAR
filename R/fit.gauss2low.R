fit.gauss2low<-function (x, bw = diff(quantile(x,pnorm(c(-2,2))))/1e4, from=median(x)-10*mad(x), to=median(x)+10*mad(x), plot = FALSE, fixat = NULL, warnOnly = TRUE, ...) 
{
    if (length(which(is.finite(x)))==0) { 
      return(list(mu=NA,sd=NA,amp=NA,mode=NA,fit=NA))
    }
    if (any(!is.finite(x))) { 
      x<-x[which(is.finite(x))]
    }

    dens <- density(x, bw = bw, kernel = "rect", from=from, to=to, ...)
    if (!is.null(fixat)) {
        ind <- which(dens$x <= fixat)
        if (plot) {
            magplot(dens)
            lines(dens$x, max(dens$y)/dnorm(fixat,fixat,mad(x))*dnorm(dens$x, fixat, mad(x)), col = "grey", lty = 2)
        }
        fit <- try(nls(y ~ a1 * 1/(sqrt(2 * pi) * s1) * exp(-(x - 
            fixat)^2/(2 * s1^2)), start = list(a1 = max(dens$y)/dnorm(fixat,fixat,mad(x)), 
            s1 = mad(x)), data = list(x = dens$x[ind], y = dens$y[ind], 
            fixat = fixat),control=list(warnOnly=FALSE)))
        if (class(fit)=='try-error') { 
          #try again with more restrictive limits 
          ind <- which(dens$x <= fixat & dens$x >= fixat-mad(x))
          fit <- nls(y ~ a1 * 1/(sqrt(2 * pi) * s1) * exp(-(x - 
              fixat)^2/(2 * s1^2)), start = list(a1 = max(dens$y)/dnorm(fixat,fixat,mad(x)), 
              s1 = mad(x)), data = list(x = dens$x[ind], y = dens$y[ind], 
              fixat = fixat),control=list(warnOnly=warnOnly))
        } 
        if (plot) {
            lines(dens$x, coef(fit)[1] * 
                dnorm(dens$x, fixat, coef(fit)[2]), col = "blue", lwd = 2)
            lines(dens$x, dens$y - coef(fit)[1] *
                dnorm(dens$x, fixat, coef(fit)[2]), col = "red", lwd = 2)
            lines(x = dens$x[ind], y = dens$y[ind], col = "orange",lwd=2)
            abline(h=0,lty=3)
        }
        summ<-try(coef(summary(fit)))
        if (class(summ)=='try-error') { 
          mu<-sd<-amp<-muerr<-sderr<-amperr<-NA
        } else { 
          mu = fixat 
          sd =summ['s1','Estimate']
          amp=summ['a1','Estimate']
          muerr = NA 
          sderr =summ['s1','Std. Error']
          amperr=summ['a1','Std. Error']
        }
    } else {
        mode <- dens$x[which.max(dens$y)]
        ind <- which(dens$x <= mode + mad(x)/2)
        if (plot) {
            magplot(dens)
            lines(dens$x, max(dens$y)/dnorm(mode,mode,mad(x))*dnorm(dens$x, mode, mad(x)), col = "grey", lty = 2)
        }
        fit <- try(nls(y ~ a1 * 1/(sqrt(2 * pi) * s1) * exp(-(x - 
            m1)^2/(2 * s1^2)), start = list(a1 = max(dens$y)/dnorm(mode,mode,mad(x)), 
            m1 = mode, s1 = mad(x)), data = list(x = dens$x[ind], 
            y = dens$y[ind]),control=list(warnOnly=FALSE)))
        if (class(fit)=='try-error') { 
          #try again with more restrictive limits 
          ind <- which(dens$x <= mode + mad(x)/2 & dens$x >= mode -mad(x))
          fit <- try(nls(y ~ a1 * 1/(sqrt(2 * pi) * s1) * exp(-(x - 
              m1)^2/(2 * s1^2)), start = list(a1 = max(dens$y)/dnorm(mode,mode,mad(x)), 
              m1 = mode, s1 = mad(x)), data = list(x = dens$x[ind], 
              y = dens$y[ind]),control=list(warnOnly=warnOnly)))
        }
        if (plot) {
            lines(dens$x, coef(fit)[1]*
                  dnorm(dens$x, coef(fit)[2], coef(fit)[3]), col = "blue", lwd = 2)
            lines(dens$x, dens$y - coef(fit)[1]*
                  dnorm(dens$x, coef(fit)[2], coef(fit)[3]), col = "red", lwd = 2)
            lines(x = dens$x[ind], y = dens$y[ind], col = "orange",lwd=2)
            abline(h=0,lty=3)
        }
        summ<-try(coef(summary(fit)))
        if (class(summ)=='try-error') { 
          mu<-sd<-amp<-muerr<-sderr<-amperr<-NA
        } else { 
          summ<-coef(summary(fit))
          mu =summ['m1','Estimate']
          sd =summ['s1','Estimate']
          amp=summ['a1','Estimate']
          muerr =summ['m1','Std. Error']
          sderr =summ['s1','Std. Error']
          amperr=summ['a1','Std. Error']
        }
    }
    return(list(mu=mu,muerr=muerr,sd=sd,sderr=sderr,amp=amp,amperr=amperr,mode=mode,fit=fit))
}
