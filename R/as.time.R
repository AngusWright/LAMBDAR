#As Time functions /*fold*/ {{{
as.time<-function(sec,digits) {
  if (sec > 60*60*24) {
    day<-floor(sec/(60*60*24))
    hr<-floor((sec-day*(60*60*24))/(60*60))
    min<-round((sec-day*(60*60*24)-hr*(60*60))/60,digits=digits)
    timestr<-paste(day,'day',hr,'hr',min,'min')
  } else if (sec > 60*60) {
    hr<-floor(sec/(60*60))
    min<-floor((sec-hr*(60*60))/60)
    sec<-round(sec-hr*(60*60)-min*60,digits=digits)
    timestr<-paste(hr,'hr',min,'min',sec,'sec')
  } else if (sec > 60 ) {
    min<-floor(sec/60)
    sec<-round(sec-min*60,digits=digits)
    timestr<-paste(min,'min',sec,'sec')
  } else {
    sec<-round(sec,digits=digits)
    timestr<-paste(sec,'sec')
  }
  return=timestr
}
#/*fend*/ }}}
