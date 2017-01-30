#
# Simple wrapper to CairoPNG and png that 
# determines which to use:
#   Preferentially CairoPNG if available.
#

PlotPNG<-function(...) {
   
  #If Cairo is available
  if (requireNamespace('Cairo',quietly=TRUE)) {
  #Open CairoPNG
    opened=Cairo::CairoPNG(...)
  } else {
  #Open png
    opened=png(...)
  }
  return=opened
}
