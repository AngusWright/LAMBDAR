convolvepsf <-
function(arr1, arr2, normalise=TRUE, mag2=TRUE, zapdig=NULL) {
  #Details {{{
  # Perform the Convolution of arr1 with arr2 using FFT
  # This convolution ignores all phase information from
  # the first image (arr1) - all positional information
  # is determined by the secondary image.
  # If noMag2 is TRUE, then all magnitude information
  # is determined by the first array
  # If normalise is true, then the final Fourier arrays
  # are normalised.
  # Use of normalise is recommended, as it allows for
  # additional error handling by scaling the zapsmall
  # function to the degree of overflow in the final
  # aperture }}}
  #If wanted, remove array 2 magnitude information {{{
  if (!mag2) {
     #Disregard Magnitude information
     #by normalisation
     arr2<-arr2*sum(arr2)/length(arr2)
  }#}}}
  #Perform Convolution {{{
  step<-Mod(fft(arr1))*fft(arr2)
  conv<-Re(fft(step,inverse=TRUE))
  #}}}
  #Normalise convolved array {{{
  if (normalise) {
    conv=conv/length(step)
  }#}}}
  #If needed, get the default zapdigit {{{
  if (is.null(zapdig)) { zapdig<-getOption("digits") }
  #}}}
  #'zap' small values to avoid numeric errors from fft, and return {{{
  return=zapsmall(conv, digits=zapdig)
  #}}}
}
