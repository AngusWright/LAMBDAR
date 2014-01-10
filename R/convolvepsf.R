convolvepsf <-
function(arr1, arr2, normalise=TRUE, nomag2=FALSE, zapdig=NULL) {
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
  # aperture
  if (nomag2) {
     #Disregard Magnitude information
     #by normalisation
     arr2<-arr2*sum(arr2)/length(arr2)
  }
  step<-Mod(fft(arr1))*fft(arr2)
  if (is.null(zapdig)) { zapdig<-getOption("digits") }
  conv<-Re(fft(step,inverse=TRUE))
  #message(paste("Max values of Input Array 1 and Convolved Array:", max(arr1),max(conv)))
  #if (abs(max(conv))-abs(max(arr1)) > 10^(-zapdig)) {
  #  zapdig<-floor(-log10(abs(max(conv))-abs(max(arr1))))
    #print(paste("Rezapping with Zapdigit:",zapdig,"(",max(conv),max(arr1),")"))
  #}
  if (normalise) {
    conv=conv/length(step)
  }
  return=zapsmall(conv, digits=zapdig)
}
