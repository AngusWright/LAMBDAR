convolvepsf <-
function(arr1, arr2, normalise=TRUE, nomag2=FALSE) {
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
  if (normalise) {
    if (!nomag2) {
      step<-Mod(fft(arr1))*fft(arr2)
      conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step)))
      zapdig<-getOption("digits")
      if (abs(max(conv))-abs(max(arr1)) > 10^(-zapdig)) {
        zapdig<-floor(-log10(abs(max(conv))-abs(max(arr1))))
        #print(paste("Rezapping with Zapdigit:",zapdig,"(",max(conv),max(arr1),")"))
        conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step), digits=zapdig))
      }
    } else {
      step<-Mod(fft(arr1))*(fft(arr2)/Mod(fft(arr2)))
      conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step)))
      zapdig<-getOption("digits")
      if (abs(max(conv))-abs(max(arr1)) > 10^(-zapdig)) {
        zapdig<-floor(-log10(abs(max(conv))-abs(max(arr1))))
        #print(paste("Rezapping with Zapdigit:",zapdig,"(",(abs(max(conv))-abs(max(arr1))),")"))
        conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step), digits=zapdig))
      }
    }
  } else {
    if (!nomag2) {
      conv<-Re(zapsmall(fft(Mod(fft(arr1))*fft(arr2),inverse=TRUE)))
    } else {
      conv<-Re(zapsmall(fft(Mod(fft(arr1))*(fft(arr2)/Mod(fft(arr2))),inverse=TRUE)))
    }
  }
  return=conv
}
