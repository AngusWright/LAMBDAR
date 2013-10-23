convolvepsf <-
function(arr1, arr2, normalise=TRUE, nomag2=FALSE) {
  # Perform the Convolution of arr1 with arr2 using FFT
  # This convolution ignores all phase information from 
  # the first image (arr1) - all positional information 
  # is determined by the secondary image.
  # If noMag2 is TRUE, then all magnitude information 
  # is determined by the first array
  # normalise
  if (normalise) {
    if (!nomag2) {
      step<-Mod(fft(arr1))*fft(arr2)
      conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step)))
    } else {
      step<-Mod(fft(arr1))*(fft(arr2)/Mod(fft(arr2)))
      conv<-Re(zapsmall(fft(step,inverse=TRUE)/length(step)))
    }
  } else {
    if (!nomag2) {
      conv<-Re(zapsmall(fft(Mod(fft(arr1))*fft(arr2),inverse=TRUE)))
    } else {
      conv<-Re(zapsmall(fft(Mod(fft(arr1))*(fft(arr2)/Mod(fft(arr2))),inverse=TRUE)))
    }
  }
  return(conv)
}
