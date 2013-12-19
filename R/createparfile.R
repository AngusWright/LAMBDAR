createparfile <-
  function(
    RootDirectory="./",
    Catalogue="./catalogue.fits",
    PSFMap="NONE",
    DataMap="./datamap.fits",
    ErrorMap="NONE",
    MaskMap="NONE",
    MapJyPerBeam=0,
    Confusion_Jy=0.0,
    Gauss_FWHM_AS=0.0,
    BeamArea_SqAS=0.0,
    FluxCorr=1.0,
    EFactor=1.0,
    PSFheight=1.0,
    NormalisePSF=1,
    DataExtn=0,
    ErrorExtn=0,
    MaskExtn=0,
    CropImage=0,
    CropImRA0=-999,
    CropImDec0=-999,
    CropImRad=NA,
    CropFitsName="cropped_im",
    RemoveContam=0,
    PSFConvolve=1,
    AngularOffset=0,
    UseMask=0,
    PointSources=0,
    SmoothAper=1,
    ResamplingRes=3,
    ResamplingIters=4,
    Verbose=0,
    ShowTime=1,
    PlotSample=0,
    Diagnostic=0,
    Interactive=0,
    Magnitudes=0,
    ABVegaFlux=1.0,
    MagZeroPoint=-999,
    MagZPLabel="MagZP",
    DoSkyEst=0,
    OutputDirectory="./Output/",
    LogFile="LAMBDARlog.txt",
    SourceMaskOnly=0,
    WriteTable=1,
    TableName="dfaResults.fits",
    WriteAAMask=1,
    AllApersFile="AllAps_Mask.fits",
    OverlayEllipse=1,
    OverlayFile="OverlayEllipse.fits",
    WriteSourceMask=1,
    SourceMaskFile="SourceMask.fits",
    WriteDFAMask=1,
    DelbConvApersFile="AllDbCnvAps_Mask.fits",
    WriteFAMask=1,
    ConvApersFile="AllConvAps_Mask.fits",
    WriteResidMap=1,
    ResidImageFile="ResidualImage.fits",
    NoContamImage="NoContamResIm.fits",
    nProcessors=1,
    ApStampWidth=1.05,
    StampMult=3 ) {

sink(file="Lambdar_default.par")
cat(paste("#-------------------------------------------------------------------------------------
#
#                   Default Parameter File for use by LAMBDAR
#
#   NB: All Parameters are provided here.
#       In some cases, certain parameters will not be required.
#
# Written by:  A.H. Wright
#
#-------------------------------------------------------------------------------------

#                #---------Working Directory Path-----------#
RootDirectory          ",RootDirectory     ,"
#                #-------File Paths: Relative to Root-------#
Catalogue              ",Catalogue         ,"    #Catalogue Filename and Path
PSFMap                 ",PSFMap            ,"    #PSF Filename and Path                      (use 'NONE' if not wanted)
DataMap                ",DataMap           ,"    #Image Filename and Path
ErrorMap               ",ErrorMap          ,"    #Error Map Filename and Path                (use 'NONE' if not wanted)
MaskMap                ",MaskMap           ,"    #Mask Filename and Path                     (use 'NONE' if not wanted)
#                #-------Telescope Specific Variables-------#
MapJyPerBeam           ",MapJyPerBeam      ,"    #Is the Image Flux measured in Jyperbeam (as opposed to JyPerPixel)?            [1/0]
Confusion_Jy           ",Confusion_Jy      ,"    #Confusion in Jy
Gauss_FWHM_AS          ",Gauss_FWHM_AS     ,"    #Gaussian FWHM of Seeing in ArcSec
BeamArea_SqAS          ",BeamArea_SqAS     ,"    #Beam area in Square Arcsec
FluxCorr               ",FluxCorr          ,"    #Flux Correction Factor
EFactor                ",EFactor           ,"    #Error Map Scale Factor
NormalisePSF           ",NormalisePSF      ,"    #Normalise the PSF to Unity?                                                    [1/0]
PSFheight              ",PSFheight         ,"    #PSF Scale-height
#                #-------Fits Extension Values--------------#
DataExtn               ",DataExtn          ,"    #Fits header extension of Data in Image File
ErrorExtn              ",ErrorExtn         ,"    #Fits header extension of Data in Error File
MaskExtn               ",MaskExtn          ,"    #Fits header extension of Data in Mask File
#                #--------Image Cropping Options------------#
CropImage             ",CropImage          ,"    #Do we want to crop the input image to a particular subsection?                 [1/0]
CropImRA0             ",CropImRA0          ,"    #RA that crop will centre around (deg; NA will focus crop on image centre)
CropImDec0            ",CropImDec0         ,"    #Dec that crop will centre around (deg; NA will focus crop on image centre)
CropImRad             ",CropImRad          ,"    #Radius of crop (deg; NA will default to 1deg radius)
CropFitsName          ",CropFitsName       ,"    #Name of the output cropped image (without '.fits' extension)
#                #--------Additional Options----------------#
RemoveContam           ",RemoveContam      ,"    #Remove Contaminants (Catalogue must contain 'CONTAMFLAG' column) from the image?         [1/0]
PSFConvolve            ",PSFConvolve       ,"    #Convolve the Aperture with a PSF?                                              [1/0]
AngularOffset          ",AngularOffset     ,"    #0 if the catalogue is in N90E0 angular coords, 1 if it is in N0E90 coords
UseMask                ",UseMask           ,"    #Use the Source Mask?                                                           [1/0]
PointSources           ",PointSources      ,"    #Force point sources to be used?                                                [1/0]
SmoothAper             ",SmoothAper        ,"    #Smooth Apertures by resampling? (Improves aperture surface integral accuracy)  [1/0]
ResamplingRes          ",ResamplingRes     ,"    #Rate of Upscaling in Resolution during Aperture Resampling (final res is n^[nIters] higher)
ResamplingIters        ",ResamplingIters   ,"    #Number of iterations to perform in resampling
Verbose                ",Verbose           ,"    #Verbose Output?
ShowTime               ",ShowTime          ,"    #Display execution & total elapsed time during run?
PlotSample             ",PlotSample        ,"    #Plot a sample of the Single Apertures for Inspection?
Diagnostic             ",Diagnostic        ,"    #Diagnostic Output of Variable Values During Computation - Helpful in Understanding Code
Interactive            ",Interactive       ,"    #Interactive Mode - Allows user access to entire parameter-space after final calculations
Magnitudes             ",Magnitudes        ,"    #Output Fluxes as ABVega Magnitudes
ABVegaFlux             ",ABVegaFlux        ,"    #Flux of ABVega in this band - only used if outputting Magnitudes
MagZeroPoint           ",MagZeroPoint      ,"    #Magnitude of the Image Zero Point. If -999, ZP will be read from FITS header
MagZPLabel             ",MagZPLabel        ,"    #Label used for the Zero Point Keyword in the FITS header
DoSkyEst               ",DoSkyEst          ,"    #Perform estimate of local sky-background for each object, and subtract it from the flux
#                #---------------Outputs--------------------#
OutputDirectory        ",OutputDirectory   ,"    #Output directory Name and Path
LogFile                ",LogFile           ,"    #Filename for Log
SourceMaskOnly         ",SourceMaskOnly    ,"    #Only Output Source Mask image to file? (Overrules all other Flags)            [1/0]
WriteTable             ",WriteTable        ,"    #Output Results Table?
TableName              ",TableName         ,"    #Filename for Results Table
WriteAAMask            ",WriteAAMask       ,"    #Write the All Apertures Mask image to file?                                   [1/0]
AllApersFile           ",AllApersFile      ,"    #Filename for All Apertures Mask
OverlayEllipse         ",OverlayEllipse    ,"    #Write Ellipse Overlaid Image to file?                                         [1/0]
OverlayFile            ",OverlayFile       ,"    #Filename for All Apertures Mask
WriteSourceMask        ",WriteSourceMask   ,"    #Write the Source Mask image to file?                                          [1/0]
SourceMaskFile         ",SourceMaskFile    ,"    #Filename for Source Mask
WriteDFAMask           ",WriteDFAMask      ,"    #Write the Convolved Apertures Mask image to file?                             [1/0]
DelbConvApersFile      ",DelbConvApersFile ,"    #Filename for All Convolved Apertures Mask
WriteFAMask            ",WriteFAMask       ,"    #Write the Convolved Apertures Mask image to file?                             [1/0]
ConvApersFile          ",ConvApersFile     ,"    #Filename for All Convolved Apertures Mask
WriteResidMap          ",WriteResidMap     ,"    #Write the Residual Image Map to file?                                         [1/0]
ResidImageFile         ",ResidImageFile    ,"    #Filename for Residual Image
NoContamImage          ",NoContamImage     ,"    #Filename for Contaminant Subtracted Residual Image
#                #-------Computational Management-----------#
nProcessors            ",nProcessors       ,"    #Number of Processors Available for Use in Computations
ApStampWidth           ",ApStampWidth      ,"    #Width of the Aperture stamps in multiples of aperture major axis; Can be changed with caution if memory issues arise.
StampMult              ",StampMult         ,"    #Number of PSF FWHM's added to the Aperture Stamp Widths, to wings on apertures after PSF convolution."))

sink()
}
