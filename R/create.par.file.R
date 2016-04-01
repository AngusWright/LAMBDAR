create.par.file <-
  function(
    filename="Lambdar_default.par",
    RootDirectory="./",
    WorkingDirectory="./",
    OutputDirectory="./",
    Catalogue="catalogue.csv",
    DataMap="image.fits",
    BeamArea_SqAS=0,
    PSFConvolve=0,
    PSFWeighted=1,
    PSFMap='NONE',
    Gauss_FWHM_AS=0,
    SimGauss_AS=0,
    RemoveContam=0,
    CheckContam=0,
    nNearestCheck=10,
    NoContamImageFile='NoContamResidualImage.fits',
    CatIDColumnLabel='CATAID',
    RAColumnLabel='RAdeg',
    DecColumnLabel='DECdeg',
    ThetaColumnLabel='rotN2E',
    SemiMajColumnLabel='radmajasec',
    SemiMinColumnLabel='radminasec',
    ContamColumnLabel='CONTAM',
    FluxWgtColumnLabel='FLUXWEIGHT',
    ErrorMap='NONE',
    MaskMap='NONE',
    WeightMap='NONE',
    WeightMapZP=0,
    Saturation=Inf,
    SaturationLabel='SATUR',
    GainLabel='GAIN',
    DataExtn=0,
    ErrorExtn=0,
    MaskExtn=0,
    WeightExtn=0,
    PointSources=0,
    EFactor=0,
    FluxCorr=1,
    CropImage=0,
    CropFitsName='croppedimage',
    CropImRA0=-999,
    CropImDec0=-999,
    CropImRad=0.5,
    Confusion_units=0,
    nProcessors=1,
    AngularOffset=0,
    MapUnitsPerBeam=0,
    ResampleAper=0,
    ResamplingRes=3,
    ResamplingIters=2,
    PSFConfidence=1,
    ApStampWidth=1.05,
    SourceMaskOnly=0,
    WriteSourceMask=0,
    WriteAAMask=0,
    AllApersFile='AllApertures_Mask.fits',
    WriteFAMask=0,
    ConvApersFile='AllConvolvedApertures_Mask.fits',
    WriteDFAMask=0,
    DeblConvApersFile='AllDeblConvolvedApertures_Mask.fits',
    WriteResidMap=0,
    ResidImageFile='ResidualImage.fits',
    WriteTable=1,
    TableName='LAMBDAR_Results',
    ShowTime=0,
    Interactive=0,
    UseMaskLim=0.2,
    Diagnostic=0,
    Verbose=0,
    PlotSample=0,
    PlotAll=0,
    Magnitudes=1,
    MagZPLabel='MagZP',
    ABVegaFlux=1.0,
    MagZeroPoint=8.9,
    BlankCor=0,
    nBlanks=10,
    RanCor=0,
    nRandoms=10,
    DoSkyEst=0,
    GetSkyRMS=0,
    SkyEstIters=5,
    SkyEstProbCut=3,
    SkyCorrelNoise=1,
    SkyDefault='median',
    SourceMaskFile='SourceMask.fits',
    TransmissionMap=0,
    GetDeblFrac=0,
    SourceMaskConfLim=0.95,
    MinApRad=0,
    MemorySafe=0,
    ApertureConfLimit=0.9,
    IterateFluxes=0,
    nIterations=2,
    FluxWgtType='scale',
    UsePixelFluxWgts=0,
    LogFile='LAMBDAR_Log.txt') {

#Sink Output to File {{{
sink(file=filename)
#}}}
#Write File {{{
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
WorkingDirectory       ",WorkingDirectory  ,"
#                #-------File Paths: Relative to Root-------#
Catalogue              ",Catalogue         ,"    #Catalogue Filename and Path
PSFMap                 ",PSFMap            ,"    #PSF Filename and Path                      (use 'NONE' if not wanted)
DataMap                ",DataMap           ,"    #Image Filename and Path
ErrorMap               ",ErrorMap          ,"    #Error Map Filename and Path                (use 'NONE' if not wanted)
MaskMap                ",MaskMap           ,"    #Mask Filename and Path                     (use 'NONE' if not wanted)
WeightMap              ",WeightMap         ,"    #Weight Map Filename and Path               (use 'NONE' if not wanted)
#                #-------Catalogue Specific Variables-------#
CatIDColumnLabel       ",CatIDColumnLabel  ,"    #Label for designating the Catalogue ID Column of the supplied catalogue
RAColumnLabel          ",RAColumnLabel     ,"    #Label for designating the RA Column of the supplied catalogue
DecColumnLabel         ",DecColumnLabel    ,"    #Label for designating the Dec Column of the supplied catalogue
ThetaColumnLabel       ",ThetaColumnLabel  ,"    #Label for designating the Theta (on-sky aperture angle) Column of the supplied catalogue
SemiMinColumnLabel     ",SemiMinColumnLabel,"    #Label for designating the SemiMinor Axis Column of the supplied catalogue
SemiMajColumnLabel     ",SemiMajColumnLabel,"    #Label for designating the SemiMajor Axis Column of the supplied catalogue
ContamColumnLabel      ",ContamColumnLabel ,"    #Label for designating the Contaminants ID Column of the supplied catalogue
FluxWgtColumnLabel     ",FluxWgtColumnLabel,"    #Label for designating the FluxWeights Column of the supplied catalogue
#                #-------Telescope Specific Variables-------#
MapUnitsPerBeam        ",MapUnitsPerBeam      ,"    #Is the Image in units per beam (as opposed to UnitsPerPixel)?            [1/0]
Confusion_units        ",Confusion_units   ,"    #Image Confusion in units
Gauss_FWHM_AS          ",Gauss_FWHM_AS     ,"    #Gaussian FWHM of Seeing in ArcSec
SimGauss_AS            ",SimGauss_AS       ,"    #Gaussian FWHM of Gaussianisation Kernel (used for convolution of Simulation Image Noise), in ArcSec
BeamArea_SqAS          ",BeamArea_SqAS     ,"    #Beam area in Square Arcsec
FluxCorr               ",FluxCorr          ,"    #Flux Correction Factor
EFactor                ",EFactor           ,"    #Error Map Scale Factor
WeightMapZP            ",WeightMapZP       ,"    #Zero Point of the Weight Map - used for masking Data if no mask is supplied
Saturation             ",Saturation        ,"    #Saturation value of the Map - used for Flagging bad fluxes. If unknown, use Inf
SaturationLabel        ",SaturationLabel   ,"    #Saturation Label in FITS header - used for reading saturation value if none supplied
GainLabel              ",GainLabel         ,"    #Gain Label in FITS header - used for generating error map if no error.map/gain supplied
#                #-------Fits Extension Values--------------#
DataExtn               ",DataExtn          ,"    #Fits header extension of Data in Image File
ErrorExtn              ",ErrorExtn         ,"    #Fits header extension of Data in Error File
MaskExtn               ",MaskExtn          ,"    #Fits header extension of Data in Mask File
WeightExtn             ",WeightExtn        ,"    #Fits header extension of Data in Weight Map File
#                #--------Image Cropping Options------------#
CropImage              ",CropImage          ,"    #Do we want to crop the input image to a particular subsection?                 [1/0]
CropImRA0              ",CropImRA0          ,"    #RA that crop will centre around (deg; NA will focus crop on image centre)
CropImDec0             ",CropImDec0         ,"    #Dec that crop will centre around (deg; NA will focus crop on image centre)
CropImRad              ",CropImRad          ,"    #Radius of crop (deg; NA will default to 1deg radius)
CropFitsName           ",CropFitsName       ,"    #Name of the output cropped image (without '.fits' extension)
#                #--------Additional Options----------------#
RemoveContam           ",RemoveContam      ,"    #Remove Contaminants from the image?         [1/0]
CheckContam            ",CheckContam       ,"    #Check Contaminants for ones that are irrelevant (do not overlap with main galaxies)?  [1/0]
nNearestCheck          ",nNearestCheck     ,"    #Number of Nearest Neighbours to check when looking for Contaminant-Galaxy Overlaps
PSFConvolve            ",PSFConvolve       ,"    #Convolve the Aperture with a PSF?                                              [1/0]
ApertureConfLimit      ",ApertureConfLimit ,"    #Confidence limit used when converting PSF convolved apertures to binary apertures
PSFWeighted            ",PSFWeighted       ,"    #Do you want PSFWeighted Photometry [1]? or Aperture Photometry [0]?            [1/0]
AngularOffset          ",AngularOffset     ,"    #0 if the catalogue is in N0E90 angular coords, 1 if it is in N90E0 coords
PointSources           ",PointSources      ,"    #Force point sources to be used?                                                [1/0]
MinApRad               ",MinApRad          ,"    #State minimum aperture to use for sources
ResampleAper           ",ResampleAper      ,"    #Smooth Apertures by resampling? (Improves aperture surface integral accuracy)  [1/0]
ResamplingRes          ",ResamplingRes     ,"    #Rate of Upscaling in Resolution during Aperture Resampling (final res is n^[nIters] higher)
ResamplingIters        ",ResamplingIters   ,"    #Number of iterations to perform in resampling
UseMaskLim             ",UseMaskLim        ,"    #Limit for determining whether or not to use an object overlapping the mask edge
Verbose                ",Verbose           ,"    #Verbose Output?
ShowTime               ",ShowTime          ,"    #Display execution & total elapsed time during run?
PlotSample             ",PlotSample        ,"    #Plot a sample of the Object Apertures & Fluxes for Inspection? Includes Aperture Images and object COGs
PlotAll                ",PlotAll           ,"    #Plot All of the Object Apertures & Fluxes for Inspection? Includes Aperture Images and object COGs
Diagnostic             ",Diagnostic        ,"    #Diagnostic Output of Variable Values During Computation - Helpful in Understanding Code
Interactive            ",Interactive       ,"    #Interactive Mode - Allows user access to entire parameter-space after final calculations
Magnitudes             ",Magnitudes        ,"    #Output Fluxes as magnitudes
ABVegaFlux             ",ABVegaFlux        ,"    #Flux of ABVega in this band - only used if outputting magnitudes
MagZeroPoint           ",MagZeroPoint      ,"    #Magnitude of the Image Zero Point. If -999, ZP will be read from FITS header
MagZPLabel             ",MagZPLabel        ,"    #Label used for the Zero Point Keyword in the FITS header
DoSkyEst               ",DoSkyEst          ,"    #Perform estimate of local sky-background for each object, and subtract it from the flux
SkyEstIters            ",SkyEstIters       ,"    #Number of iterations of sigma-cutting used in sky estimation
SkyEstProbCut          ",SkyEstProbCut     ,"    #Sigma Level used in sigma-cutting of sky pixels
SkyDefault             ",SkyDefault        ,"    #Default Value to use for local sky if estimation fails. May be numeric, 'median', or 'mean'. Otherwise will be NA.
SkyCorrelNoise         ",SkyCorrelNoise    ,"    #Level of Correlated noise in the image, if known (factor is multiplicative - 1 == no correlated noise)
GetSkyRMS              ",GetSkyRMS         ,"    #As above without subtraction, and output the local sky RMS (if doing sky estimate, this happens automatically)
RanCor                 ",RanCor            ,"    #Do you want to calculate mean & median sky-flux using randoms around every aperture? [1/0]
nRandoms               ",nRandoms          ,"    #Number of randoms calculated *per aperture*
BlankCor               ",BlankCor          ,"    #Do you want to calculate mean & median confusion using blanks around every aperture? [1/0]
nBlanks                ",nBlanks           ,"    #Number of blanks calculated *per aperture*
UsePixelFluxWgts       ",UsePixelFluxWgts  ,"    #Do you want the pixel flux at the object RA DEC to be used for relative flux.weighting? [1/0]
FluxWgtType            ",FluxWgtType       ,"    #What is the form of the input flux.weights ('flux', 'mag', or 'scale')?
IterateFluxes          ",IterateFluxes     ,"    #Do you want to iterate the flux determination to improve deblending?
nIterations            ",nIterations       ,"    #How many iterations do you want to do?
#                #---------------Outputs--------------------#
OutputDirectory        ",OutputDirectory   ,"    #Output directory Name and Path
LogFile                ",LogFile           ,"    #Filename for Log
GetDeblFrac            ",GetDeblFrac       ,"    #Only Output Deblend Fractions for each aperture to file? (Overrules all other Flags)            [1/0]
SourceMaskOnly         ",SourceMaskOnly    ,"    #Only Output Source Mask image to file? (Overrules all other Flags except the above)            [1/0]
WriteTable             ",WriteTable        ,"    #Output Results Table?
TableName              ",TableName         ,"    #Filename for Results Table
WriteAAMask            ",WriteAAMask       ,"    #Write the All Apertures Mask image to file?                                   [1/0]
AllApersFile           ",AllApersFile      ,"    #Filename for All Apertures Mask
WriteFAMask            ",WriteFAMask       ,"    #Write the Convolved Apertures Mask image to file?                             [1/0]
ConvApersFile          ",ConvApersFile     ,"    #Filename for All Convolved Apertures Mask
WriteDFAMask           ",WriteDFAMask      ,"    #Write the Deblended Convolved Apertures Mask image to file?                   [1/0]
DeblConvApersFile      ",DeblConvApersFile ,"    #Filename for Deblended Convolved Apertures Mask
WriteSourceMask        ",WriteSourceMask   ,"    #Write the Source Mask image to file?                                          [1/0]
SourceMaskConfLim      ",SourceMaskConfLim ,"    #When convolving with PSF, to what confidence should the Sourcemask extend when masking Apertures?
SourceMaskFile         ",SourceMaskFile    ,"    #Filename for Source Mask
TransmissionMap        ",TransmissionMap   ,"    #Do you want to produce a Transmission Map instead of a Sourcemask (transparent over sources instead of opaque)? [0/1]
WriteResidMap          ",WriteResidMap     ,"    #Write the Residual Image Map to file?                                         [1/0]
ResidImageFile         ",ResidImageFile    ,"    #Filename for Residual Image
NoContamImageFile      ",NoContamImageFile ,"    #Filename for Contaminant Subtracted Residual Image
#                #-------Computational Management-----------#
MemorySafe             ",MemorySafe        ,"    #Do we want to perform checks to ensure the program does not use more memory than is available? [1/0]
nProcessors            ",nProcessors       ,"    #Number of Processors Available for Use in Computations
ApStampWidth           ",ApStampWidth      ,"    #Width of the Aperture stamps in multiples of aperture major axis; Can be changed with caution if memory issues arise.
PSFConfidence          ",PSFConfidence     ,"    #PSF Confidence Value used in buffering the Aperture Stamp Widths; PSF integrated out to this width, and that radii is added."))
#}}}
#Close Sink and return NULL {{{
sink()
return=NULL
#}}}
}
