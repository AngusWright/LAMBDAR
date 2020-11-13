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
    RadialTolerance=25,
    MaxNumPSF=500,
    PSFCheck=0,
    PSFWeighted=1,
    PSFMap='NONE',
    Gauss_FWHM_AS=0,
    PSFLabel='SEEING',
    PSFLabelType='FWHM.ARCSEC',
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
    ApertureType=1,
    EFactor=1,
    FluxCorr=1,
    CropImage=0,
    CropFitsName='croppedimage',
    CropImRA0=-999,
    CropImDec0=-999,
    CropImRad=0.5,
    Confusion_units=0,
    nProcessors=1,
    CarefulWithEdges=1,
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
    PlotDevice='png',
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
    LDACBinary='ldactoasc',
    PyFITSRead=0, 
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
RootDirectory          ",paste(RootDirectory              ,collapse=' '),"
WorkingDirectory       ",paste(WorkingDirectory           ,collapse=' '),"
#                #-------File Paths: Relative to Root-------#
Catalogue              ",paste(Catalogue                  ,collapse=' '),"    #Catalogue Filename and Path
PSFMap                 ",paste(PSFMap                     ,collapse=' '),"    #PSF Filename and Path                      (use 'NONE' if not wanted)
DataMap                ",paste(DataMap                    ,collapse=' '),"    #Image Filename and Path
ErrorMap               ",paste(ErrorMap                   ,collapse=' '),"    #Error Map Filename and Path                (use 'NONE' if not wanted)
MaskMap                ",paste(MaskMap                    ,collapse=' '),"    #Mask Filename and Path                     (use 'NONE' if not wanted)
WeightMap              ",paste(WeightMap                  ,collapse=' '),"    #Weight Map Filename and Path               (use 'NONE' if not wanted)
#                #-------Catalogue Specific Variables-------#
CatIDColumnLabel       ",paste(CatIDColumnLabel           ,collapse=' '),"    #Label for designating the Catalogue ID Column of the supplied catalogue
RAColumnLabel          ",paste(RAColumnLabel              ,collapse=' '),"    #Label for designating the RA Column of the supplied catalogue
DecColumnLabel         ",paste(DecColumnLabel             ,collapse=' '),"    #Label for designating the Dec Column of the supplied catalogue
ThetaColumnLabel       ",paste(ThetaColumnLabel           ,collapse=' '),"    #Label for designating the Theta (on-sky aperture angle) Column of the supplied catalogue
SemiMajColumnLabel     ",paste(SemiMajColumnLabel         ,collapse=' '),"    #Label for designating the SemiMajor Axis Column of the supplied catalogue
SemiMinColumnLabel     ",paste(SemiMinColumnLabel         ,collapse=' '),"    #Label for designating the SemiMinor Axis Column of the supplied catalogue
ContamColumnLabel      ",paste(ContamColumnLabel          ,collapse=' '),"    #Label for designating the Contaminants ID Column of the supplied catalogue
FluxWgtColumnLabel     ",paste(FluxWgtColumnLabel         ,collapse=' '),"    #Label for designating the FluxWeights Column of the supplied catalogue
#                #-------Telescope Specific Variables-------#
MapUnitsPerBeam        ",paste(as.numeric(MapUnitsPerBeam),collapse=' '),"    #Is the Image in units per beam (as opposed to UnitsPerPixel)?            [1/0]
Confusion_units        ",paste(Confusion_units            ,collapse=' '),"    #Image Confusion in units
Gauss_FWHM_AS          ",paste(Gauss_FWHM_AS              ,collapse=' '),"    #Gaussian FWHM of Seeing in ArcSec
PSFLabel               ",paste(PSFLabel                   ,collapse=' '),"    #If Gaussian FWHM above is -ve, look for the PSF Sigma/FWHM in header with this keyword
PSFLabelType           ",paste(PSFLabelType               ,collapse=' '),"    #Is the PSFLabel keyword a Sigma or FWHM [SIGMA/FWHM] and in Pixels or Arcsec [PIX/AS]
SimGauss_AS            ",paste(SimGauss_AS                ,collapse=' '),"    #Gaussian FWHM of Gaussianisation Kernel (used for convolution of Simulation Image Noise), in ArcSec
BeamArea_SqAS          ",paste(BeamArea_SqAS              ,collapse=' '),"    #Beam area in Square Arcsec
FluxCorr               ",paste(FluxCorr                   ,collapse=' '),"    #Flux Correction Factor
EFactor                ",paste(EFactor                    ,collapse=' '),"    #Error Map Scale Factor
WeightMapZP            ",paste(WeightMapZP                ,collapse=' '),"    #Zero Point of the Weight Map - used for masking Data if no mask is supplied
Saturation             ",paste(Saturation                 ,collapse=' '),"    #Saturation value of the Map - used for Flagging bad fluxes. If unknown, use Inf
SaturationLabel        ",paste(SaturationLabel            ,collapse=' '),"    #Saturation Label in FITS header - used for reading saturation value if none supplied
GainLabel              ",paste(GainLabel                  ,collapse=' '),"    #Gain Label in FITS header - used for generating error map if no error.map/gain supplied
#                #-------Fits Extension Values--------------#
DataExtn               ",paste(DataExtn                   ,collapse=' '),"    #Fits header extension of Data in Image File
ErrorExtn              ",paste(ErrorExtn                  ,collapse=' '),"    #Fits header extension of Data in Error File
MaskExtn               ",paste(MaskExtn                   ,collapse=' '),"    #Fits header extension of Data in Mask File
WeightExtn             ",paste(WeightExtn                 ,collapse=' '),"    #Fits header extension of Data in Weight Map File
#                #--------Image Cropping Options------------#
CropImage              ",paste(as.numeric(CropImage)      ,collapse=' '),"    #Do we want to crop the input image to a particular subsection?                 [1/0]
CropImRA0              ",paste(CropImRA0                  ,collapse=' '),"    #RA that crop will centre around (deg; -999 will focus crop on image centre)
CropImDec0             ",paste(CropImDec0                 ,collapse=' '),"    #Dec that crop will centre around (deg; -999 will focus crop on image centre)
CropImRad              ",paste(CropImRad                  ,collapse=' '),"    #Radius of crop (deg; Negative will default to 1deg radius)
CropFitsName           ",paste(CropFitsName               ,collapse=' '),"    #Name of the output cropped image (without '.fits' extension)
#                #--------Additional Options----------------#
RemoveContam           ",paste(as.numeric(RemoveContam)   ,collapse=' '),"    #Remove Contaminants from the image?         [1/0]
CheckContam            ",paste(as.numeric(CheckContam)    ,collapse=' '),"    #Check Contaminants for ones that are irrelevant (do not overlap with main galaxies)?  [1/0]
nNearestCheck          ",paste(nNearestCheck              ,collapse=' '),"    #Number of Nearest Neighbours to check when looking for Contaminant-Galaxy Overlaps
RadialTolerance        ",paste(as.numeric(RadialTolerance),collapse=' '),"    #Radial tolerance (in pixels) for accepting sources in PSF estimation (default: 25pix)
MaxNPSF                ",paste(as.numeric(MaxNPSF)        ,collapse=' '),"    #Maximum Number of Point Sources used in PSF estimation. Lower numbers allow iterative cleaning of sources 
PSFCheck               ",paste(as.numeric(PSFCheck)       ,collapse=' '),"    #Check PSF by comparing to a stacked PSF made from the image?                   [1/0]
PSFConvolve            ",paste(as.numeric(PSFConvolve)    ,collapse=' '),"    #Convolve the Aperture with a PSF?                                              [1/0]
ApertureConfLimit      ",paste(ApertureConfLimit          ,collapse=' '),"    #Confidence limit used when converting PSF convolved apertures to binary apertures
PSFWeighted            ",paste(as.numeric(PSFWeighted)    ,collapse=' '),"    #Do you want PSFWeighted Photometry [1]? or Aperture Photometry [0]?            [1/0]
AngularOffset          ",paste(as.numeric(AngularOffset)  ,collapse=' '),"    #0 if the catalogue is in N0E90 angular coords, 1 if it is in N90E0 coords
ApertureType           ",paste(ApertureType               ,collapse=' '),"    #Choice of initial aperture type: 0 - all as point sources; 1 - boxcar aps; 2 - gaussian aps
MinApRad               ",paste(MinApRad                   ,collapse=' '),"    #State minimum aperture to use for sources
ResampleAper           ",paste(as.numeric(ResampleAper)   ,collapse=' '),"    #Smooth Boxcar Apertures by resampling? (Improves aperture surface integral accuracy)  [1/0]
ResamplingRes          ",paste(ResamplingRes              ,collapse=' '),"    #Rate of Upscaling in Resolution during Aperture Resampling (final res is n^[nIters] higher)
ResamplingIters        ",paste(ResamplingIters            ,collapse=' '),"    #Number of iterations to perform in resampling
UseMaskLim             ",paste(UseMaskLim                 ,collapse=' '),"    #Limit for determining whether or not to use an object overlapping the mask edge
Verbose                ",paste(as.numeric(Verbose)        ,collapse=' '),"    #Verbose Output?
ShowTime               ",paste(as.numeric(ShowTime)       ,collapse=' '),"    #Display execution & total elapsed time during run?
PlotSample             ",paste(as.numeric(PlotSample)     ,collapse=' '),"    #Plot a sample of the Object Apertures & Fluxes for Inspection? Includes Aperture Images and object COGs
PlotAll                ",paste(as.numeric(PlotAll)        ,collapse=' '),"    #Plot All of the Object Apertures & Fluxes for Inspection? [1] Or just the non-contaminant sources? [2] Includes Aperture Images and object COGs
Diagnostic             ",paste(as.numeric(Diagnostic)     ,collapse=' '),"    #Diagnostic Output of Variable Values During Computation - Helpful in Understanding Code
Interactive            ",paste(as.numeric(Interactive)    ,collapse=' '),"    #Interactive Mode - Allows user access to entire parameter-space after final calculations
Magnitudes             ",paste(as.numeric(Magnitudes)     ,collapse=' '),"    #Output Fluxes as magnitudes
ABVegaFlux             ",paste(ABVegaFlux                 ,collapse=' '),"    #Flux of ABVega in this band - only used if outputting magnitudes
MagZeroPoint           ",paste(MagZeroPoint               ,collapse=' '),"    #Magnitude of the Image Zero Point. If -999, ZP will be read from FITS header
MagZPLabel             ",paste(MagZPLabel                 ,collapse=' '),"    #Label used for the Zero Point Keyword in the FITS header
DoSkyEst               ",paste(as.numeric(DoSkyEst)       ,collapse=' '),"    #Perform estimate of local sky-background for each object, and subtract it from the flux
SkyEstIters            ",paste(SkyEstIters                ,collapse=' '),"    #Number of iterations of sigma-cutting used in sky estimation
SkyEstProbCut          ",paste(SkyEstProbCut              ,collapse=' '),"    #Sigma Level used in sigma-cutting of sky pixels
SkyDefault             ",paste(SkyDefault                 ,collapse=' '),"    #Default Value to use for local sky if estimation fails. May be numeric, 'median', or 'mean'. Otherwise will be NA.
SkyCorrelNoise         ",paste(SkyCorrelNoise             ,collapse=' '),"    #Level of Correlated noise in the image, if known (factor is multiplicative - 1 == no correlated noise)
GetSkyRMS              ",paste(as.numeric(GetSkyRMS)      ,collapse=' '),"    #As above without subtraction, and output the local sky RMS (if doing sky estimate, this happens automatically)
RanCor                 ",paste(as.numeric(RanCor)         ,collapse=' '),"    #Do you want to calculate mean & median sky-flux using randoms around every aperture? [1/0]
nRandoms               ",paste(nRandoms                   ,collapse=' '),"    #Number of randoms calculated *per aperture*
BlankCor               ",paste(as.numeric(BlankCor)       ,collapse=' '),"    #Do you want to calculate mean & median confusion using blanks around every aperture? [1/0]
nBlanks                ",paste(nBlanks                    ,collapse=' '),"    #Number of blanks calculated *per aperture*
UsePixelFluxWgts       ",paste(UsePixelFluxWgts           ,collapse=' '),"    #Do you want the pixel flux at the object RA DEC to be used for relative flux.weighting? [1/0]
FluxWgtType            ",paste(FluxWgtType                ,collapse=' '),"    #What is the form of the input flux.weights ('flux', 'mag', or 'scale')?
IterateFluxes          ",paste(as.numeric(IterateFluxes)  ,collapse=' '),"    #Do you want to iterate the flux determination to improve deblending?
nIterations            ",paste(nIterations                ,collapse=' '),"    #How many iterations do you want to do?
#                #---------------Outputs--------------------#
OutputDirectory        ",paste(OutputDirectory            ,collapse=' '),"    #Output directory Name and Path
LogFile                ",paste(LogFile                    ,collapse=' '),"    #Filename for Log
PlotDevice             ",paste(PlotDevice                 ,collapse=' '),"    #The device to use for plotting: can be one of: PNG, PDF, or X11
GetDeblFrac            ",paste(as.numeric(GetDeblFrac)    ,collapse=' '),"    #Only Output Deblend Fractions for each aperture to file? (Overrules all other Flags)            [1/0]
SourceMaskOnly         ",paste(as.numeric(SourceMaskOnly) ,collapse=' '),"    #Only Output Source Mask image to file? (Overrules all other Flags except the above)            [1/0]
WriteTable             ",paste(as.numeric(WriteTable)     ,collapse=' '),"    #Output Results Table?
TableName              ",paste(TableName                  ,collapse=' '),"    #Filename for Results Table
WriteAAMask            ",paste(as.numeric(WriteAAMask)    ,collapse=' '),"    #Write the All Apertures Mask image to file?                                   [1/0]
AllApersFile           ",paste(AllApersFile               ,collapse=' '),"    #Filename for All Apertures Mask
WriteFAMask            ",paste(as.numeric(WriteFAMask)    ,collapse=' '),"    #Write the Convolved Apertures Mask image to file?                             [1/0]
ConvApersFile          ",paste(ConvApersFile              ,collapse=' '),"    #Filename for All Convolved Apertures Mask
WriteDFAMask           ",paste(as.numeric(WriteDFAMask)   ,collapse=' '),"    #Write the Deblended Convolved Apertures Mask image to file?                   [1/0]
DeblConvApersFile      ",paste(DeblConvApersFile          ,collapse=' '),"    #Filename for Deblended Convolved Apertures Mask
WriteSourceMask        ",paste(as.numeric(WriteSourceMask),collapse=' '),"    #Write the Source Mask image to file?                                          [1/0]
SourceMaskConfLim      ",paste(SourceMaskConfLim          ,collapse=' '),"    #When convolving with PSF, to what confidence should the Sourcemask extend when masking Apertures?
SourceMaskFile         ",paste(SourceMaskFile             ,collapse=' '),"    #Filename for Source Mask
TransmissionMap        ",paste(as.numeric(TransmissionMap),collapse=' '),"    #Do you want to produce a Transmission Map instead of a Sourcemask (transparent over sources instead of opaque)? [0/1]
WriteResidMap          ",paste(as.numeric(WriteResidMap)  ,collapse=' '),"    #Write the Residual Image Map to file?                                         [1/0]
ResidImageFile         ",paste(ResidImageFile             ,collapse=' '),"    #Filename for Residual Image
NoContamImageFile      ",paste(NoContamImageFile          ,collapse=' '),"    #Filename for Contaminant Subtracted Residual Image
#                #-------Computational Management-----------#
MemorySafe             ",paste(as.numeric(MemorySafe)     ,collapse=' '),"    #Do we want to perform checks to ensure the program does not use more memory than is available? [1/0]
nProcessors            ",paste(nProcessors                ,collapse=' '),"    #Number of Processors Available for Use in Computations
CarefulWithEdges      ",paste(as.numeric(CarefulWithEdges),collapse=' '),"    #Do we want to be careful with sources at the edge of the image? Things too close to image edges are removed
ApStampWidth           ",paste(ApStampWidth               ,collapse=' '),"    #Width of the Aperture stamps in multiples of aperture major axis; Can be changed with caution if memory issues arise.
PSFConfidence          ",paste(PSFConfidence              ,collapse=' '),"    #PSF Confidence Value used in buffering the Aperture Stamp Widths; PSF integrated out to this width, and that radii is added.
LDACBinary             ",paste(LDACBinary                 ,collapse=' '),"    #If using LDAC catalogue format, specify here the absolute path to the LDAC-to-ASCII conversion binary.
PyFITSRead             ",paste(PyFITSRead                 ,collapse=' '),"    #Do you want to read images via PyFITS instead of via the FITSio package (with large images can be ~2x faster."))
#}}}
#Close Sink and return NULL {{{
sink()
return=NULL
#}}}
}
