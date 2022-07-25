version = '0.5.1'
date    = '2022-07-14 12:00:00 GMT'
yaml_descriptions = """
- version : 0.5.1
  date : 2022-07-12
  comment : Update of DLC server URL to scopesim.univie.ac.at
  changes :
  - Changed URL in defaults.yaml file

- version : 0.5.0
  date : 2022-04-22
  comment : IFU Spectroscopy mode for METIS
  changes :
  - The IFU effects for the METIS LMS mode
  - Effects for including extra FITS keywords for the MICADO-ESO delivery
  - Minor change to the OpticalTrain.readout method to allow custom FITS keywords to be added
  - #-strings for accessing the .meta dict contents of Effect objects
  - added the tests.mocks.basic_instrument package for test purposes
  - refactored package downloading functions in server.database 
  - Packages can now be downloaded directly from a git commit 
  - new RTDs structure for docs based on ipynb files
  - change to SkyCalcTERCurve to use local files to avoid calling the skycalc server  
  - New Effects:
    - MetisLMSSpectralTraceList(SpectralTraceList)
    - MetisLMSSpectralTrace(SpectralTrace)
    - MetisLMSImageSlicer(ApertureMask)
    - MetisLMSEfficiency(TERCurve)
    - ExtraFitsKeywords(Effect)
    - EffectsMetaKeywords(ExtraFitsKeywords)
    - SourceDescriptionFitsKeywords(ExtraFitsKeywords)
    - SimulationConfigFitsKeywords(ExtraFitsKeywords)
    - SpectralTraceListWheel(Effect)
    - Bias(Effect)
    
- version : 0.4.1rc1
  date : 2022-03-25
  comment : Updates since METIS science team release
  changes :
  - New Effects:
    - TopHatFilterCurve   
    - TopHatFilterWheel   
    - SpanishVOFilterWheel
  - DetectorList x(y)_size columns now accept units of pixel and mm
  - warnings and errors now handled using python logging package
  - minor bug fixes

- version : 0.4.0
  date : 2022-03-03
  comment : Version released for the METIS science team
  changes :
  - release of new spectroscopy effect SpectalTraceList
  - moved individual spectral trace code to utils file
  - rewritten FovManager class
  - added make_cube, make_image, make_spectrum to FieldOfView class
  - removed fov_grid from Effects
  - added new detector array z_order section (900-class effects)
  - wavelength-dependent PSF in spectroscopic modes
  - proper handling of cube sources
  - headers for output files
  
"""
