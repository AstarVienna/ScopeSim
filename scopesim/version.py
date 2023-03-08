version = '0.5.5'
date    = '2023-03-08 18:00:00 GMT'
yaml_descriptions = """
- version : 0.5.5
  date : 2023-03-08
  comment : Hotfix for header keyword generators
  changes :
  - Return to § for incremental extension keywords #168
  - thin slit confusing dispersion direction #169
  - Adds unequal (i.e. 2x1) binning and option to rotate the CCD by integer multiples of 90 degrees #170
  - add filters and slits to wheels #176
  - psf_utils.rescale_kernel: fix for negative shifts #177
  - Fix bug where the ._meta_dicts can become longer than the .fields #178
  - Add test that Source() is additive identity #179
  - Allow astropy Units to be values in FITS headers. #180

- version : 0.5.4
  date : 2022-10-06
  comment : Hotfix for header keyword generators
  github_pr_url : https://github.com/AstarVienna/ScopeSim/pull/166
  changes :
  - incremental special characters for header keywords changed from `§`to `++`
  - source object function calls are now given their own FITS header keyword FNSRCn (function-call source N) due to astropy not liking the combination of HIERARCH and CONTINUE keywords

- version : 0.5.3
  date : 2022-09-29
  comment : Minor upgrade to Spec modes and to FITS keywords
  github_pr_url : https://github.com/AstarVienna/ScopeSim/pull/165
  changes :
  - Effect object ExtraFitsKeywords now has the ability to add keywords with incrementing index numbers based on the extension number
  - FOV + FOVManager + FOVVolumes classes now accept aperture_id as an argument
  - ApertureList effects object now has an apply_to function which splits the FOVVolumeList accordingly

- version : 0.5.2
  date : 2022-08-25
  comment : Update of DLC server URL to scopesim.univie.ac.at
  changes :
  - Updated MANIFEST.in to include all the files needed by the basic_instrument test optical train
  - Small update to allow iterative extension specific FITS header keywords. E.g. EXTNAME = DETn.DATA

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
