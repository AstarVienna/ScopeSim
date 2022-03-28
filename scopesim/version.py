version = '0.4.1rc1'
date    = '2022-03-25 13:00:00 GMT'
yaml_descriptions = """
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
