version = '0.4.0'
date    = '2022-03-03 13:00:00 GMT'
yaml_descriptions = """
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
