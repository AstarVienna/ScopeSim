# Version 0.7.1
**2023-11-07**

Bug fixes and better error report

## What's Changed
- Improve the output of `scopsim.bug_report()` and automatically log that report when an unhandled error occurs: https://github.com/AstarVienna/ScopeSim/pull/287
- Some bug fixes related to that same bug_report: https://github.com/AstarVienna/ScopeSim/pull/290, https://github.com/AstarVienna/ScopeSim/pull/291
- Improve ScopeSim's README file: https://github.com/AstarVienna/ScopeSim/pull/294
- Deal with warnings from the latest Python version 3.12: https://github.com/AstarVienna/ScopeSim/pull/289
- Internally restructure and clean the test suite, make sure individual tests are not influencing each other: https://github.com/AstarVienna/ScopeSim/pull/285

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.7.0...v0.7.1

# Version 0.7.0
**2023-10-18**

Off-by-one fix.

## What's Changed
- Fix a long-standing bug regarding the internal implementation of WCS coordinates, which had multiple consequences, see https://github.com/AstarVienna/ScopeSim/pull/276 for details.
- This fix might break some existing codes using work-arounds for the bug described above.

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.6.2...v0.7.0

# Version 0.6.2
**2023-09-14**

Patch with bugfixes and code improvements

## What's Changed
- Fix documentation on readthedocs.io: https://github.com/AstarVienna/ScopeSim/pull/269
- Bug fixes related to plotting methods: https://github.com/AstarVienna/ScopeSim/pull/270
- Formatting of code and documentation to meet community standards: https://github.com/AstarVienna/ScopeSim/pull/271
- General refactoring, some in preparation of more substantial bug fixes: https://github.com/AstarVienna/ScopeSim/pull/272

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.6.1...v0.6.2

# Version 0.6.1
**2023-09-06**

Patch with visualisation improvements and general refactoring

## What's Changed
- Improvements to console representation and plot methods of various classes: https://github.com/AstarVienna/ScopeSim/pull/252 , https://github.com/AstarVienna/ScopeSim/pull/260 , https://github.com/AstarVienna/ScopeSim/pull/263 , https://github.com/AstarVienna/ScopeSim/pull/266
- Refactoring and changes to the inheritance of some classes: https://github.com/AstarVienna/ScopeSim/pull/261 , https://github.com/AstarVienna/ScopeSim/pull/264 , https://github.com/AstarVienna/ScopeSim/pull/265
- Changes to the CI configuration to reduce fails caused by web request timeouts https://github.com/AstarVienna/ScopeSim/pull/255 , https://github.com/AstarVienna/ScopeSim/pull/262

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.6.0...v0.6.1

# Version 0.6.0
**2023-07-10**

Summer 2023

## What's Changed
- Rename MAORY to MORFEO https://github.com/AstarVienna/ScopeSim/pull/195
- Fix NCPA and PSF affecting spectroscopy https://github.com/AstarVienna/ScopeSim/pull/238
- Fix line widths bug https://github.com/AstarVienna/ScopeSim/pull/213
- Add rectification utilities https://github.com/AstarVienna/ScopeSim/pull/237
- Include grating efficiencies https://github.com/AstarVienna/ScopeSim/pull/215
- Improve downloading of IRDB https://github.com/AstarVienna/ScopeSim/pull/234
- Improve Windows support

## New Contributors
* @teutoburg made their first contribution in https://github.com/AstarVienna/ScopeSim/pull/216

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.5.6...v0.6.0

# Version 0.5.6
**2023-03-13**

Hotfix to include minimal set of SVO data

## What's Changed
- Run notebooks in CI https://github.com/AstarVienna/ScopeSim/pull/183
- Add SVO data because SVO is down https://github.com/AstarVienna/ScopeSim/pull/185
- Fix OpticalTrain shared cmds attribute and fix docstring https://github.com/AstarVienna/ScopeSim/pull/186

# Version 0.5.5
**2023-03-08**

v0.5.5 is the first release by the 2023 A*Vienna team

## What's Changed
- Return to § for incremental extension keywords https://github.com/AstarVienna/ScopeSim/pull/168
- thin slit confusing dispersion direction https://github.com/AstarVienna/ScopeSim/pull/169
- Adds unequal (i.e. 2x1) binning and option to rotate the CCD by integer multiples of 90 degrees https://github.com/AstarVienna/ScopeSim/pull/170
- add filters and slits to wheels https://github.com/AstarVienna/ScopeSim/pull/176
- psf_utils.rescale_kernel: fix for negative shifts https://github.com/AstarVienna/ScopeSim/pull/177
- Fix bug where the ._meta_dicts can become longer than the .fields https://github.com/AstarVienna/ScopeSim/pull/178
- Add test that Source() is additive identity https://github.com/AstarVienna/ScopeSim/pull/179
- Allow astropy Units to be values in FITS headers. https://github.com/AstarVienna/ScopeSim/pull/180

# Version 0.5.4
**2022-10-06**

Hotfix for header keyword generators

## What's Changed
- incremental special characters for header keywords changed from `§`to `++`
- source object function calls are now given their own FITS header keyword FNSRCn (function-call source N) due to astropy not liking the combination of HIERARCH and CONTINUE keywords

# Version 0.5.3
**2022-09-29**

Minor upgrade to Spec modes and to FITS keywords

## What's Changed
- Effect object ExtraFitsKeywords now has the ability to add keywords with incrementing index numbers based on the extension number
- FOV + FOVManager + FOVVolumes classes now accept aperture_id as an argument
- ApertureList effects object now has an apply_to function which splits the FOVVolumeList accordingly

# Version 0.5.2
**2022-08-25**

Update of DLC server URL to scopesim.univie.ac.at

## What's Changed
- Updated MANIFEST.in to include all the files needed by the basic_instrument test optical train
- Small update to allow iterative extension specific FITS header keywords. E.g. EXTNAME = DETn.DATA

# Version 0.5.1
**2022-07-12**

Update of DLC server URL to scopesim.univie.ac.at

## What's Changed
  - Changed URL in defaults.yaml file

# Version 0.5.0
**2022-04-22**

IFU Spectroscopy mode for METIS

## What's Changed
- The IFU effects for the METIS LMS mode
- Effects for including extra FITS keywords for the MICADO-ESO delivery
- Minor change to the OpticalTrain.readout method to allow custom FITS keywords to be added
- #-strings for accessing the .meta dict contents of Effect objects
- added the tests.mocks.basic_instrument package for test purposes
- refactored package downloading functions in server.database
- Packages can now be downloaded directly from a git commit
- new RTDs structure for docs based on ipynb files
- change to SkyCalcTERCurve to use local files to avoid calling the skycalc server

### New Effects:
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

# Version 0.4.1rc1
**2022-03-25**

Updates since METIS science team release

## What's Changed
- DetectorList x(y)_size columns now accept units of pixel and mm
- warnings and errors now handled using python logging package
- minor bug fixes

### New Effects:
- TopHatFilterCurve
- TopHatFilterWheel
- SpanishVOFilterWheel

# Version 0.4.0
**2022-03-03**

Version released for the METIS science team

## What's Changed
- release of new spectroscopy effect SpectalTraceList
- moved individual spectral trace code to utils file
- rewritten FovManager class
- added make_cube, make_image, make_spectrum to FieldOfView class
- removed fov_grid from Effects
- added new detector array z_order section (900-class effects)
- wavelength-dependent PSF in spectroscopic modes
- proper handling of cube sources
- headers for output files
