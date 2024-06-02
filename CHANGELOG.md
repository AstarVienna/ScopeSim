# Version 0.8.3
**2024-06-02**

Mostly a small hotfix to allow chaning of exptime in AutoExposure, plus some housekeeping.

## What's Changed
### Dependency Changes
* Bump requests from 2.31.0 to 2.32.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/417
* Update astropy version by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/418
### Other Changes
* Install ScopeSim_Data in the poetry environment by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/416
* Do not ignore PytestRemovedIn8Warning. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/420
* Fix configure_logging by using function scope. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/422
* Fix exptime dit ndit bug by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/424

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.8.2...v0.8.3

# Version 0.8.2
**2024-05-14**

## What's Changed
### New Features or Improvements
* Cap negative values below 0 before quantifying to an unsigned int. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/414

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.8.1...v0.8.2

# Version 0.8.1
**2024-05-13**

Small changes required for the first METIS Simulated data release.

## What's Changed
### New Features or Improvements
* Do something sensible when a trace falls outside the FoV. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/407
* Fix specref not always an integer by @JenniferKarr in https://github.com/AstarVienna/ScopeSim/pull/411
### Dependency Changes
* Bump tqdm from 4.66.1 to 4.66.3 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/409
* Bump jinja2 from 3.1.3 to 3.1.4 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/410
### Documentation Improvements
* Fix dev_master -> main in readme by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/406

## New Contributors
* @JenniferKarr made their first contribution in https://github.com/AstarVienna/ScopeSim/pull/411

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.8.0...v0.8.1

# Version 0.8.0
**2024-04-15**

Many small fixes, some new effects, some important fixes related to coordinates, and lots of cleanup.

## What's Changed

### API Changes
* Add more useful error message to download functions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/309
* Migrate from `requests` ➡️ `httpx` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/312
* Improve logging by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/339
* Use ChainMap for UserCommands by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/375

### Changes to or addition of Effects
* Add Shutter effect by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/304
* Add Quantization effect by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/308
* Implement apply decision for quantization by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/396

### New Features or Improvements
* Add basic sky coordinates (WCS) to ScopeSim output by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/307
* Include more progress bars by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/311
* Further improvements to logging by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/349
* Add IFU cube rectification and more by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/258
* Resolve recursive bang-strings by @astronomyk in https://github.com/AstarVienna/ScopeSim/pull/351
* Also show scopsim version in bug_report by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/394

### Dependency Changes
* Improve CI run for notebooks by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/300
* Migrate to Poetry by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/314
* Bump nbconvert from 6.4.5 to 6.5.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/315
* Bump jupyter-server from 1.13.5 to 2.11.2 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/316
* Bump requests from 2.28.2 to 2.31.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/317
* Upgrade numpy to 1.26.3 and some other dependencies by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/336
* Drop support for Python 3.8 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/327
* Bump jinja2 from 3.1.2 to 3.1.3 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/338
* Bump jupyter-lsp from 2.2.1 to 2.2.2 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/344
* Bump notebook from 7.0.6 to 7.0.7 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/347
* Bump jupyterlab from 4.0.10 to 4.0.11 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/346
* Bump pillow from 10.1.0 to 10.2.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/350
* Some small dependency- and version-related changes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/363
* Bump pillow from 10.2.0 to 10.3.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/393
* Bump idna from 3.6 to 3.7 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/398

### Documentation Improvements
* Add config file for auto-generated release notes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/301
* Move changelog to dedicated file, add more readme badges by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/302
* Use PyPI badge for Python versions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/326
* Replace "Telescopy" with "ScopeSim" in README by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/348
* Also make pdf and epub by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/370
* Fix RTD Poetry configuration by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/379

### Other Changes
* Some refactoring of z_order fuctionality by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/303
* Fix notebook tests by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/320
* Include ScopeSim_Data in notebook tests by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/324
* Allow runnotebooks.sh to run without arguments by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/330
* Revert "Also use poetry for calling jupytext" by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/334
* Use new linkchecker action by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/335
* Properly stack stars by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/337
* Remove obsolete files by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/340
* Rearrange CI tests by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/341
* Add test to see whether all Python files can be imported. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/343
* Ensure logging messages don't reach the root logger by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/345
* Sort corner pixels to deal with negative CDELTs by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/321
* Use logger instead of print by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/353
* Minor logging fixes in download module by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/354
* Add DeprecationWarnings for fov_grid methods by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/313
* Delete redundant vesion.py by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/355
* Minor formatting changes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/358
* Add more debug logging by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/356
* Further harmonize `filename` kwarg by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/361
* Minor cleanup in `user_commands.py` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/362
* Replace printing with logging by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/360
* Removing currsys as global parameter by @astronomyk in https://github.com/AstarVienna/ScopeSim/pull/364
* Additional hotfix for the removed currsys by @astronomyk in https://github.com/AstarVienna/ScopeSim/pull/368
* More cmds and kwargs stuff by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/369
* Remove unused and broken `make_imagehdu_from_table()` by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/371
* Fixes needed for IFU/LMS mode by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/376
* Do not set user commands as rc.__currsys__ by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/377
* Weed out unused utils functions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/381
* Refactor some rarely-used utils functions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/382
* Make `required_keys` always a `set` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/383
* Some minor improvements and refactoring in the FOVManager by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/384
* Slightly more sophisticated use of numpy by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/385
* Remove unsatisfied assert by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/386

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.7.1...v0.8.0

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
