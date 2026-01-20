# Version 0.11.2
**2026-01-20**

Maintenance and bugfix release containing some important fixes regarding flux values.

## What's Changed
### Bugs fixed
* Clean up FieldVaryingPSF by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/833
* Rename lamp off state from 'none' to 'off' by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/835
* Remove numpy deprecation warning by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/841
* Fix CUNIT not respected in FVPSF effect by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/851
* Transmission curves to extrapolate with value 0 by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/859
* Extrapolate all transmissions with 0 by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/863
* Fix Astropy warning about HIERARCH keywords by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/864
* Fix WCU flux mismatch by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/860
### Changes to or addition of Effects
* PSFs associated to PupilMaskWheel, download from server by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/838
* `ReferencePixelBorder` newly implemented by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/840
* Interactive updating of SkycalcTERCurve by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/853
### Dependency Changes
* Bump the dev-dependencies group with 2 updates by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/831
* Bump pyyaml from 6.0.2 to 6.0.3 in the main-dependencies group by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/830
* Bump actions/checkout from 5 to 6 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/837
* Bump the dev-dependencies group with 2 updates by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/843
* Bump fonttools from 4.58.4 to 4.61.0 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/844
* Bump urllib3 from 2.5.0 to 2.6.0 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/849
* Bump synphot from 1.6.0 to 1.6.1 in the main-dependencies group by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/854
* Bump the dev-dependencies group with 2 updates by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/855
* Bump pytest-cov from 6.3.0 to 7.0.0 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/856
* Bump urllib3 from 2.6.0 to 2.6.3 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/858
* Sync dependencies with ecosystem by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/861
* Allow Python 3.14, sync dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/865
### Documentation Improvements
* Fix email domain for Kieran Leschinski by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/857
* Fix formatting and docstrings, add change messages by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/862
### Other Changes
* Uppercase HIERARCH keywords by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/866

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.11.1...v0.11.2


# Version 0.11.1
**2025-10-26**

## What's Changed
### API Changes
* Allow log level to be lower case by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/811
### Bugs fixed
* Rename `Effect.__getitem__()` ➡️ `.get_from_meta()` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/815
* Improve handling of BUNIT in FOV by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/822
* Fix BUNIT not added for Source from image + flux by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/823
* Catch extra special units by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/824
* No ABmag in BUNIT by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/826
* Use lower case photlam for headers by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/827
### Changes to or addition of Effects
* Inter-pixel capacitance effect by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/806
* Compute METIS/LMS predisperser angle for given central wavelength  by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/812
* Local TER library by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/765
### Dependency Changes
* Bump jupyterlab from 4.2.5 to 4.4.8 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/805
* Bump sphinx from 7.3.7 to 8.1.3 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/809
### Documentation Improvements
* Fix some docstring formatting, add version info by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/828
### Other Changes
* Refactor data collections download by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/810
* Replace `np.trapz` with `scipy.integrate.trapezoid` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/813
* No longer ignore some warnings by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/814

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.11.0...v0.11.1


# Version 0.11.0
**2025-09-21**

Large release in preparation of the ScopeSim Workshop happening next week. See below for important changes and deprecations.

## Highlights
* Finally add support for Python 3.13 and Numpy 2.
* Partially resolve circular dependency between ScopeSim and ScopeSim_Templates
* Several new effects in preparation for MOSAIC simulations (#796).
* New and improved plotting methods on some classes for easier and nicer visualisation of results.

## What's Changed
### Deprecations
#### Download functionality
* Deprecate downloading IRDB packages directly from GitHub by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/727

Installing IRDB packages directly from GitHub using `download_packages()` is deprecated.
We're not aware that this functionality saw any significant use by the community.
If you do need to install unreleased IRDB packages from the GitHub repository (and you know _why_ you need that!), the same result can be achieved by cloning the repository (or your own fork) and then use `sim.link_irdb(<path_to_clone>)` to point ScopeSim there.
That solution has the added benefit that you can keep thing in sync through git alone, which is more reliable anyway.
The relevant functionality will be removed in version 0.12.0, at which point attempting to use `download_packages()` with a GitHub ref will result in an error.

Passing a list as the first argument to `download_example_data()` has been deprecated since version 0.8.4 and will now result in a `TypeError` being raised.
This is to catch any remaining uses of the old call signature of this function.
From version 0.12.0 onwards, arguments will be silently forwarded to the file retriever, meaning if a list is still passed at that point, it will result in an error further downstream.

#### Source Fields
* Deprecate direct attribute access in `SourceField` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/733

After the `SourceField` refactor in version 0.9.0, it was still possible to access any item of the underlying data structur (`Table` or `ImageHDU`) directly through the `SourceField` via a standard dict-like syntax.
This was retained for backwards compatibility, because some code still relied on this functionality.
It does however have the potential to leade to some hard-to-track bugs and is generally an unclean solution, which is why we have been working behind the scenes to slowly remove all internal code that still uses this feature.
Now because this is related to how the `Source` object might be put together in user code, we're using a more carefuly approach in deprecating this, to avoid breaking existing code for no urgent reason.
We're thus activating a `PendingDeprecationWarning` now, which should give users an advanced warning that we will probably fully deprecate this in the future.
If you do encounter this warning, please try to use the `.field` attribute to explicitly access the underlying data structur before using dict-like lookup on it.

### API Changes
* Improve handling of DeprecationWarnings by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/730
* Remove obsolete code by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/737
* Accept filename and relative path as well as absolute path by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/758
### Bugs fixed
* Fix creation of cube source from file name by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/559
* Fix WCU radiometry by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/754
* Use `RegularGridInterpolator` for proper extrapolation of small PSF images by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/798
### Changes to or addition of Effects
* Functionality for MOSAIC minimum viable product by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/796
### New Features or Improvements
* Add plotter utils functions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/734
* Add `.plot()` method to `ImagePlane` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/735
* Improve cube plotter by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/736
### Dependency Changes
* Bump synphot from 1.5.0 to 1.6.0 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/713
* Minor bumps on a few dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/718
* Bump actions/checkout from 4 to 5 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/756
* Update internal dependencies more often by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/769
* Bump pyyaml from 6.0.1 to 6.0.2 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/771
* Bump lxml from 5.2.2 to 5.4.0 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/772
* Refine dependabot configuration by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/780
* Bump the dev-dependencies group with 8 updates by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/782
* Update to Poetry 2 based pyproject.toml, update workflows by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/783
* Bump actions/setup-python from 5 to 6 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/786
* Bump packaging from 23.2 to 25.0 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/794
* Bump astar-utils from 0.3.3 to 0.3.4 by @dependabot[bot] in https://github.com/AstarVienna/ScopeSim/pull/791
* Remove semicircular dependency on ScopeSim_Templates by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/799
* Allow numpy 2.x by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/784
* Add support for Python 3.13 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/800
### Documentation Improvements
* Remove mention of deprecated GitHub download from docs by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/745
* Fix toml classifiers and URLs by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/795
* Docstrings, formatting, linting by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/802
### Other Changes
* Ignore a random DeprecationWarning from Pillow by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/715
* Various typing fixes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/728
* Fix locally failing mock path test by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/729
* Refactor `flatten_dict()` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/725
* Remove empty `spectrograph.py` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/738
* Split `FieldOfView` into subclasses by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/739
* Remove obsolete `fov_manager_utils.py` and tests by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/740
* Remove obsolete `PoorMansHeader` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/741
* Remove some unused code by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/744
* Remove unused spectroscopy integration tests by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/746
* Remove `Effect.fov_grid()` in base class by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/747
* Fix CI comments by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/748
* Add basic pylint config in pyproject.toml by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/760
* Update workflows by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/801

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.10.0...v0.11.0


# Version 0.10.0
**2025-06-27**

Large release before the summer, with many improvements and also breaking changes, please make sure to also update your IRDB packages accordingly!

## Highlights
* Many new features, including the "simplified IFU cube readout" mode.
* More physical treatment of the "electrons to ADU" process, including different gain values (#604).
* Extensive refactoring and many many bug fixes.
* More consistency across different classes with similar attributes.
* Improved treatment of flux units (still WIP, see #662).
* Many improvements to visual outputs (especially in notebooks), a few new progress bars in long-running loops.

## What's Changed
### API Changes
* Refactor FOV, add `BackgroundSourceField` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/610
* Deprecate `FieldOfView.spectra` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/653
* Deprecate adding a `Table` directly to an `ImagePlane` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/657
* Implement simplified IFU mode with cube readout by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/654
* Rename `CubeSourceField.wave` to `.waveset` for consistency by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/666
* Rename `SummedExposure` ➡️ `ExposureIntegration` by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/701
* Change  `.wave` to `.waveset` by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/707
### Bugs fixed
* Fix BUNIT in `FieldOfView` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/659
* Fix source spectral resolution by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/669
* Fix METIS LSF axes orientation by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/673
* Remove "per pixel" from BUNIT for consistency by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/684
* Fix search path when updating inst pkgs path by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/710
### Changes to or addition of Effects
* ADConversion effect by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/455
* Add pretty repr to `Effect` base class to avoid cmds spam by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/661
* Allow PSF interpolation to be switched off by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/671
* Implement line spread function for IFU cube mode by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/670
* Convert to signed float before chopnodcombine by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/687
### New Features or Improvements
* Add `zeros_from_header()` utils function by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/656
* Add `._repr_html_()` magic method to `OpticalTrain` for nice notebook display by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/667
* Add version requirement support for IRDB packages by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/660
* Use progress bar for spectral traces by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/674
* More logging and another progress bar by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/697
### Dependency Changes
* Bump actions/download-artifact from 4.2.1 to 4.3.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/647
* Bump numpy, astropy, httpx and synphot by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/652
* Bump setuptools from 75.3.0 to 78.1.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/663
* Bump tornado to get rid of security warning by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/664
* Bump versions of internal dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/672
* Bring dependency versions back in sync with lock file by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/675
* Bump a few minor dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/681
* Bump requests from 2.32.2 to 2.32.4 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/690
* Bump spextra from 0.41.3 to 0.41.4 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/695
* Bump urllib3 from 1.26.19 to 2.5.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/703
* Bump astar-utils from 0.3.2 to 0.3.3 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/704
* Bump internal dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/711
### Documentation Improvements
* Add proper error message for removed `Quantization` effect by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/676
* Add version info in various docstrings by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/678
* Minor cleanup of docs by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/677
* Add units flow page to docs by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/689
### Other Changes
* Make better use of `SourceField` in `FieldOfView` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/651
* Remove obsolete `FovVolumeList.detector_limits` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/655
* Add replaceholder workflow by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/688
* Remove unused private attributes from `OpticsManager` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/694
* Refactor and cleanup in `OpticalElement` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/700

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.9.3...v0.10.0


# Version 0.9.3
**2025-04-25**

Addition of average exposure effect by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/603, otherwise mostly bugfixes, small quality-of-life improvements and some internal refactoring.

## What's Changed
### Bugs fixed
* Fix broken flux scaling for extended sources in spectroscopy by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/638
### Changes to or addition of Effects
* Implement average exposure effect by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/603
* Remove spurious base classes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/621
### New Features or Improvements
* Make setting the `inst_pkgs` path more convenient by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/627
* Check if all ref keys are found in spectra by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/633
* Ignore non-package directories for bug_report by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/631
* Switch off spammy INFO logs from httpx by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/642
### Dependency Changes
* Bump beautifulsoup to 4.13.3 to get rid of deprecation warning by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/600
* Disable GitHub-based webtests by default by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/601
* Bump actions/download-artifact from 4.1.9 to 4.2.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/602
* Bump python-dateutil from 2.8.2 to 2.9.0.post0 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/628
* Bump h11 from 0.14.0 to 0.16.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/645
### Other Changes
* Remove base class check in `Effect.apply_to()` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/620
* Fixup adding mock path to UserCommands by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/640
* Fixup adding mock path to UserCommands (again) by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/641
* Use `.to_value(<unit>)` instead of `.to(<unit>).value` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/643

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.9.2...v0.9.3


# Version 0.9.2
**2025-03-13**

Mainly the inclusion of METIS WCU effects in https://github.com/AstarVienna/ScopeSim/pull/494 plus a few bugfixes and quality-of-life improvements.

## What's Changed
### Bugs fixed
* Fix detectorwindow by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/507
* Don't fail subpixels for edge sources by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/549
* Ensure at least one pixel in header by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/586
### Changes to or addition of Effects
* METIS WCU effects by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/494
### Dependency Changes
* Bump abatilo/actions-poetry from 2 to 4 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/520
* Bump actions/download-artifact from 4.1.7 to 4.1.8 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/521
* Bump jinja2 from 3.1.4 to 3.1.5 by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/519
* Bump sphinx from 5.3.0 to 7.3.7 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/525
* Bump jupyter-sphinx from 0.2.3 to 0.5.3 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/528
* Bump sphinxcontrib-apidoc from 0.4.0 to 0.5.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/522
* Bump nbsphinx from 0.9.3 to 0.9.6 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/526
* Bump skycalc-ipy from 0.5.1 to 0.5.2 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/523
* Bump jupyter from 1.0.0 to 1.1.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/530
* Bump ipympl from 0.9.4 to 0.9.6 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/538
* Bump jupytext from 1.16.0 to 1.16.6 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/541
* Bump myst-nb from 1.1.2 to 1.2.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/558
* Bump matplotlib from 3.8.2 to 3.10.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/577
* Bump actions/download-artifact from 4.1.8 to 4.1.9 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/576
* Bump jinja2 from 3.1.5 to 3.1.6 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/582
### Documentation Improvements
* Add bugfix label to auto release notes by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/562
### Other Changes
* make test_applies_dark_current_with_level_of_dit deterministic by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/517
* Add dependabot configuration by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/518
* Fix runnotebooks.sh by sorting the notebooks by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/535
* Add a few more `@pytest.mark.slow` markers by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/560
* Generatorify `stringify_dict` function by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/563
* Update dependabot.yml by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/578
* fix fp mask emission by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/583
* Fine-tune dependabot config by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/591

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.9.1...v0.9.2


# Version 0.9.1
**2024-11-26**

Many small but important fixes and improvements that were missing in v0.9.0.

## What's Changed
### Changes to or addition of Effects
* Improve z_order by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/500
* Fix FieldVaryingPSF by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/506
### New Features or Improvements
* Enable DeprecationWarnings by default by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/508
* Minor fixes for ScopeSimple by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/510
* Add support for package and mode status keywords by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/509
### Dependency Changes
* Require Templates 0.6.0 and sync from there by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/501
* Bump tornado from 6.4.1 to 6.4.2 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/511
### Documentation Improvements
* Cleanup example notebooks by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/496
### Other Changes
* Added source scaling test by @janusbrink in https://github.com/AstarVienna/ScopeSim/pull/495
* Move logging config to separate yaml, configure individual loggers by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/502
* Ignore strip_cdata warning by lxml, introduced by bs4 by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/504
* Fix #490: make rescale_imagehdu more robust against dimension mismatches by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/503
* Make random test deterministic by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/512
* Fix two more issues with ScopeSimple by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/514

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.9.0...v0.9.1


# Version 0.9.0
**2024-11-10**

> [!IMPORTANT]
> The minimum required Python version for this package is now **3.10** (see Dependency Changes).
> Python 3.13 is currently not supported as it causes issues on some platforms that are not yet fully understood. We are currently working on fixing 3.13 support.


## What's Changed
### API Changes
* ScopeSimple by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/426
* Make `DataContainer` an attribute of `Effect` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/482
### New Features or Improvements
* Replace NaNs in images, log warning by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/466
* Make observe work with no source (empty field) by @oczoske in https://github.com/AstarVienna/ScopeSim/pull/483
* Add a new top level to CMDS nested mapping to store current simulation run settings by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/493
### Dependency Changes
* Bump actions/download-artifact from 3 to 4.1.7 in /.github/workflows by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/465
* Drop support for Python 3.9 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/471
* Bump some dependency versions by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/472
* Use sphinx-book-theme for RTD by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/473
* Bump two dependencies, use tqdm.auto by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/477
* Bump scipy from 1.11.4 to 1.14.1 and httpx from 0.23.0 to 0.26.0 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/480
* Remove skycalc_cli from dependencies by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/481
* Bump notebook from 7.0.7 to 7.2.2 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/486
* Limit supported Python version to 3.10 <= x < 3.13 by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/497
### Documentation Improvements
* Update source_from_images.ipynb by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/469
* Update example notebook to properly show the same observation by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/476
### Other Changes
* Add `SourceField` and subclasses, rework `Source` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/405
* Some small improvement, mostly focused on numpy by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/467
* Don't fail on symlinks for packages path by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/474
* Cleanup `DataContainer` in preparation for Effect refactoring by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/478
* Fix Poisson NaN bug and clarify by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/484
* Source object data scaling fix by @janusbrink in https://github.com/AstarVienna/ScopeSim/pull/485
* Workaround for #491 metadata.requires can return None by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/492
* Fixed cdelt calculation when scaling an imageHDU by @janusbrink in https://github.com/AstarVienna/ScopeSim/pull/490

## New Contributors
* @janusbrink made their first contribution in https://github.com/AstarVienna/ScopeSim/pull/485

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.8.4...v0.9.0


# Version 0.8.4
**2024-08-29**

*Last version to support Python 3.9*

This includes many small and not-so-small fixes and improvements, which are made available here to Python 3.9, while a few more major changes will be released soon as Version 0.9.0, but without support for Python 3.9.

## What's Changed
### Changes to or addition of Effects
* Test all cases of DIT & NDIT vs. exptime and AutoExposure by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/428
* Avoid in-place operations because of dtype conflicts by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/432
* Split PSF and Electronic effects into subpackages by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/434
* Improve SVO filter handling by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/454
* Added function to round the edges of square PSF kernels by @astronomyk in https://github.com/AstarVienna/ScopeSim/pull/459
### Dependency Changes
* Bump tornado from 6.4 to 6.4.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/427
* Bump urllib3 from 1.26.18 to 1.26.19 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/429
* Bump certifi from 2023.11.17 to 2024.7.4 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/435
* Bump zipp from 3.17.0 to 3.19.1 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/436
* Bump setuptools from 69.0.2 to 70.0.0 by @dependabot in https://github.com/AstarVienna/ScopeSim/pull/437
* Bump some dependencies to get rid of some indirect ones by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/445
* Use `pooch` for `download_example_data` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/446
* Bump matplotlib, synphot, pyerfa by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/450
### Other Changes
* Formatting and Refactoring of `nghxrg.py` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/431
* Rename `DetectorArray` ➡️ `DetectorManager` plus Docstrings and Refactoring by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/423
* Enlarge initial Field of View to suppord wide-field imagers by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/433
* Fix some exposure things by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/440
* Deepcopy readout, to prevent it being overwritten. Closes #439 by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/441
* Allow star to be shifted by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/443
* Attempt to catch crypto warning by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/444
* Use new main branch for DevOps workflows by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/448
* Refactor `example_data_utils` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/447
* Remove unused function `return_latest_github_actions_jobs_status()` by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/449
* Remove remaining master references by @teutoburg in https://github.com/AstarVienna/ScopeSim/pull/460
* Add units to DataContainer table directly. by @hugobuddel in https://github.com/AstarVienna/ScopeSim/pull/461

**Full Changelog**: https://github.com/AstarVienna/ScopeSim/compare/v0.8.3...v0.8.4

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
