z-orders
========
Each ``Effect``-object has one or more ``z`` values. This value tells
ScopeSim where to use the ``Effect`` when creating an ``OpticalTrain`` and
running a simulations.

There are 8 major sections during a ScopeSim simulation: 4 while setting up the
model and 4 while running the simulation. In brief they are

* Setup

  1. Generate a global system throughput
  2. Generate a list of monochromatic "field of view" slices
  3. Generate the image plane for the flux expectation image (ph/s/pix)
  4. Generate the detectors

* Simulation

  1. Apply effects to pure ``Source`` objects
  2. Apply effects to ``FieldOfView`` slices from the ``Source`` objects
  3. Apply effects to ``ImagePlane``
  4. Apply effects to `` Detector`` cut-outs from ``ImagePlane``



0..99 Base classes
------------------
Source

* 10 TERCurve               !!! Add scaling for mags
* 20 SurfaceList

FOVs

* 30 Shift3D
* 40 PSF
* 41 AnalyticalPSF
* 42 SemiAnalyticalPSF
* 43 DiscretePSF
* 50 PupilPlaneEffect       !!! Write   (Integrated rotation etc)

ImagePlane

* 60 FieldPlaneEffect       !!! Write   (Distortion, Vignetting)
* 70 SpectralTraceList
* 80 ApertureMask
* 81 ApertureList

Detector

* 90 DetectorList


100..199 Make SystemThroughput (1D)
-----------------------------------
``<OpticsManager>.surfaces_table``

TERCurves

* 110 TERCurve              !!! Add scaling for mags
* 111 AtmosphericTERCurve
* 112 SkycalcTERCurve
* 113 QuantumEfficiencyCurve
* 114 FilterCurve

SurfaceLists

* 120 SurfaceList
* 124 FilterWheel


200..299 Make FOVs (3D)
-----------------------
``<OpticsManager>.fov_setup_effects``

TERCurves

* 220 SurfaceList
* 224 FilterCurve

SurfaceLists

* 221 FilterWheel

Shift3Ds

* 231 AtmosphericDispersion
* 232 AtmosphericDispersionCorrection       !!! If quick_adc=True

Analytical PSFs

* 241 NonCommonPathAberration
* 242 GaussianDiffractionPSF
* 243 SeeingPSF (PSF)
* 244 Vibration (PSF)

Semi Analytical PSFs

* 251 PoppyFieldVaryingPSF
* 252 PoppyFieldConstantPSF

Discrete PSFs

* 261 FieldVaryingPSF
* 262 FieldConstantPSF

Spectroscopic Trace maps

* 270 SpectralTraceList

Apertures

* 280 ApertureMask
* 281 ApertureList

Detectors

* 290 DetectorList


300..399 Make ImagePlane (2D)
-----------------------------
``<OpticsManager>.image_plane_setup_effects``

* 380 ApertureMask
* 390 DetectorList


400..499 Make Detector (0D)
---------------------------
``<OpticsManager>.detector_setup_effects``

* 490 DetectorList


500..599 apply-to(Source)
-------------------------
``<OpticsManager>.source_effects``

TERCurves

* 510 TERCurve              !!! Add scaling for mags
* 511 AtmosphericTERCurve
* 512 SkycalcTERCurve
* 513 QuantumEfficiencyCurve
* 514 FilterCurve

SurfaceLists

* 520 SurfaceList
* 524 FilterWheel


600..699 apply-to(FieldOfView)
------------------------------
``<OpticsManager>.fov_effects``

* 632 AtmosphericDispersionCorrection
* 640 PSF
    * in all variations
* 650 PupilPlaneEffect      !!! Write   (Integrated rotation etc)
* 651 IntegratedPupilRotation ! Write
* 652 NonSiderialTracking


700..799 apply-to(ImagePlane)
-----------------------------
``<OpticsManager>.image_plane_effects``

* 721 MasterSurfaceList     !!! Write  (bg emission)
* 744 Vibration
* 761 Vignetting            !!! Write
* 762 Distortion            !!! Write
* 780 ReferencePixelBorder
* 790 AutoDitty


800..899 apply-to(Detector)
---------------------------
``<OpticsManager>.detector_effects``

Noises

* 810 ReadNoise
* 811 BasicReadNoise
* 812 HawaiiReadNoise       !!! Write
* 813 AquariusReadNoise     !!! Write
* 820 ShotNoise
* 863 Chopping              !!! Write

Extra flux

* 830 DarkCurrent

Other phenomena

* 840 LinearityCurve
* 841 PixelCrossTalk        !!! Write
* 842 PixelLeakage          !!! Write
* 850 BadPixelMask          !!! Write
* 851 GainMask              !!! Write
* 852 PedestalMask          !!! Write

Exposures

* 860 SummedExposure
* 870 BinnedImage