z-orders
========

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
* 70 SpectralTraceList      !!! Write
* 80 ApertureMask
* 81 ApertureList           !!! Write

Detector
* 90 DetectorList


100..199 Make SystemThroughput (1D)
-----------------------------------
<OpticsManager>.surfaces_table

TERCurves
* 110 TERCurve              !!! Add scaling for mags
* 111 AtmosphericTERCurve
* 112 SkycalcTERCurve
* 113 QuantumEfficiencyCurve
* 114 FilterCurve

SurfaceLists
* 120 SurfaceList
* 121 MasterSurfaceList     !!! Write


200..299 Make FOVs (3D)
-----------------------
<OpticsManager>.fov_setup_effects

TERCurves
* 214 FilterCurve

SurfaceLists
* 221 MasterSurfaceList     !!! Write

Shift3Ds
* 231 AtmosphericDispersion
* 232 AtmosphericDispersionCorrection       !!! If quick_adc=True

Analytical PSFs
* 241 NonCommonPathAberration
* 242 GaussianDiffractionPSF
* 243 Seeing (PSF)
* 244 Vibration (PSF)

Semi Analytical PSFs
* 251 PoppyFieldVaryingPSF
* 252 PoppyFieldConstantPSF

Discrete PSFs
* 261 FieldVaryingPSF
* 262 FieldConstantPSF

Spectroscopic Trace maps
* 271 LongSlitTraceMap      !!! Write
* 272 IfuTraceMap           !!! Write
* 273 MosTraceMap           !!! Write

Apertures
* 280 ApertureMask
* 281 ApertureList          !!! Write
* 282 SquareApertureList    !!! Write
* 283 RoundApertureList     !!! Write
* 284 PolygonApertureList   !!! Write

Detectors
* 290 DetectorList


300..399 Make ImagePlane (2D)
-----------------------------
<OpticsManager>.image_plane_setup_effects

* 370 SpectralTraceList
* 380 ApertureMask
* 390 DetectorList


400..499 Make Detector (0D)
---------------------------
<OpticsManager>.detector_setup_effects

* 490 DetectorList


500..599 apply-to(Source)
-------------------------
<OpticsManager>.source_effects

* 521 MasterSurfaceList     !!! Write   (system throughput)


600..699 apply-to(FieldOfView)
------------------------------
<OpticsManager>.fov_effects

* 632 AtmosphericDispersionCorrection
* 640 PSF
    * in all variations
* 650 PupilPlaneEffect      !!! Write   (Integrated rotation etc)
* 651 IntegratedPupilRotation ! Write
* 652 NonSiderialTracking
* 670 SpectralTraceList
    * in all variations


700..799 apply-to(ImagePlane)
-----------------------------
<OpticsManager>.image_plane_effects

* 721 MasterSurfaceList     !!! Write  (bg emission)
* 744 Vibration
* 761 Vignetting            !!! Write
* 762 Distortion            !!! Write
* 763 Chopping              !!! Write
* 780 ReferencePixelBorder


800..899 apply-to(Detector)
---------------------------
<OpticsManager>.detector_effects

Noises
* 810 ReadNoise
* 811 BasicReadNoise
* 812 HawaiiReadNoise       !!! Write
* 813 AquariusReadNoise     !!! Write
* 820 ShotNoise

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