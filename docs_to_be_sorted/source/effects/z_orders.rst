z-orders
========

0..99 - Make a FOV list
-----------------------
TERCurves
* 10 TERCurve
* 11 AtmosphericTERCurve
* 12 SkycalcTERCurve
* 13 QuantumEfficiencyCurve
* 14 FilterCurve

SurfaceLists
* 20 SurfaceList

Shift3Ds

* 30 Shift3D
* 31 AtmosphericDispersion
* 32 AtmosphericDispersionCorrection

Analytical PSFs

* 40 AnalyticalPSF
* 41 NonCommonPathAberration
* 42 GaussianDiffractionPSF
* 43 Seeing (PSF)
* 44 Vibration (PSF)

Semi Analytical PSFs

* 50 SemiAnalyticalPSF
* 51 PoppyFieldVaryingPSF
* 52 PoppyFieldConstantPSF

Discrete PSFs

* 60 DiscretePSF
* 61 FieldVaryingPSF
* 62 FieldConstantPSF

Apertures

* 80 ApertureMask
* 81 ApertureList
* 82 SquareApertureList
* 83 RoundApertureList
* 84 PolygonApertureList

Misc

* 70 SpectralTraceList
* 90 DetectorList

100..199 - Make an image plane
------------------------------
* 70 SpectralTraceList
* 110 ApertureMask
* 111 ApertureList
* 112 SquareApertureList
* 113 RoundApertureList
* 114 PolygonApertureList
* 120 DetectorList

200..299 - Source altering effects
----------------------------------
* 210 TERCurve
* 211 AtmosphericTERCurve
* 212 SkycalcTERCurve
* 213 QuantumEfficiencyCurve
* 220 SurfaceList


300..399 - (3D) FOV specific effects
------------------------------------

Shift3D
* 330 Shift3D
* 331 AtmosphericDispersion
* 332 AtmosphericDispersionCorrection

Analytical PSFs
* 340 AnalyticalPSF
* 341 NonCommonPathAberration
* 342 GaussianDiffractionPSF
* 343 Seeing (PSF)

Semi Analytical PSFs
* 350 SemiAnalyticalPSF
* 351 PoppyFieldVaryingPSF
* 352 PoppyFieldConstantPSF

Discrete PSFs
* 360 DiscretePSF
* 361 FieldVaryingPSF
* 362 FieldConstantPSF


400..499 - (2D) FOV-independent effects
---------------------------------------
* 444 VibrationPSF

500..599 - Electronic effects
-----------------------------
* 500 DetectorList

* 510 ReadNoise
* 511 RandomReadNoise
* 512 HawaiiReadNoise
* 513 AquariusReadNoise

* 520 ShotNoise

* 530 DarkCurrent

* 540 LinearityCurve
* 541 PixelCrossTalk
* 542 PixelLeakage

* 550 BadPixelMask
* 551 GainMask
* 552 PedestalMask




