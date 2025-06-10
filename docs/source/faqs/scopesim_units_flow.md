# Approximate sketch of "units flow" in ScopeSim

Effects that change units are shown in hexagonal boxes.
Rectangular boxes are basically `ImageHDU` wrappers.
Parallelogram boxes are lists of other entities.

```{mermaid}
%%{init: {"theme": "dark"} }%%
flowchart TB
    Source(["Source [PHOTLAM(/arcsec2)]"])
    FOV1[/"FOV.fields"/]
    FOV2["FOV.hdu [ph/s]"]
    FOV3[/"FOV.fields"/]
    FOV4["FOV.hdu [ph/s/um/arcsec2]"]
    FOV5["FOV.hdu [ph/s]"]
    SPT{{"`SpectralTrace:
          dispersion, sum(sky)`"}}
    IMP["ImagePlane [ph/s]"]
    Det1["Detector [ph/s]"]
    QECurve{{"QECurve: ph/s -> e-/s"}}
    Det2["Detector [e-/s]"]
    SE{{"ExposureIntegration: sum(time)"}}
    Det3["Detector [e-]"]
    ADC{{"ADConversion: e- -> ADU"}}
    Det4["Detector [ADU]"]
    Output(["Output [ADU]"])

    Source-- extract -->FOV1
    subgraph Imaging
        FOV1-- "sum(wave), *area" -->FOV2
    end
    Source-- extract -->FOV3
    subgraph Spectroscopy
        FOV3-- *area -->FOV4
        FOV4-->SPT
        SPT-->FOV5
    end
    FOV2-- project -->IMP
    FOV5-- project -->IMP
    IMP-- "extract" -->Det1
    Det1-->QECurve
    QECurve-->Det2
    Det2-->SE
    SE-->Det3
    Det3-->ADC
    ADC-->Det4
    Det4-->Output
```
