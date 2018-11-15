"""
Include the Class structures that are relevant for an optical train

What data does an Optical Train need to contain? 


What common functionality do all optical trains have?




- 1D effects (spectral)
    - transmission per surface
    - reflection per surface
    - emission per surface
    
- 2D effects (spatial)
    - rotation offset
    - linear offset
    - distortion
    
    - convolution
    - rotation smear
    - linear smear
    
    - unequal illumination
    - increase in background
    
- 3D effects (spectro-spatial)
    - all of the above but for each of the wavelength slices

    
Spectral effects
----------------

As always, the purely spatial effects are taken into consideration first

- transmission
- emission

these should be input from a radiometric description from a single file for all
surfaces in the optical train. 
The radiometric treatment of each surface is independent of the position. 
If units aren't given the thermal emission is calculated with the temperature
area and etandue of the optic


SpectralSurface
lam, transmissivity, reflectivity, emissivity, temperature, units

- the main part of the SpectralSurface will be a series of pysynphot 
curves
- the self.__dict__ can take all the parameters
- functionality includes:
    - photons_in_range()
    - multiplication and addition
    - scaling to a certain magnitude or flux level
    
    
    


Spatial Effects
---------------
Can be broken up into the following categories:    
    
    
SpatialShift
lam, x, y, dx, dy

SpatialSpread
lam, x, y, dx, dy

SpatialConvolution
lam, x, y, kernel [FITS extension]

SpatialIntensity
lam, x, y, intensity

RotationalShift
lam, x, y, angle

RotationalSpread
lam, x, y, angle


contains a table of lam, x, y, angle
if only one entry then it is a 2D effect and is relevant to all slices
if the table has more that one entry, the effect is 3D and should be interpolated
for the required wavelength

a SpatialShift for example could be used just to move the fov, or to hold the
positions of a slit image on the FPA
a SpatialShift can also be used for both the atmospheric dispersion and a 
non-standard ADC correction
a SpatialShift can be used for implementing distortions if the 

a SpatialSpread can be used to simulate telescope tracking errors

PSFs are a 3D SpatialConvolution and also follow the lam, x, y, kernel format

NCPAs could also follow a 3D SpatialConvolution format, but for a single x and y

Flatfielding can be a instance of SpatialIntensity 

How to comine all the effects
-----------------------------
An optical train should work out for each slice
- how many combinations of effects are there, and
- the region borders on the focal plane where each combination is valid



Schedule of operations
----------------------
The spectral effects come first. Make a dictionary with all the spectral curves
- keep them ordered - ordered.dict
- run through and do the radiometry for background photons
- keep the filter curve and the atmospheric profile seperate
    - it's the question of whether we recompile the system transmission curve
    everytime we run the simulation, or just once when the optical train is made
    
A compromise would be to make the full optical train system curve, except for
the filters. 

All the backkground photons will arrive before the filter, thus a telescope plus
instument curve can be generated and kept. thus this is the middle curve
put a flag on the SpectralSurface object with "dynamic" or "static". Dynamic
curves are recomputed everytime the "observe" command is called. 
If they are static, they are then included in the system transmission and emission
This way we can also include the atmosphere as a dynamic part.



    
    
    

"""

def optics_optical_train():
    pass
    
    
    
    
    
    