import math
import numpy as np

from astropy.modeling.models import Moffat2D

from astropy.coordinates import Angle
from astropy import units as u
from matplotlib.path import Path

def calculate_FWHM(wavelength,airmass,config):
    median_seeing=config['median_seeing']
    median_seeing_wl=config['median_seeing_wl']
    L0=config['wavefront_outerscale']
    D=config['telescope_diameter']
    
    r0=0.1*median_seeing**(-1)*(wavelength/median_seeing_wl)**1.2*airmass**(-0.6)  
    F_kolb=1/(1+300*(D/L0))-1
    
    FWHM_atm=median_seeing*airmass**(0.6)*(wavelength/median_seeing_wl)**(-0.2)*np.sqrt(1+F_kolb*2.183*(r0/L0)**0.356)
    FWHM_dl=0.212*wavelength/D
    
    FWHM_total=np.sqrt(FWHM_atm**2+FWHM_dl**2)
    
    return FWHM_total

def make_moffat_PSFs(wavelengths,offsets,airmass,diameter,config,beta=2.5):
    scale=config['sim_scale']
    
    wavelengths,offsets = np.array(wavelengths),np.array(offsets)
    boundary=math.ceil(diameter/2/scale) #radius of aperture in pixels

    FWHMs = calculate_FWHM(wavelengths,airmass,config)

    x = np.arange(-boundary, boundary+1)
    y = np.arange(-boundary, boundary+1)
    x, y = np.meshgrid(x, y)
          
    PSFs=np.zeros((len(wavelengths),boundary*2+1,boundary*2+1))
    for count in range(0,len(wavelengths)):
        alpha=FWHMs[count]/scale/(2*np.sqrt(2**(1/beta)-1))
        moffat_total=(np.pi*alpha**2)/(beta-1)
        x_pos=offsets[count][0]/scale  
        y_pos=offsets[count][1]/scale

        Moffat=Moffat2D(1,x_pos,y_pos,alpha,beta)
        Moffat_data=Moffat(x,y)
        PSFs[count]=Moffat_data/moffat_total

    return PSFs

def parallatic_angle(HA,dec,lat):
    HA=Angle(HA*u.hour).rad
    q = np.arctan2(np.sin(HA),(np.cos(dec)*np.tan(lat)-np.sin(dec)*np.cos(HA)))
    return q

def diff_shift(wave, airmass, atm_ref_wav,config):
    T = config['temperature']
    HR = config['humidity']
    P = config['pressure']
    Lambda0 = atm_ref_wav
    wave = wave
    
    ZD_deg = airmass_to_zenith_dist(airmass)
    ZD = np.deg2rad(ZD_deg)

    # saturation pressure Ps (millibars)
    PS = -10474.0 + 116.43*T - 0.43284*T**2 + 0.00053840*T**3

    # water vapour pressure
    Pw = HR * PS
    # dry air pressure
    Pa = P - Pw

    #dry air density
    Da = (1 + Pa * (57.90*1.0e-8 - 0.0009325/T + 0.25844/T**2)) * Pa/T
    
    #water vapour density
    Dw = (1 + Pw * (1 + 3.7 * 1E-4 * Pw) * (- 2.37321 * 1E-3 + 2.23366/T
                                            - 710.792/T**2
                                            + 77514.1/T**3)) * Pw/T

    S0 = 1.0/Lambda0
    S = 1.0/wave

    N0_1 = (1.0E-8*((2371.34+683939.7/(130.0-S0**2)+4547.3/(38.9-S0**2))*Da
            + (6487.31+58.058*S0**2-0.71150*S0**4+0.08851*S0**6)*Dw))

    N_1 = 1.0E-8*((2371.34+683939.7/(130.0-S**2)+4547.3/(38.9-S**2))*Da
                  + (6487.31+58.058*S**2-0.71150*S**4+0.08851*S**6)*Dw)

    DR = np.tan(ZD)*(N0_1-N_1) * u.rad

    return (DR.to(u.arcsec)).value

def airmass_to_zenith_dist(airmass):
    return np.rad2deg(np.arccos(1/airmass))

def line(A,B):
    m=(A[1]-B[1])/(A[0]-B[0])
    c=A[1]-m*A[0] 
    return m,c

def make_aperture(band,major_axis,scale):
    scale=scale
    boundary=math.ceil(major_axis/2/scale) #radius of aperture in pixels

    if band[0:6]=="LR_VIS" or band[0:6]=="LR_NIR" or band[0:6]=="HR_NIR":
        sampling = major_axis/3/scale
        
        triangle_side=sampling*np.sqrt(3)/3
        aperture_centre=[boundary,boundary]
        
        centre_0=aperture_centre
        centre_1=[centre_0[0],centre_0[1]+sampling]
        centre_2=[centre_0[0]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6),centre_0[1]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6)]
        centre_3=[centre_0[0]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6),centre_0[1]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6)]
        centre_4=[centre_0[0],centre_0[1]-sampling]
        centre_5=[centre_0[0]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6),centre_0[1]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6)]
        centre_6=[centre_0[0]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6),centre_0[1]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6)]

        centres=[centre_0,centre_1,centre_2,centre_3,centre_4,centre_5,centre_6]

        apertures=np.zeros((7,boundary*2+1,boundary*2+1))
        for count,centre in enumerate(centres):
            P1=[centre[0]+triangle_side*np.cos(np.pi*1/3),centre[1]+triangle_side*np.sin(np.pi/3)]
            P2=[centre[0]+triangle_side,centre[1]]
            P3=[centre[0]+triangle_side*np.cos(np.pi/3),centre[1]-triangle_side*np.sin(np.pi/3)]
            P4=[centre[0]-triangle_side*np.cos(np.pi/3),centre[1]-triangle_side*np.sin(np.pi/3)]
            P5=[centre[0]-triangle_side,centre[1]]
            P6=[centre[0]-triangle_side*np.cos(np.pi/3),centre[1]+triangle_side*np.sin(np.pi/3)]

            polygon=[P1,P2,P3,P4,P5,P6]
            height=boundary*2+1
            width=boundary*2+1
            poly_path=Path(polygon)

            x, y = np.mgrid[:height, :width]
            coors=np.hstack((x.reshape(-1, 1), y.reshape(-1,1))) # coors.shape is (4000000,2)

            mask = poly_path.contains_points(coors)
            mask=mask.reshape(height, width)
            
            apertures[count]=mask
            
    return apertures
