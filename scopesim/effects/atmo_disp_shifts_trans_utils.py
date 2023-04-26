
import math
import numpy as np

from configobj import ConfigObj

from astropy.modeling.models import Moffat2D
from astropy.modeling.functional_models import Disk2D

import matplotlib.pyplot as plt

from astropy.coordinates import Angle
from astropy import units as u

def calculate_FWHM(wavelength,airmass):
    Config = ConfigObj('conf.ini')

    D = float(Config['Parameters']['telescope_diameter']) 
    L0 = float(Config['Parameters']['wavefront_outer_scale'])
    median_seeing = float(Config['Seeing']['median_seeing'])
    median_seeing_wl = float(Config['Seeing']['median_seeing_wl'])
       
    r0=0.1*median_seeing**(-1)*(wavelength/median_seeing_wl)**1.2*airmass**(-0.6)  
    F_kolb=1/(1+300*(D/L0))-1
    
    FWHM_atm=median_seeing*airmass**(0.6)*(wavelength/median_seeing_wl)**(-0.2)*np.sqrt(1+F_kolb*2.183*(r0/L0)**0.356)
    FWHM_dl=0.212*wavelength/D
    
    FWHM_total=np.sqrt(FWHM_atm**2+FWHM_dl**2)
    
    return FWHM_total

def make_moffat_PSFs(wavelengths,offsets,airmass,diameter,beta=2.5):
    
    wavelengths,offsets = np.array(wavelengths),np.array(offsets)
    Config = ConfigObj('conf.ini')
    scale = float(Config['Sim_Parameters']['scale'])
     
    boundary=math.ceil(diameter/2/scale) #radius of aperture in pixels

    FWHMs = calculate_FWHM(wavelengths,airmass)

    x = np.arange(-boundary, boundary+1)
    y = np.arange(-boundary, boundary+1)
    x, y = np.meshgrid(x, y)
          
    PSFs=np.zeros((len(wavelengths),boundary*2+1,boundary*2+1))
    total_vals=[]

    for count in range(0,len(wavelengths)):
        alpha=FWHMs[count]/scale/(2*np.sqrt(2**(1/beta)-1))
        x_pos=offsets[count][0]/scale  
        y_pos=offsets[count][1]/scale

        Moffat=Moffat2D(1,x_pos,y_pos,alpha,beta)
        Moffat_data=Moffat(x,y)
  
        PSFs[count]=Moffat_data
        
        moffat_total=(np.pi*alpha**2)/(beta-1)
        total_vals.append(moffat_total) 

    return PSFs,total_vals

def parallatic_angle(HA,dec,lat):
    HA=Angle(HA*u.hour).rad
    q = np.arctan2(np.sin(HA),(np.cos(dec)*np.tan(lat)-np.sin(dec)*np.cos(HA)))
    return q

def diff_shift(wave, airmass, atm_ref_wav):
    Lambda0 = atm_ref_wav
    wave = wave
    
    config=ConfigObj('conf.ini')

    T = float(config["Environment"]["temperature"])+273.15
    HR = float(config["Environment"]["humidity"])/100
    P = float(config["Environment"]["pressure"])
    
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

    #1 - P instead of 1 + Pa here? Why? Makes minimal affect of actual values...
    
    #water vapour density ?
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

def make_aperture(type,major_axis,hex_rotation=0):
    config=ConfigObj('conf.ini')
    scale=float(config['Sim_Parameters']['scale'])
    boundary=math.ceil(major_axis/2/scale) #radius of aperture in pixels
    if type == "circle":
    
        x = np.arange(-boundary, boundary+1)
        y = np.arange(-boundary, boundary+1)
        x, y = np.meshgrid(x, y)
 
        Disk=Disk2D(1,0,0,major_axis/2/scale)
        aperture=Disk(x,y)    

        return aperture
    
    if type == "hexagons":
        sampling = major_axis/3/scale
        aperture_array=np.zeros([boundary*2+1,boundary*2+1])

        triangle_side=sampling*np.sqrt(3)/3
        core = 2 * triangle_side * np.cos(np.pi/4)
        aperture_centre=[boundary,boundary]
        alpha = hex_rotation
        
        centre_0=aperture_centre
        centre_1=[centre_0[0]+sampling*np.cos(np.pi/2-alpha),centre_0[1]+sampling-sampling*(1-np.sin(np.pi/2-alpha))]
        centre_2=[centre_0[0]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6-alpha),centre_0[1]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6-alpha)]
        centre_3=[centre_0[0]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6+alpha),centre_0[1]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6+alpha)]
        centre_4=[centre_0[0]-sampling*np.cos(np.pi/2-alpha),centre_0[1]-sampling+sampling*(1-np.sin(np.pi/2-alpha))]
        centre_5=[centre_0[0]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6-alpha),centre_0[1]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6-alpha)]
        centre_6=[centre_0[0]-np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.cos(np.pi/6+alpha),centre_0[1]+np.sqrt((triangle_side*3/2)**2+(sampling/2)**2)*np.sin(np.pi/6+alpha)]

        centres=[centre_0,centre_1,centre_2,centre_3,centre_4,centre_5,centre_6]

        for centre in centres:
            if alpha == 0:
                P1=[centre[0]+triangle_side*np.cos(np.pi*1/3-alpha),centre[1]+triangle_side*np.sin(np.pi/3-alpha)]
                P2=[centre[0]+triangle_side*np.cos(alpha),centre[1]-triangle_side*np.sin(alpha)]
                P3=[centre[0]+triangle_side*np.cos(np.pi/3+alpha),centre[1]-triangle_side*np.sin(np.pi/3+alpha)]
                P4=[centre[0]-triangle_side*np.cos(np.pi/3-alpha),centre[1]-triangle_side*np.sin(np.pi/3-alpha)]
                P5=[centre[0]-triangle_side*np.cos(-alpha),centre[1]-triangle_side*np.sin(-alpha)]
                P6=[centre[0]-triangle_side*np.cos(np.pi/3+alpha),centre[1]+triangle_side*np.sin(np.pi/3+alpha)]

                L12_m,L12_c=line(P1,P2)
                L23_m,L23_c=line(P2,P3)
                L34_m,L34_c=line(P3,P4)
                L45_m,L45_c=line(P4,P5)
                L56_m,L56_c=line(P5,P6)
                L61_m,L61_c=line(P6,P1)
                           
                for y in range(0,len(aperture_array)):
                    for x in range(0,len(aperture_array)):           
                        if x < centre_0[0] + triangle_side * 2 and x > centre_0[0] - triangle_side * 2 and y < centre_0[1] + sampling and y > centre_0[1] - sampling: 
                            aperture_array[y][x]=1         
                        elif y < L61_m*x + L61_c and y > L34_m*x + L34_c and y < L12_m*x + L12_c and y > L23_m*x + L23_c and y > L45_m*x + L45_c and y < L56_m*x + L56_c:
                            aperture_array[y][x]=1  
                            
            elif alpha != np.pi/6 and alpha != -np.pi/6:
                P1=[centre[0]+triangle_side*np.cos(np.pi*1/3-alpha),centre[1]+triangle_side*np.sin(np.pi/3-alpha)]
                P2=[centre[0]+triangle_side*np.cos(alpha),centre[1]-triangle_side*np.sin(alpha)]
                P3=[centre[0]+triangle_side*np.cos(np.pi/3+alpha),centre[1]-triangle_side*np.sin(np.pi/3+alpha)]
                P4=[centre[0]-triangle_side*np.cos(np.pi/3-alpha),centre[1]-triangle_side*np.sin(np.pi/3-alpha)]
                P5=[centre[0]-triangle_side*np.cos(-alpha),centre[1]-triangle_side*np.sin(-alpha)]
                P6=[centre[0]-triangle_side*np.cos(np.pi/3+alpha),centre[1]+triangle_side*np.sin(np.pi/3+alpha)]

                L12_m,L12_c=line(P1,P2)
                L23_m,L23_c=line(P2,P3)
                L34_m,L34_c=line(P3,P4)
                L45_m,L45_c=line(P4,P5)
                L56_m,L56_c=line(P5,P6)
                L61_m,L61_c=line(P6,P1)
                
                for y in range(0,len(aperture_array)):
                    for x in range(0,len(aperture_array)):
                        if y > centre_0[1] - core and y < centre_0[1] + core and x > centre_0[1] - core and x < centre_0[1] + core:
                            aperture_array[y][x]=1
                        elif y < L61_m*x + L61_c and y > L34_m*x + L34_c and y < L12_m*x + L12_c and y > L23_m*x + L23_c and y > L45_m*x + L45_c and y < L56_m*x + L56_c:
                            aperture_array[y][x]=1
                                
            elif alpha == np.pi/6 or alpha == - np.pi/6:
                P1=[centre[0]+triangle_side*np.cos(np.pi/3-alpha),centre[1]+triangle_side*np.sin(np.pi/3-alpha)]
                P2=[centre[0]+triangle_side*np.cos(alpha),centre[1]-triangle_side*np.sin(alpha)]
                P3=[centre[0]+triangle_side*np.cos(np.pi/3+alpha),centre[1]-triangle_side*np.sin(np.pi/3+alpha)]
                P4=[centre[0]-triangle_side*np.cos(np.pi/3-alpha),centre[1]-triangle_side*np.sin(np.pi/3-alpha)]
                P5=[centre[0]-triangle_side*np.cos(-alpha),centre[1]-triangle_side*np.sin(-alpha)]
                P6=[centre[0]-triangle_side*np.cos(np.pi/3+alpha),centre[1]+triangle_side*np.sin(np.pi/3+alpha)]
                
                L23_m,L23_c=line(P2,P3)
                L34_m,L34_c=line(P3,P4)
                L56_m,L56_c=line(P5,P6)
                L61_m,L61_c=line(P6,P1)
                
                for y in range(0,len(aperture_array)):
                    for x in range(0,len(aperture_array)):
                        if y < centre_0[1] + triangle_side * 2 and y > centre_0[1] - triangle_side * 2 and x < centre_0[0] + sampling and x > centre_0[0] - sampling: 
                            aperture_array[y][x]=1    
                        elif y < L61_m*x + L61_c and y > L34_m*x + L34_c and y > L23_m*x + L23_c and  y < L56_m*x + L56_c and x > centre[0] - sampling/2 and x < centre[0] + sampling/2:
                            aperture_array[y][x]=1   
                             
        return aperture_array
    
def line(A,B):
    m=(A[1]-B[1])/(A[0]-B[0])
    c=A[1]-m*A[0] 
    return m,c
