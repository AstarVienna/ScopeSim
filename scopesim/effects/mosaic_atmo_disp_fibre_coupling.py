import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import matplotlib as mpl

from . import mosaic_atmo_disp_fibre_coupling_utils as disp_utils
from . import TERCurve


class MosaicAtmosDispFibreCoupling(TERCurve):
    def __init__(self, **kwargs):
        params = {"HA_start": 0,
                  "HA_end": 0,
                  "declination": 0,
                  "band": "HR_NIR_H",
                  "sampling": 0.01,     # um
                  "guide_waveref": -1,
                  "aperture_waveref": -1}

        super(TERCurve, self).__init__()
        self.meta.update(params)
        self.meta.update(kwargs)

        self.coupling_fractions = AtmosDispFibreCoupling(**kwargs)

    def apply_to(self, obj, **kwargs):

        if isinstance(obj, SourceBase):
            # return objects are: 1D wavelength, 2D transmission curves
            waves, trans_arr = self.coupling_fractions.run()


        return base


class AtmosDispFibreCoupling:
    """
    Class to quantify atmospheric dispersion effects on a MOSAIC integration.
    See .run for inputs and outputs
    
    Author: Jay Stephan
    
    Bugs/Issues:
    HR Wavelength and aperture major axis values need to be confirmed
    """
    def __init__(self, **kwargs):
        #Simulation parameters for MOSAIC
        params = {'telescope_diameter':39, #m, diameter of telescope (ELT)
                  'wavefront_outerscale':46, #m, wavefront outer scale for FWHM change with airmass/wavelength. Value to be confirmed.
                  'median_seeing':.68, #arcsec, median seeing at Paranal
                  'median_seeing_wl':.5, #um, wavelength median seeing corresponds to

                  'latitude':-24.6272, #deg, MUST be negative (southern hemisphere) for now
                  'temperature':10+273.15, #K, temperature at Paranal
                  'humidity':14.5/100, #fraction, relative humidity at Paranal
                  'pressure':750, #mbar, pressure at Paranal
                  
                  'LR_VIS_major_axis':.702, #arcsec, major axis of observing modes apertures
                  'LR_VIS_hexagon_radius': 2, #radius of aperture hexagon array in hexagons
                  'LR_NIR_major_axis':.57, 
                  'LR_NIR_hexagon_radius': 2,
                  'HR_VIS_major_axis':.700, 
                  'HR_VIS_hexagon_radius': 3,
                  'HR_NIR_major_axis':.57, #This needs to be confirmed, but should be correct
                  'HR_NIR_hexagon_radius': 2,
                  
                  'custom_hexagon_radius': 0, #Change to be non-zero to overide aperture hexagon array
                           
                  'LR_VIS_B':[.390,.458], #um, MOSAIC bands
                  'LR_VIS_G':[.450,.591], 
                  'LR_VIS_R':[.586,.770], 
                  'LR_VIS_All':[.390,.770], 
                  'LR_NIR_IY':[.770,1.063], 
                  'LR_NIR_J':[1.01,1.395], 
                  'LR_NIR_H':[1.420,1.857], 
                  'LR_NIR_All':[.770,1.857], 
                  'HR_VIS_G':[.510,.568], #HR values need to be confirmed, especially VIS
                  'HR_VIS_R':[.610,.680],
                  'HR_NIR_IY':[.770,.907],
                  'HR_NIR_H':[1.523,1.658],

                  'sim_scale':.005, #arcsec/pixel, scale to carry out the simulation - smaller is slower, more accurate. Do not put above 0.01
                  'sim_HA_samples':21, #number of instantaneous snapshots to average over for the integration                 
                  'relative_plate_PA_angle':0, #deg, relative angle of the plate/apertures and PA=0. For PA=0 along semi major axis, set to 90 deg
                  }
    
        # for kwarg in kwargs.items():
        #     params[kwarg[0]]=kwarg[1]
        params.update(kwargs)
        self.config = params
        self.input={}
        self.output={}
        
    def load_MOSAIC_band(self,band,sampling):
        """
        Generates wavelengths for the simulation from the chosen band, 
        and acquires relevant diameter for the band (changes between LR/HR, and VIS/NIR)
        """
        self.input['sampling']=sampling
        self.input['band']=band
        
        wave_min,wave_max= self.config[band][0],self.config[band][1]

        self.output['wavelengths'] = np.arange(wave_min,wave_max,sampling)   
        self.output['major_axis']=self.config[band[0:6]+"_major_axis"]
        
    def load_HA(self,HA_start,HA_end,declination):
        """
        Calculates the airmasses to use in the simulation using the observation parameters (HA and dec)
        """
        HA_range=np.linspace(HA_start,HA_end,self.config['sim_HA_samples'])
        self.input['HA_range']=HA_range
        self.input['targ_dec']=declination

        #latitude needs to be negative for now
        lat = self.config['latitude'] * np.pi/180
        dec = declination*np.pi/180
        
        #Need to check if the target is below the horizon for the given list of HA, and if so remove it.
        LHA_below_horizon=np.rad2deg(np.arccos(-np.tan(lat)*np.tan(dec)))/15 
        if str(LHA_below_horizon) != 'nan': 
            for val in HA_range.copy(): 
                if abs(val) > abs(LHA_below_horizon):
                    print("At HA %2.2fh, target goes below horizon - removing this from HA range" % (val))
                    HA_range.remove(val)        
        if dec > np.pi/2 + lat: #If the target has a too high declination, it will never be seen at Cerro Paranal
            print("Target always below Horizon")
            return

        airmasses=1/(np.sin(lat)*np.sin(dec)+np.cos(lat)*np.cos(dec)*np.cos(Angle(HA_range*u.hour).rad))
        self.output['airmasses']=np.array(airmasses)
        
        para_angles=disp_utils.parallatic_angle(np.array(HA_range),dec,lat)
        self.output['raw_para_angles']=np.array(para_angles) #actual PAs
        
    def calculate_shifts(self, guide_waveref, aperture_waveref):
        """
        Calculates shifts of the wavelengths at the different airmasses using Fillipenko
        """  
        self.input['guide_waveref']=guide_waveref
        self.input['aperture_waveref']=aperture_waveref

        airmasses=self.output['airmasses']
        wavelengths=self.output['wavelengths']

        #centring refers to the index of the hour angles at which we centre the aperture/guiding on a wavelength
        centring_index=int((len(airmasses)-1)/2)

        centre_shift=disp_utils.diff_shift(aperture_waveref,airmasses[centring_index],guide_waveref,self.config) #shift of the original aperture centre wavelength from guide wavelength
        centring_q=self.output['raw_para_angles'][centring_index]

        raw_para_angles=self.output['raw_para_angles']
        para_angles=self.output['raw_para_angles'].copy()
        for i in range(0,len(para_angles)): #change in PAs from centring index
            para_angles[i]=para_angles[i]-self.output['raw_para_angles'][centring_index]

        shifts_para=[]
        phi=np.deg2rad(self.config['relative_plate_PA_angle'])
        for count,airmass in enumerate(airmasses): #for each airmass, calculate AD shift
            shift_vals=disp_utils.diff_shift(wavelengths,airmass,guide_waveref,self.config)  
            airmass_shifts=[]

            for i in range(0,len(shift_vals)):
                x=(shift_vals[i])*np.sin(raw_para_angles[count])-centre_shift*np.sin(centring_q)
                y=(shift_vals[i])*np.cos(raw_para_angles[count])-centre_shift*np.cos(centring_q)
                airmass_shifts.append([x*np.cos(phi)-y*np.sin(phi),y*np.cos(phi)+x*np.sin(phi)])
                
            shifts_para.append(airmass_shifts)

        self.output['shifts']=np.array(shifts_para)
        centre_shift_para=[-centre_shift*np.sin(centring_q),-centre_shift*np.cos(centring_q)]
        centre_shift_para=[centre_shift_para[0]*np.cos(phi)-centre_shift_para[1]*np.sin(phi),
                           centre_shift_para[1]*np.cos(phi)+centre_shift_para[0]*np.sin(phi)]
        self.output['centre_shift']=centre_shift_para
    
    def load_PSFs(self):
        """
        Generates the shifted moffat PSFs - one for for wavelength at each airmass
        """
        all_PSFs=[]
       
        for i in range(0,len(self.output['airmasses'])):
            PSFs=disp_utils.make_moffat_PSFs(self.output['wavelengths'],self.output['shifts'][i],
                                                         self.output['airmasses'][i],self.output['major_axis'],self.config)
            all_PSFs.append(PSFs)
            
        self.output['PSFs']=all_PSFs
        
    def load_aperture(self,band):
        """
        Generates apertures, one for each fibre in the bundle
        """
        apertures,apertures_table=disp_utils.make_aperture(band,self.output['major_axis'],self.config)
        self.output['apertures']=apertures
        self.apertures_table=apertures_table
           
    def calculate_transmissions(self):
        """
        Calculates the transmission of the PSFs through the apertures.
        """
        apertures_resized = np.repeat(self.output['apertures'][:,np.newaxis], np.shape(self.output['PSFs'])[1], axis=1)
        apertures_resized = np.repeat(apertures_resized[:,np.newaxis], np.shape(self.output['PSFs'])[0], axis=1)

        PSFs_through_apertures=apertures_resized*self.output['PSFs']
        self.output['PSFs_through_apertures']=PSFs_through_apertures
        
        transmissions=np.sum(PSFs_through_apertures,axis=(-1,-2))
        self.output['raw_transmissions']=transmissions
        
        integration_transmissions=np.mean(transmissions,axis=1)
        self.output['raw_integration_transmissions']=integration_transmissions #collapsed along airmass axis
        
    def plots(self):
        """
        Function to illustrate simulation results
        1) Transmission vs wavelength curves for individual fibres and entire bundle
        2) Track plot of monochromatic spot PSFs on the aperture over an integration
        """
        plt.style.use('bmh')

        fig,ax=plt.subplots(figsize=[7,5])
        plt.ylim(0,0.1)
        ax.set_ylabel("Individual Fibre Transmission")
        for fibre_trans in self.output['raw_integration_transmissions']:
            ax.plot(self.output['wavelengths'],fibre_trans,color='black',linewidth=0.5) 
        ax2=ax.twinx()
        ax2.axhline(0,label="Individual Fibre",color='black',linewidth=0.5)
        ax2.plot(self.output['wavelengths'],np.sum(self.output['raw_integration_transmissions'],axis=0),label="Bundle = {}um".format(self.input['aperture_waveref']),color='red')
        plt.axvline(self.input['guide_waveref'],label="Guide = {}um".format(self.input['guide_waveref']),color='black',linestyle='--',linewidth=0.5)
        ax2.set_ylabel("Bundle Transmission")
        ax2.set_ylim(0,1)
        ax.set_xlabel("Wavelength (um)")
        plt.legend()
        
        fig, ax = plt.subplots(figsize=[5,5]) 
        weights = np.linspace(0, len(self.output['wavelengths'])-1,4)
        norm = mpl.colors.Normalize(vmin=min(weights), vmax=max(weights))
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='seismic')
        circle1 = plt.Circle((0, 0), self.output['major_axis']/2, color='black', fill=False, label='~Aperture')
        ax.add_patch(circle1)    
        plt.axvline(0,color='black',linestyle='--',linewidth=0.7,label="PA = {}".format(self.config['relative_plate_PA_angle']))
        plt.scatter(self.output['centre_shift'][0],self.output['centre_shift'][1],label='Guide = {}um'.format(self.input['guide_waveref']),color='black',marker='+')
        plt.xlim(-0.4,0.4)
        plt.ylim(-0.4,0.4)
        shifts=self.output['shifts']
        for i in weights:
            xs,ys=[],[]
            for o in range(0,len(shifts)):
                xs.append(shifts[o][int(i)][0])
                ys.append(shifts[o][int(i)][1]) 
            plt.plot(xs,ys,marker='x',color=cmap.to_rgba(int(i)),label="{}um".format(round(self.output['wavelengths'][int(i)],4)))
        plt.legend()
        plt.xlabel("x (arcsec)")
        plt.ylabel("y (arcsec)")
            
    def run(self, HA_start, HA_end, declination, band,
            sampling=0.01, guide_waveref=-1, aperture_waveref=-1):
        """
        Function to run MOSAIC AD simulation for individual fibre transmission curves.
        Goes through multiple steps:
        1) Finds observation airmasses based on provided hour angles and declination
        2) Generates wavelengths to calculate transmission for based on the chosen MOSAIC observing band
        3) Calculates the shifts of these wavelengths at each airmass
        4) Generates aperture (for each fibre)
        5) Generates PSFs of the shifted moffats
        6) Calculates transmissions using apertures and PSFs
        
        Parameters
        ----------
        HA_start,HA_end: float, hours
            HA values to start and end the simulated observation at
            
        declination: float, degrees
            Declination of the target in the simulated observation
            
        band: string
            Which MOSAIC observing band to carry out the simulation in.
            Must be in the form "X_X_X", such as "LR_VIS_B" with the following options:
            {LR: {VIS: {B, G, R, All}, NIR: {IY, J, H, All}}, HR: {VIS: {G, R}, NIR: {IY, H}}} 
            
        sampling: float, um
            What interval in um to sample the wavelengths at, default 0.01um
            
        guide_waveref, aperture_waveref: float, um
            What wavelength to guide the telescope on and centre the apertures on halfway through the integration respectively
            defaults (when values set to -1) are wavelengths halfway through the chosen band 
        
        Returns
        -------
        wavelengths: 1D array of floats
            [um] wavelength samples the transmissions have been calculated for
            
        fibre_transmissions_dic: dictionary, containing 1D arrays of floats
            [0..1] Transmissions of each of the fibres for each different
            simulated wavelength. Dictionary key is the fibre id, {0,1...}.
            The 1D array axis is the wavelength
        """

        self.load_HA(HA_start,HA_end,declination)
        self.load_MOSAIC_band(band,sampling)
        
        if guide_waveref==-1 or aperture_waveref ==-1:
            band_mid=round((self.output['wavelengths'][0]+self.output['wavelengths'][-1])/2,3)
            guide_waveref,aperture_waveref=band_mid,band_mid
        self.input['guide_waveref']=guide_waveref
        self.input['aperture_waveref']=aperture_waveref
        
        self.calculate_shifts(guide_waveref,aperture_waveref)
        self.load_aperture(band)
        self.load_PSFs()
        self.calculate_transmissions()
        
        wavelengths=self.output['wavelengths']
        fibre_transmissions=self.output['raw_integration_transmissions'] 
        
        fibre_transmissions_dict={}
        for count,fibre_transmission in enumerate(fibre_transmissions):
            fibre_transmissions_dict[count]=fibre_transmission
        
        return wavelengths, fibre_transmissions_dict
