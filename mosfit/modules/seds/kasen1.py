"""Definitions for the `Kasen1` class."""
import numpy as np
from mosfit.modules.seds.sed import SED
import pickle 
from astropy import constants as c
from astropy import units as u
import os



# Important: Only define one ``Module`` class per file.


class Kasen1(SED):
    '''
    Defining the Kasen-simulation based SED

    FOR TYPE 1 == TIDAL TAIL
    #Frankencode
    Kamile Lukosiute Spooktober 2018
    '''

    # Kasen-calculated parameters
    MASS = np.array([.001, .0025, .005, .01, .02, .025, .03, .04, .05, .075, .1])
    VKIN = np.array([.03, .05, .1, .2, .3])
    XLAN = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-9])

    MASS_S = np.array(['0.001', '0.0025', '0.005', '0.010' ,'0.020', '0.025', '0.030', '0.040', '0.050', '0.075', '0.1' ])
    VKIN_S = np.array(['0.03', '0.05', '0.10', '0.20', '0.30'])
    XLAN_S = np.array(['1e-1', '1e-2', '1e-3', '1e-4', '1e-5', '1e-9'])

    C_CONST = c.c.cgs.value

    def __init__(self, **kwargs):
        super(Kasen1, self).__init__(**kwargs)

        # Read in times and frequencies arrays (same for all SEDs)
        self._dir_path = '/data/des51.b/data/kamile/'
        self._kasen_wavs = pickle.load( open(os.path.join(self._dir_path, 'kasen_seds/wavelength_angstroms.p'), "rb"))
        self._kasen_times = pickle.load( open(os.path.join(self._dir_path, 'kasen_seds/times_days.p'), "rb"))

        # create in memory the kasen_seds to later pick from 
        self._all_kasen_seds = []
        for m in self.MASS_S:
            for v in self.VKIN_S:
                for x in self.XLAN_S:
                            fname = 'kasen_seds/knova_d1_n10_m' + m + '_vk' + v + '_fd1.0_Xlan' + x + '.0.p'
                            if fname == 'kasen_seds/knova_d1_n10_m0.1_vk0.30_fd1.0_Xlan1e-1.0.p':
                                continue
                            kasen_sed = pickle.load( open(os.path.join(self._dir_path, fname) , "rb" ))
                            self._all_kasen_seds.append(kasen_sed)


    def weight(self, phi, theta):
        # viewing angle heta, half opening angle phi
        if phi + theta > np.pi/2.: 
            x = ( (np.sin(phi)**2. - np.cos(theta)**2.)**.5/
                np.sin(theta) )
        else:
            x = 0

        weight0 = ( (np.pi*(np.sin(self._phi)**2.)*(np.cos(self._theta)) + 
            2*(1-np.cos(self._theta))*(np.arcsin(x) - x*(1 - x**2.)**.5))/np.pi )

        # TYPE 1 == TIDAL TAIL SO GEOMETRIC FACTOR IS 1 - geometrical_weight
        return 1. -  weight0


    def process(self, **kwargs):
        kwargs = self.prepare_input(self.key('luminosities'), **kwargs)
        self._luminosities = kwargs[self.key('luminosities')]
        self._times = kwargs[self.key('rest_times')]
    
        self._band_indices = kwargs['all_band_indices']
        self._frequencies = kwargs['all_frequencies']

        # Physical parameters to index Kasen simulations
        self._vk = kwargs[self.key('vk')]
        self._xlan = kwargs[self.key('xlan')]
        self._mass = kwargs[self.key('Msph')]
        
        # Total weight function
        self._phi = kwargs[self.key('phi')] # half opening
        self._theta = kwargs[self.key('theta')] # viewing

        # mass fractional weight to convert between spherical and conical mass
        self._mass_weight =  ( 0.5 * ( (2+np.cos(self._phi) ) *
            ( 1 - np.cos(self._phi) )**2 +
            ( np.sin(self._phi)**2 ) * np.cos(self._phi) ) )


        self._weight_geom = self.weight(self._phi, self._theta)

        weight = 1./(1-self._mass_weight) * self._weight_geom
 
        # Some temp vars for speed.
        cc = self.C_CONST
        zp1 = 1.0 + kwargs[self.key('redshift')]
        czp1 = cc / zp1

        seds = []
        lums1 = []
        rest_wavs_dict = {}

        # Find nearest neighbors to the Kasen-calculated simulation
        m_closest = self.MASS_S[(np.abs(self.MASS-self._mass)).argmin()]
        v_closest = self.VKIN_S[(np.abs(self.VKIN-self._vk)).argmin()]
        x_closest = self.XLAN_S[(np.abs(self.XLAN-self._xlan)).argmin()]

        # if it's that pesky missing one
        if x_closest == '1e-1' and m_closest == '0.1' and v_closest == '0.30' :
            m_closest = '0.075'
            
        kasen_seds = {}
        # Find dictionary that corresponds to m_closest, v_closest, x_closest
        for i in self._all_kasen_seds:
            if i['Xlan'] == float(x_closest) and  i['mass'] == float(m_closest) and i['vk'] == float(v_closest):
                kasen_seds['SEDs'] = i['SEDs']
                break

        # For each time (luminosities as proxy)
        for li, lum in enumerate(self._luminosities):
            bi = self._band_indices[li]
            if bi >= 0:
                rest_wavs = rest_wavs_dict.setdefault(
                    bi, self._sample_wavelengths[bi] / zp1)
            else:
                rest_wavs = np.array(  # noqa: F841
                    [czp1 / self._frequencies[li]])

            # Find corresponding closest time
            t_closest_i = (np.abs(self._kasen_times-self._times[li])).argmin()

            # Evaluate the SED at the rest frame wavelengths 
            sed = np.array([])
            for w in rest_wavs:
                # find index of closest wav
                w_closest_i = np.abs(self._kasen_wavs-w).argmin()
                sed = np.append(sed, weight * kasen_seds['SEDs'][t_closest_i][w_closest_i] )


            # replace array w/very small val if t = 0 (hacky fix but whatver I have like two thesis weeks left)
            if self._times[li] == 0.0:
                sed[sed >= 0.] = 1.e30 # should be all values (can't have neg luminosity


	    # Calculate luminosity from sed
            L_t = np.trapz(weight * kasen_seds['SEDs'][t_closest_i], x=self._kasen_wavs)
            self._luminosities[li] = self._luminosities[li] +  L_t
            lums1.append(L_t)
            seds.append(sed)
            
        #This line turns all the nans to 0s
        seds[-1][np.isnan(seds[-1])] = 0.0
        

        seds = self.add_to_existing_seds(seds, **kwargs)
        return {'sample_wavelengths': self._sample_wavelengths, 'seds': seds, 'lum_1': lums1 }

